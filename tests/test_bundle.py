from unittest import mock
import pytest

import numpy as np
import pandas as pd
from knoten import bundle
from csmapi import csmapi

@pytest.fixture
def control_network():
    df_dict = {
        'id': ['bob', 'bob', 'bob', 'tim', 'tim', 'sally', 'sally', 'sally', 'sally'],
        'serialnumber': ['a', 'b', 'c', 'd', 'b', 'a', 'b', 'c', 'd'],
        'line': np.arange(9),
        'sample': np.arange(9)[::-1],
        'aprioriX': np.zeros(9),
        'aprioriX': np.zeros(9),
        'aprioriY': np.zeros(9),
        'aprioriZ': np.zeros(9),
        'adjustedX': np.zeros(9),
        'adjustedY': np.zeros(9),
        'adjustedZ': np.zeros(9),
        'pointType': [2, 2, 2, 3, 3, 2, 2, 2, 2]
    }
    return pd.DataFrame(df_dict)

@pytest.fixture
def sensors():
    sensors = {
        'a': mock.MagicMock(spec=csmapi.RasterGM),
        'b': mock.MagicMock(spec=csmapi.RasterGM),
        'c': mock.MagicMock(spec=csmapi.RasterGM),
        'd': mock.MagicMock(spec=csmapi.RasterGM)
    }
    return sensors

def test_closest_approach_intersect():
    points = np.array([[-1, 1, 2], [0, 2, 2], [0, 1, 3]])
    directions = np.array([[1, 0, 0], [0, 2, 0], [0, 0, -1]])
    res = bundle.closest_approach(points, directions)

    np.testing.assert_allclose(res, [0, 1, 2])

def test_closest_approach_no_intersect():
    points = np.array([[-1, 1, 2], [0.5, 1-np.sqrt(3)/2.0, 2], [0.5, 1+np.sqrt(3)/2.0, 4]])
    directions = np.array([[0, 1, 0], [np.sqrt(3)/2.0, 0.5, 0], [0, 0, 1]])
    res = bundle.closest_approach(points, directions)

    np.testing.assert_allclose(res, [0, 1, 2], atol=1e-12)

def test_compute_ground_points(control_network, sensors):
    expected_bob = np.array([1.0, 7.0, 0.0])
    expected_sally = np.array([6.5, 1.5, 0.0])

    with mock.patch('knoten.bundle.closest_approach', side_effect=[expected_bob, expected_sally]) as mock_closest:
        out_df = bundle.compute_apriori_ground_points(control_network, sensors)
        mock_closest.assert_called()

    np.testing.assert_array_equal(
        out_df[out_df.id == "bob"][["aprioriX", "aprioriY", "aprioriZ"]].values,
        np.repeat(expected_bob[:, None], 3, axis=1).T)
    np.testing.assert_array_equal(
        out_df[out_df.id == "bob"][["adjustedX", "adjustedY", "adjustedZ"]].values,
        np.repeat(expected_bob[:, None], 3, axis=1).T)
    np.testing.assert_array_equal(
        out_df[out_df.id == "tim"][["aprioriX", "aprioriY", "aprioriZ"]].values,
        np.zeros((2, 3)))
    np.testing.assert_array_equal(
        out_df[out_df.id == "tim"][["adjustedX", "adjustedY", "adjustedZ"]].values,
        np.zeros((2, 3)))
    np.testing.assert_array_equal(
        out_df[out_df.id == "sally"][["aprioriX", "aprioriY", "aprioriZ"]].values,
        np.repeat(expected_sally[:, None], 4, axis=1).T)
    np.testing.assert_array_equal(
        out_df[out_df.id == "sally"][["adjustedX", "adjustedY", "adjustedZ"]].values,
        np.repeat(expected_sally[:, None], 4, axis=1).T)

def test_get_sensor_parameter():
    mock_sensor = mock.MagicMock(spec=csmapi.RasterGM)
    mock_sensor.getParameterSetIndices.return_value = [0, 1, 2]
    parameters = bundle.get_sensor_parameters(mock_sensor)

    mock_sensor.getParameterSetIndices.assert_called()
    assert len(parameters) == 3

def test_csm_parameters():
    mock_sensor = mock.MagicMock(spec=csmapi.RasterGM)
    test_parameter = bundle.CsmParameter(mock_sensor, 5)

    mock_sensor.getParameterName.assert_called_with(5)
    mock_sensor.getParameterType.assert_called_with(5)
    mock_sensor.getParameterUnits.assert_called_with(5)
    mock_sensor.getParameterValue.assert_called_with(5)
    assert test_parameter.index == 5
    assert test_parameter.name == mock_sensor.getParameterName.return_value
    assert test_parameter.type == mock_sensor.getParameterType.return_value
    assert test_parameter.units == mock_sensor.getParameterUnits.return_value
    assert test_parameter.value == mock_sensor.getParameterValue.return_value

def test_compute_sensor_partials():
    ground_pt = [9, 8, 10]
    sensor = mock.MagicMock(spec=csmapi.RasterGM)
    sensor.computeSensorPartials.side_effect = [(5, 3), (4, 2), (6, 1)]
    parameters = [mock.MagicMock(), mock.MagicMock(), mock.MagicMock()]
    partials = bundle.compute_sensor_partials(sensor, parameters, ground_pt)

    np.testing.assert_array_equal(partials, [(5, 4, 6), (3, 2, 1)])

def test_compute_ground_partials():
    ground_pt = [9, 8, 10]
    sensor = mock.MagicMock(spec=csmapi.RasterGM)
    sensor.computeGroundPartials.return_value = (1, 2, 3, 4, 5, 6)
    partials = bundle.compute_ground_partials(sensor, ground_pt)
    np.testing.assert_array_equal(partials, [[1, 2, 3], [4, 5, 6]])

def test_compute_jacobian(control_network, sensors):
    parameters = {sn: [mock.MagicMock()]*2 for sn in sensors}
    sensor_partials = [(i+1) * np.ones((2, 2)) for i in range(9)]
    ground_partials = [-(i+1) * np.ones((2, 3)) for i in range(9)]
    with mock.patch('knoten.bundle.compute_sensor_partials', side_effect=sensor_partials) as sensor_par_mock, \
         mock.patch('knoten.bundle.compute_ground_partials', side_effect=ground_partials) as ground_par_mock:
        J, coefficient_columns = bundle.compute_jacobian(control_network, sensors, parameters)

    expected_J = [
    [1, 1, 0, 0, 0, 0, 0, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0],
    [1, 1, 0, 0, 0, 0, 0, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0],
    [0, 0, 2, 2, 0, 0, 0, 0, -2, -2, -2, 0, 0, 0, 0, 0, 0],
    [0, 0, 2, 2, 0, 0, 0, 0, -2, -2, -2, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 3, 3, 0, 0, -3, -3, -3, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 3, 3, 0, 0, -3, -3, -3, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 4, 4, 0, 0, 0, -4, -4, -4, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 4, 4, 0, 0, 0, -4, -4, -4, 0, 0, 0],
    [0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 0, -5, -5, -5, 0, 0, 0],
    [0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 0, -5, -5, -5, 0, 0, 0],
    [6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, -6, -6],
    [6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, -6, -6],
    [0, 0, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -7, -7, -7],
    [0, 0, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -7, -7, -7],
    [0, 0, 0, 0, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, -8, -8, -8],
    [0, 0, 0, 0, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, -8, -8, -8],
    [0, 0, 0, 0, 0, 0, 9, 9, 0, 0, 0, 0, 0, 0, -9, -9, -9],
    [0, 0, 0, 0, 0, 0, 9, 9, 0, 0, 0, 0, 0, 0, -9, -9, -9]]
    np.testing.assert_array_equal(J, expected_J)

def test_compute_residuals(control_network, sensors):
    # sensor.groundToImage.side_effect = [csmapi.ImageCoord(-0.1, 7.8), csmapi.ImageCoord(0.7, 6.6), csmapi.ImageCoord(1.5, 5.4),
    #                                     csmapi.ImageCoord(2.3, 4.2), csmapi.ImageCoord(3.1, 4.9), csmapi.ImageCoord(5.8, 3.7),
    #                                     csmapi.ImageCoord(6.6, 2.5), csmapi.ImageCoord(7.4, 1.3), csmapi.ImageCoord(8.2, 0.1)]

    sensors['a'].groundToImage.side_effect = [csmapi.ImageCoord(-0.1, 7.8), csmapi.ImageCoord(5.8, 3.7)]
    sensors['b'].groundToImage.side_effect = [csmapi.ImageCoord(0.7, 6.6), csmapi.ImageCoord(3.1, 4.9), csmapi.ImageCoord(6.6, 2.5)]
    sensors['c'].groundToImage.side_effect = [csmapi.ImageCoord(1.5, 5.4), csmapi.ImageCoord(7.4, 1.3)]
    sensors['d'].groundToImage.side_effect = [csmapi.ImageCoord(2.3, 4.2), csmapi.ImageCoord(8.2, 0.1)]

    V = bundle.compute_residuals(control_network, sensors)
    assert V.shape == (18,)
    np.testing.assert_allclose(V, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1])
