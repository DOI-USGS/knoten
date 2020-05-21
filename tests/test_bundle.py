from unittest import mock
import pytest

import numpy as np
import pandas as pd
from knoten import bundle
from collections import OrderedDict
from csmapi import csmapi

@pytest.fixture
def control_network():
    df_dict = {
        'id': ['bob', 'bob', 'bob', 'tim', 'tim', 'sally', 'sally', 'sally', 'sally'],
        'serialnumber': ['a', 'b', 'c', 'd', 'b', 'a', 'b', 'c', 'd'],
        'line': np.arange(9),
        'sample': np.arange(9)[::-1],
        'aprioriCovar': [[], [], [], [], [], [], [], [], []],
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
    res, covar = bundle.closest_approach(points, directions)

    np.testing.assert_allclose(res, [0, 1, 2])

def test_closest_approach_no_intersect():
    points = np.array([[-1, 1, 2], [0.5, 1-np.sqrt(3)/2.0, 2], [0.5, 1+np.sqrt(3)/2.0, 4]])
    directions = np.array([[0, 1, 0], [np.sqrt(3)/2.0, 0.5, 0], [0, 0, 1]])
    res, covar = bundle.closest_approach(points, directions)

    np.testing.assert_allclose(res, [0, 1, 2], atol=1e-12)

def test_compute_ground_points(control_network, sensors):
    expected_bob = np.array([1.0, 7.0, 0.0])
    bob_covar = np.array([[1.0, 0.1, 0.2], [0.1, 1.5, 0.15], [0.2, 0.15, 3.0]])
    expected_sally = np.array([6.5, 1.5, 0.0])
    sally_covar = np.array([[2.0, 1.1, 0.6], [1.1, 1.0, 0.45], [0.6, 0.45, 3.2]])

    with mock.patch('knoten.bundle.closest_approach', side_effect=[(expected_bob, bob_covar), (expected_sally, sally_covar)]) as mock_closest:
        out_df = bundle.compute_apriori_ground_points(control_network, sensors)
        mock_closest.assert_called()

    for _, row in out_df[out_df.id == "bob"].iterrows():
        np.testing.assert_array_equal(row[["aprioriX", "aprioriY", "aprioriZ"]].values,
                                      expected_bob)
        np.testing.assert_array_equal(row[["adjustedX", "adjustedY", "adjustedZ"]].values,
                                      expected_bob)
        np.testing.assert_array_equal(list(row["aprioriCovar"]),
                                      [bob_covar[0,0], bob_covar[0,1], bob_covar[0,2],
                                       bob_covar[1,1], bob_covar[1,2], bob_covar[2,2]])
    for _, row in out_df[out_df.id == "tim"].iterrows():
        np.testing.assert_array_equal(row[["aprioriX", "aprioriY", "aprioriZ"]].values,
                                      np.zeros(3))
        np.testing.assert_array_equal(row[["adjustedX", "adjustedY", "adjustedZ"]].values,
                                      np.zeros(3))
        assert not list(row["aprioriCovar"])
    for _, row in out_df[out_df.id == "sally"].iterrows():
        np.testing.assert_array_equal(row[["aprioriX", "aprioriY", "aprioriZ"]].values,
                                      expected_sally)
        np.testing.assert_array_equal(row[["adjustedX", "adjustedY", "adjustedZ"]].values,
                                      expected_sally)
        np.testing.assert_array_equal(list(row["aprioriCovar"]),
                                      [sally_covar[0,0], sally_covar[0,1], sally_covar[0,2],
                                       sally_covar[1,1], sally_covar[1,2], sally_covar[2,2]])

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
    coefficient_columns = OrderedDict()
    coefficient_columns['a'] = (0, 2)
    coefficient_columns['b'] = (2, 4)
    coefficient_columns['c'] = (4, 6)
    coefficient_columns['d'] = (6, 8)
    coefficient_columns['bob'] = (8, 11)
    coefficient_columns['tim'] = (11, 14)
    coefficient_columns['sally'] = (14, 17)
    with mock.patch('knoten.bundle.compute_sensor_partials', side_effect=sensor_partials) as sensor_par_mock, \
         mock.patch('knoten.bundle.compute_ground_partials', side_effect=ground_partials) as ground_par_mock:
        J = bundle.compute_jacobian(control_network, sensors, parameters, coefficient_columns)

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

def test_compute_sigma0():
    V = np.arange(0, 16) + 1
    W_obs = np.diag(np.arange(16, 0, -1))
    W_params = np.array(
        [[1,  2,  3,  0,  0,  0],
         [4,  5,  6,  0,  0,  0],
         [7,  8,  9,  0,  0,  0],
         [0,  0,  0, -1, -2, -3],
         [0,  0,  0, -4, -5, -6],
         [0,  0,  0, -7, -8, -9]]
     )
    dX = np.arange(-6, 0)
    assert bundle.compute_sigma(V, dX, W_params, W_obs) == np.sqrt(7809 / 10)

def test_compute_sigma0_sparse():
    V = np.arange(0, 16) + 1
    W_obs = np.diag(np.arange(16, 0, -1))
    W_sensors  = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    W_points = {
        "point_1" : np.array([[-1, -2, -3], [-4, -5, -6], [-7, -8, -9]])
    }
    dX = np.arange(-6, 0)
    column_dict = {
        "image_1" : (0, 3),
        "point_1" : (3, 6)
    }
    assert bundle.compute_sigma_sparse(V, dX, W_sensors, W_points, W_obs, column_dict) == np.sqrt(7809 / 10)


def test_compute_image_weight():
    def mock_covar(index_1, index_2):
        if index_1 != index_2:
            return 0
        return index_1 + 1

    sensor = mock.MagicMock(spec=csmapi.RasterGM)
    sensor.getParameterCovariance = mock_covar
    params = [mock.MagicMock(), mock.MagicMock(), mock.MagicMock()]
    params[0].index = 0
    params[1].index = 1
    params[2].index = 3

    np.testing.assert_allclose(
        bundle.compute_image_weight(sensor, params),
        [[1, 0,   0],
         [0, 1/2, 0],
         [0, 0,   1/4]])

def test_compute_point_weight(control_network):
    control_network.at[(control_network['id'] == 'bob').idxmax(), 'aprioriCovar'] = np.array([1, 0, 0, 2, 0, 3])
    np.testing.assert_allclose(
        bundle.compute_point_weight(control_network, 'bob'),
        [[1, 0,   0],
         [0, 1/2, 0],
         [0, 0,   1/3]])

def test_update_parameters(control_network, sensors):
    def mock_getParameterValue(index):
        return index + 1

    params = {}
    coefficient_columns = {}
    current_column = 0

    for sn, sensor in sensors.items():
        sensor.getParameterValue = mock_getParameterValue
        sensor_params = []
        for param_idx in range(6):
            param_mock = mock.MagicMock()
            param_mock.index = param_idx
            sensor_params.append(param_mock)
        params[sn] = sensor_params
        coefficient_columns[sn] = (current_column, current_column+6)
        current_column += 6

    for point_id in control_network['id'].unique():
        coefficient_columns[point_id] = (current_column, current_column+3)
        current_column += 3

    updates = np.arange(0, current_column)

    bundle.update_parameters(sensors, params, control_network, updates, coefficient_columns)

    for sn, sensor in sensors.items():
        coeff_start = coefficient_columns[sn][0]
        for param_idx in range(6):
            sensor.setParameterValue.assert_any_call(
                param_idx,
                param_idx + 1 + coeff_start + param_idx
            )

    for point_id, group in control_network.groupby('id'):
        coeff_start = coefficient_columns[point_id][0]
        coeff_end = coefficient_columns[point_id][1]
        num_measures = len(group)
        np.testing.assert_allclose(
            group[['adjustedX', 'adjustedY', 'adjustedZ']].values,
            np.tile(updates[coeff_start:coeff_end], num_measures).reshape((num_measures,3))
        )
