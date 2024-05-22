from unittest import mock
import pytest

import numpy as np
import csmapi
from knoten import csm, sensor_utils, utils, bundle, shape, illuminator

@pytest.fixture
def mock_sensor():
    mock_sensor = mock.MagicMock(spec=csmapi.RasterGM)
    return mock_sensor

@pytest.fixture
def pt():
    return csmapi.ImageCoord(0.0, 0.0)

@pytest.fixture
def test_shape():
    return shape.Ellipsoid(0, 0)

@pytest.fixture
def test_illuminator():
    return illuminator.Illuminator(0, 0)


def test_phase_angle(mock_sensor, pt, test_shape, test_illuminator):
    with (
        mock.patch.object(illuminator.Illuminator, 'get_position_from_csm_sensor', return_value=utils.Point(100.0, 100.0, 0.0)),
        mock.patch.object(csm, 'get_state', return_value={'sensorPos': utils.Point(100.0, 0.0, 0.0), 'lookVec': utils.Point(-1.0, 0.0, 0.0)}),
        mock.patch.object(shape.Ellipsoid, 'intersect_surface', return_value=utils.Point(0,0,0))
    ):
        phase_angle = sensor_utils.phase_angle(pt, mock_sensor, test_shape, test_illuminator)

        np.testing.assert_allclose(phase_angle, 45.0)


def test_emission_angle(mock_sensor, pt, test_shape):
    with (
        mock.patch.object(shape.Ellipsoid, 'get_surface_normal', return_value=utils.Point(-1,0,0)), 
        mock.patch.object(shape.Ellipsoid, 'intersect_surface', return_value=utils.Point(0,0,0)),
        mock.patch.object(csm, 'get_state', return_value={'sensorPos': utils.Point(100.0, 0.0, 0.0), 'lookVec': utils.Point(-1.0, 0.0, 0.0)}),
        mock.patch.object(illuminator.Illuminator, 'get_position_from_csm_sensor', return_value=utils.Point(100.0, 100.0, 0.0))
    ):
        emission_angle = sensor_utils.emission_angle(pt, mock_sensor, test_shape)

        np.testing.assert_array_equal(emission_angle, 180.0)


def test_slant_distance(mock_sensor, pt, test_shape):
    with (
        mock.patch.object(csm, 'get_state', return_value={'sensorPos': utils.Point(100.0, 0.0, 0.0), 'lookVec': utils.Point(-1.0, 0.0, 0.0)}),
        mock.patch.object(shape.Ellipsoid, 'intersect_surface', return_value=utils.Point(0,0,0))
    ):
        slant_distance = sensor_utils.slant_distance(pt, mock_sensor, test_shape)

        np.testing.assert_array_equal(slant_distance, 100.0)


@mock.patch.object(csm, 'get_state', return_value={'sensorPos': utils.Point(-100.0, 0.0, 0.0)})
def test_target_center_distance(mock_sensor, pt):
    target_center_distance = sensor_utils.target_center_distance(pt, mock_sensor)

    np.testing.assert_array_equal(target_center_distance, 100.0)


@mock.patch.object(csm, 'get_state', return_value={'sensorPos': utils.Point(0.0, 0.0, 100.0)})
def test_sub_spacecraft_point(mock_sensor, pt):
    sub_spacecraft_point = sensor_utils.sub_spacecraft_point(pt, mock_sensor)

    np.testing.assert_array_equal(sub_spacecraft_point, [90.0, 0.0])


@mock.patch.object(csm, 'get_state', return_value={'sensorPos': utils.Point(100.0, 0.0, 0.0), 'lookVec': utils.Point(-1.0, 0.0, 0.0)})
@mock.patch.object(shape.Ellipsoid, 'intersect_surface', return_value=utils.Point(10.0, 0.0, 0.0))
def test_local_radius_intersection(mock_sensor, pt, test_shape):
    local_radius = sensor_utils.local_radius(pt, mock_sensor, test_shape)

    np.testing.assert_array_equal(local_radius, 10.0)


@mock.patch.object(csm, 'get_state', return_value={'sensorPos': utils.Point(1000.0, 0.0, 0.0), 'lookVec': utils.Point(-1000.0, 0.0, 0.0)})
@mock.patch.object(shape.Ellipsoid, 'intersect_surface', return_value=utils.Point(10.0, 0.0, 0.0))
def test_local_radius_ground(mock_sensor, pt, test_shape):
    local_radius = sensor_utils.local_radius(pt, mock_sensor, test_shape)

    np.testing.assert_array_equal(local_radius, 10.0)


@mock.patch.object(csm, 'get_state', return_value={'lookVec': utils.Point(-1.0, 0.0, 0.0)})
def test_right_ascension_declination(mock_sensor, pt):
    right_ascension_declination = sensor_utils.right_ascension_declination(pt, mock_sensor)

    np.testing.assert_array_equal(right_ascension_declination, [180.0, 0.0])


def test_line_resolution(mock_sensor, pt, test_shape):
    with (
        mock.patch.object(bundle, 'compute_image_partials', return_value=np.array([2, 1, 4, 4, 4, 8])),
        mock.patch.object(shape.Ellipsoid, 'intersect_surface', return_value=utils.Point(0,0,0)),
        mock.patch.object(csm, 'get_state', return_value={'sensorPos': utils.Point(100.0, 0.0, 0.0), 'lookVec': utils.Point(-1.0, 0.0, 0.0)})
    ):
        line_resolution = sensor_utils.line_resolution(pt, mock_sensor, test_shape)

        np.testing.assert_array_equal(line_resolution, 6.0)



def test_sample_resolution(mock_sensor, pt, test_shape):
    with (
        mock.patch.object(bundle, 'compute_image_partials', return_value=np.array([2, 1, 4, 4, 4, 8])),
        mock.patch.object(shape.Ellipsoid, 'intersect_surface', return_value=utils.Point(0,0,0)),
        mock.patch.object(csm, 'get_state', return_value={'sensorPos': utils.Point(100.0, 0.0, 0.0), 'lookVec': utils.Point(-1.0, 0.0, 0.0)})
    ):
        sample_resolution = sensor_utils.sample_resolution(pt, mock_sensor, test_shape)

        np.testing.assert_array_equal(sample_resolution, 9.0)


def test_pixel_resolution(mock_sensor, pt, test_shape):
    with (
        mock.patch.object(bundle, 'compute_image_partials', return_value=np.array([2, 1, 4, 4, 4, 8])),
        mock.patch.object(shape.Ellipsoid, 'intersect_surface', return_value=utils.Point(0,0,0)),
        mock.patch.object(csm, 'get_state', return_value={'sensorPos': utils.Point(100.0, 0.0, 0.0), 'lookVec': utils.Point(-1.0, 0.0, 0.0)})
    ):
        pixel_resolution = sensor_utils.pixel_resolution(pt, mock_sensor, test_shape)

        np.testing.assert_array_equal(pixel_resolution, 7.5)