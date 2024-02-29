from unittest import mock
import pytest

from plio.io.io_gdal import GeoDataset

import csmapi
from knoten import csm

@pytest.fixture
def mock_sensor():
    mock_sensor = mock.MagicMock(spec=csmapi.RasterGM)
    return mock_sensor

@pytest.fixture
def pt():
    return csmapi.ImageCoord(0.0, 0.0)

def test_generate_ground_point_with_float(mock_sensor):
    csm.generate_ground_point(0, (0.5, 0.5), mock_sensor)
    # The internal conversion from tuple to csmapi.ImageCoord means 
    # assert_called_once_with fails due to different addresses of
    # different objects.
    mock_sensor.imageToGround.assert_called_once()

def test_generate_ground_point_with_imagecoord(mock_sensor, pt):
    height = 0.0
    csm.generate_ground_point(height, pt, mock_sensor)
    mock_sensor.imageToGround.assert_called_once_with(pt, height)

@mock.patch.object(csm, 'get_radii', return_value=(10,10))
@mock.patch('pyproj.transformer.Transformer.transform', return_value=(0,0,0))
@mock.patch.object(csm, '_compute_intersection_distance', return_value=0)
def test_generate_ground_point_with_dtm(mock_sensor, pt, _):
    # Passing the mock_dem fixture fails for some reason. The
    # isinstance(obj, GeoDataset) check fails, causing the singldispath
    # to never dispatch to the func under test.
    mock_dem = mock.MagicMock(spec_set=GeoDataset)
    mock_dem.no_data_value = 10
    mock_dem.read_array.return_value = [[100]]
    mock_dem.latlon_to_pixel.return_value = (0.5,0.5)
    csm.generate_ground_point(mock_dem, pt, mock_sensor)
    # This call is mocked so that the intitial intersection and 
    # one iteration should occur. Therefore, the call count
    # should always be 2.
    assert mock_sensor.imageToGround.call_count == 2

from collections import namedtuple

@mock.patch.object(csm, 'get_radii', return_value=(10,10))
@mock.patch('pyproj.transformer.Transformer.transform', return_value=(0,0,0))
def test_generate_ground_point_with_dtm_ndv(mock_sensor, pt):
    # Passing the mock_dem fixture fails for some reason. The
    # isinstance(obj, GeoDataset) check fails, causing the singldispath
    # to never dispatch to the func under test.
    mock_dem = mock.MagicMock(spec_set=GeoDataset)

    # If the no data value equals the height, this should raise a value error
    mock_dem.no_data_value = 100
    mock_dem.read_array.return_value = [[100]]
    mock_dem.latlon_to_pixel.return_value = (0.5,0.5)
    with pytest.raises(ValueError):
        csm.generate_ground_point(mock_dem, pt, mock_sensor)

def test__compute_intersection_distance():
    Point = namedtuple("Point", 'x, y, z')
    pt1 = Point(0,0,0)
    pt2 = Point(1,1,1)
    dist = csm._compute_intersection_distance(pt1, pt2)
    assert dist == 1