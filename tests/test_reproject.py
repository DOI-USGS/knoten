from unittest import mock
import pytest

from autocnet.transformation import spatial

def test_reproject():
    with mock.patch('pyproj.transform', return_value=[1,1,1]) as mock_pyproj:
        res = spatial.reproject([1,1,1], 10, 10, 'geocent', 'latlon')
        mock_pyproj.assert_called_once()
        assert res == (1,1,1)
