import unittest
from unittest import mock

import csmapi

from knoten import utils
from knoten.illuminator import Illuminator

class TestIlluminator(unittest.TestCase):
    
    def test_get_position_from_csm_sensor(self):
        mock_sensor = mock.MagicMock(spec=csmapi.RasterGM)
        mock_sensor.getIlluminationDirection.return_value = utils.Point(1.0, 10.0, 0.0)

        test_illum = Illuminator()

        ground_pt = utils.Point(10.0, 10.0, 0.0)
        
        position = test_illum.get_position_from_csm_sensor(mock_sensor, ground_pt)
    
        self.assertEqual(position, (9.0, 0.0, 0.0))

    def tearDown(self):
        pass