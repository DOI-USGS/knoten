import unittest
from unittest import mock

from knoten import shape, utils

class TestEllipsoid(unittest.TestCase):

    def test_get_surface_normal(self):
        test_shape = shape.Ellipsoid(10, 10)

        ground_pt = utils.Point(1, 0, 0)
        normal = test_shape.get_surface_normal(ground_pt)
        
        self.assertEqual(normal, (1.0, 0.0, 0.0))

    def test_intersect_surface(self):
        test_shape = shape.Ellipsoid(10, 10)
        sensor_pos = utils.Point(100, 0, 0)
        look_vec = utils.Point(-1, 0, 0)

        intersection = test_shape.intersect_surface(sensor_pos, look_vec)

        self.assertEqual(intersection.x, 10.0)
        self.assertEqual(intersection.y, 0.0)
        self.assertEqual(intersection.z, 0.0)

    def tearDown(self):
        pass