import numpy as np
from knoten import utils
from collections import namedtuple

Point = namedtuple("Point", 'x, y, z')
Sphere = namedtuple("Sphere", 'lat, lon, radius')

def test_sep_angle_right():
  pt1 = Point(1, 0, 0)
  pt2 = Point(0, 1, 0)
  np.testing.assert_array_equal(utils.sep_angle(pt1, pt2), np.pi / 2.0)

def test_sep_angle_acute():
  pt1 = Point(1, 0, 0)
  pt2 = Point(1, 1, 0)
  np.testing.assert_allclose(utils.sep_angle(pt1, pt2), np.pi / 4.0, atol=1e-12)

def test_sep_angle_obtuse():
  pt1 = Point(1, 0, 0)
  pt2 = Point(-1, 1, 0)
  np.testing.assert_array_equal(utils.sep_angle(pt1, pt2), 3.0 * np.pi / 4.0)

def test_sep_angle_normalization():
  pt1 = Point(1, 0, 0)
  pt2 = Point(1, 1, 0)
  pt3 = Point(100, 0, 0)
  pt4 = Point(100, 100, 0)
  np.testing.assert_array_equal(utils.sep_angle(pt1, pt2), utils.sep_angle(pt3, pt4))

def test_magnitude_unit():
  assert utils.magnitude(Point(1.0, 0.0, 0.0)) == 1.0
  assert utils.magnitude(Point(0.0, 1.0, 0.0)) == 1.0
  assert utils.magnitude(Point(0.0, 0.0, 1.0)) == 1.0

def test_magnitude_nonunit():
  assert utils.magnitude(Point(0.0, 0.0, 0.0)) == 0.0
  assert utils.magnitude(Point(2.0, 1.0, 4.0)) == np.sqrt(21.0)
  np.testing.assert_allclose(utils.magnitude(Point(0.2, 0.1, 0.4)), np.sqrt(0.21), atol=1e-12)

def test_distance():
  assert utils.distance(Point(1.0, 2.0, 3.0), Point(6.0, 5.0, 4.0)) == np.sqrt(35)

def test_spherical_to_rect():
  result = utils.spherical_to_rect(Sphere(0.0, 0.0, 1000.0))
  np.testing.assert_allclose(result.x, 1000.0, atol=1e-12)
  np.testing.assert_allclose(result.y, 0.0, atol=1e-12)
  np.testing.assert_allclose(result.z, 0.0, atol=1e-12)

  result = utils.spherical_to_rect(Sphere(0.0, np.pi, 1000.0))
  np.testing.assert_allclose( result.x, -1000.0, atol=1e-12)
  np.testing.assert_allclose( result.y, 0.0, atol=1e-12)
  np.testing.assert_allclose( result.z, 0.0, atol=1e-12)

  result = utils.spherical_to_rect(Sphere(np.pi / 2.0, 0.0, 1000.0))
  np.testing.assert_allclose( result.x, 0.0, atol=1e-12)
  np.testing.assert_allclose( result.y, 0.0, atol=1e-12)
  np.testing.assert_allclose( result.z, 1000.0, atol=1e-12)

  result = utils.spherical_to_rect(Sphere(np.pi / -2.0, 0.0, 1000.0))
  np.testing.assert_allclose( result.x, 0.0, atol=1e-12)
  np.testing.assert_allclose( result.y, 0.0, atol=1e-12)
  np.testing.assert_allclose( result.z, -1000.0, atol=1e-12)

def test_rect_to_spherical():
  result = utils.rect_to_spherical(Point(1000.0, 0.0, 0.0))
  np.testing.assert_array_equal(result, Sphere(0.0, 0.0, 1000.0))

  result = utils.rect_to_spherical(Point(-1000.0, 0.0, 0.0))
  np.testing.assert_array_equal(result, Sphere(0.0, np.pi, 1000.0))

  result = utils.rect_to_spherical(Point(0.0, 0.0, 1000.0))
  np.testing.assert_array_equal(result, Sphere(np.pi / 2.0, 0.0, 1000.0))

  result = utils.rect_to_spherical(Point(0.0, 0.0, -1000.0))
  np.testing.assert_array_equal(result,  Sphere(np.pi / -2.0, 0.0, 1000.0))

def test_ground_azimuth():
  LatLon = namedtuple("LatLon", "lat lon")

  ground_pt = LatLon(0, -180)
  subsolar_pt = LatLon(0, 90)
  np.testing.assert_array_equal(270.0, utils.ground_azimuth(ground_pt, subsolar_pt))

def test_perpendicular_vector():
  vec_a = Point(6.0, 6.0, 6.0)
  vec_b = Point(2.0, 0.0, 0.0)
  result = Point(0.0, 6.0, 6.0)
  np.testing.assert_array_equal(utils.perpendicular_vector(vec_a, vec_b), result)

def test_unit_vector():
  result = utils.unit_vector(Point(5.0, 12.0, 0.0))
  np.testing.assert_allclose(result[0], 0.384615, atol=1e-6)
  np.testing.assert_allclose(result[1], 0.923077, atol=1e-6)
  np.testing.assert_array_equal(result[2], 0.0)

def test_scale_vector():
  vec = Point(1.0, 2.0, -3.0)
  scalar = 3.0
  result = Point(3.0, 6.0, -9.0)
  np.testing.assert_array_equal(utils.scale_vector(vec, scalar), result)