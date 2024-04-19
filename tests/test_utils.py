import numpy as np
from knoten import utils

def test_sep_angle_right():
  pt1 = utils.Point(1, 0, 0)
  pt2 = utils.Point(0, 1, 0)
  np.testing.assert_array_equal(utils.sep_angle(pt1, pt2), np.pi / 2.0)

def test_sep_angle_acute():
  pt1 = utils.Point(1, 0, 0)
  pt2 = utils.Point(1, 1, 0)
  np.testing.assert_allclose(utils.sep_angle(pt1, pt2), np.pi / 4.0, atol=1e-12)

def test_sep_angle_obtuse():
  pt1 = utils.Point(1, 0, 0)
  pt2 = utils.Point(-1, 1, 0)
  np.testing.assert_array_equal(utils.sep_angle(pt1, pt2), 3.0 * np.pi / 4.0)

def test_sep_angle_normalization():
  pt1 = utils.Point(1, 0, 0)
  pt2 = utils.Point(1, 1, 0)
  pt3 = utils.Point(100, 0, 0)
  pt4 = utils.Point(100, 100, 0)
  np.testing.assert_array_equal(utils.sep_angle(pt1, pt2), utils.sep_angle(pt3, pt4))

def test_magnitude_unit():
  assert utils.magnitude(utils.Point(1.0, 0.0, 0.0)) == 1.0
  assert utils.magnitude(utils.Point(0.0, 1.0, 0.0)) == 1.0
  assert utils.magnitude(utils.Point(0.0, 0.0, 1.0)) == 1.0

def test_magnitude_nonunit():
  assert utils.magnitude(utils.Point(0.0, 0.0, 0.0)) == 0.0
  assert utils.magnitude(utils.Point(2.0, 1.0, 4.0)) == np.sqrt(21.0)
  np.testing.assert_allclose(utils.magnitude(utils.Point(0.2, 0.1, 0.4)), np.sqrt(0.21), atol=1e-12)

def test_distance():
  assert utils.distance(utils.Point(1.0, 2.0, 3.0), utils.Point(6.0, 5.0, 4.0)) == np.sqrt(35)

def test_spherical_to_rect():
  result = utils.spherical_to_rect(utils.Sphere(0.0, 0.0, 1000.0))
  np.testing.assert_allclose(result.x, 1000.0, atol=1e-12)
  np.testing.assert_allclose(result.y, 0.0, atol=1e-12)
  np.testing.assert_allclose(result.z, 0.0, atol=1e-12)

  result = utils.spherical_to_rect(utils.Sphere(0.0, np.pi, 1000.0))
  np.testing.assert_allclose( result.x, -1000.0, atol=1e-12)
  np.testing.assert_allclose( result.y, 0.0, atol=1e-12)
  np.testing.assert_allclose( result.z, 0.0, atol=1e-12)

  result = utils.spherical_to_rect(utils.Sphere(np.pi / 2.0, 0.0, 1000.0))
  np.testing.assert_allclose( result.x, 0.0, atol=1e-12)
  np.testing.assert_allclose( result.y, 0.0, atol=1e-12)
  np.testing.assert_allclose( result.z, 1000.0, atol=1e-12)

  result = utils.spherical_to_rect(utils.Sphere(np.pi / -2.0, 0.0, 1000.0))
  np.testing.assert_allclose( result.x, 0.0, atol=1e-12)
  np.testing.assert_allclose( result.y, 0.0, atol=1e-12)
  np.testing.assert_allclose( result.z, -1000.0, atol=1e-12)

def test_rect_to_spherical():
  result = utils.rect_to_spherical(utils.Point(1000.0, 0.0, 0.0))
  np.testing.assert_array_equal(result, utils.Sphere(0.0, 0.0, 1000.0))

  result = utils.rect_to_spherical(utils.Point(-1000.0, 0.0, 0.0))
  np.testing.assert_array_equal(result, utils.Sphere(0.0, np.pi, 1000.0))

  result = utils.rect_to_spherical(utils.Point(0.0, 0.0, 1000.0))
  np.testing.assert_array_equal(result, utils.Sphere(np.pi / 2.0, 0.0, 1000.0))

  result = utils.rect_to_spherical(utils.Point(0.0, 0.0, -1000.0))
  np.testing.assert_array_equal(result,  utils.Sphere(np.pi / -2.0, 0.0, 1000.0))

def test_ground_azimuth():
  ground_pt = utils.LatLon(0, -180)
  subsolar_pt = utils.LatLon(0, 90)
  np.testing.assert_array_equal(270.0, utils.ground_azimuth(ground_pt, subsolar_pt))

def test_perpendicular_vector():
  vec_a = utils.Point(6.0, 6.0, 6.0)
  vec_b = utils.Point(2.0, 0.0, 0.0)
  result = utils.Point(0.0, 6.0, 6.0)
  np.testing.assert_array_equal(utils.perpendicular_vector(vec_a, vec_b), result)

def test_unit_vector():
  result = utils.unit_vector(utils.Point(5.0, 12.0, 0.0))
  np.testing.assert_allclose(result[0], 0.384615, atol=1e-6)
  np.testing.assert_allclose(result[1], 0.923077, atol=1e-6)
  np.testing.assert_array_equal(result[2], 0.0)

def test_scale_vector():
  vec = utils.Point(1.0, 2.0, -3.0)
  scalar = 3.0
  result = utils.Point(3.0, 6.0, -9.0)
  np.testing.assert_array_equal(utils.scale_vector(vec, scalar), result)

def test_matrix_vec_product():
  vec_a = utils.Point(0.0,  1.0,  0.0)
  vec_b = utils.Point(-1.0,  0.0,  0.0)
  vec_c = utils.Point(0.0,  0.0,  1.0)
  mat = utils.Matrix(vec_a, vec_b, vec_c)
  vec = utils.Point(1.0,  2.0,  3.0)

  result = utils.Point(2.0, -1.0, 3.0)
  np.testing.assert_array_equal(result, utils.matrix_vec_product(mat, vec))