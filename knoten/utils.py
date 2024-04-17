import pyproj
import numpy as np 
from collections import namedtuple

def sep_angle(a_pt, b_pt, c_pt):
    return sep_angle(a_pt - b_pt, c_pt - b_pt)

def sep_angle(a_vec, b_vec):
    dot_prod = a_vec.x * b_vec.x + a_vec.y * b_vec.y + a_vec.z * b_vec.z
    dot_prod /= magnitude(a_vec) * magnitude(b_vec)

    if(dot_prod >= 1.0): return 0.0
    if(dot_prod <= -1.0): return np.pi

    return np.arccos(dot_prod)

def magnitude(vec):
    return np.sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z)

def distance(start, stop):
    Point = namedtuple("Point", 'x, y, z')
    diff = Point(stop.x - start.x, stop.y - start.y, stop.z - start.z)

    return magnitude(diff)

def radiansToDegrees(radian_lat_lon):
    LatLon = namedtuple("LatLon", 'lat lon')

    degree_lon = radian_lat_lon.lon
    if (degree_lon < 0):
      degree_lon += 2 * np.pi

    degree_lon = np.rad2deg(degree_lon)
    degreeLat = np.rad2deg(radian_lat_lon.lat)
    return LatLon(degreeLat, degree_lon)

def spherical_to_rect(spherical):
    Point = namedtuple("Point", 'x, y, z')

    x = spherical.radius * np.cos(spherical.lat) * np.cos(spherical.lon)
    y = spherical.radius * np.cos(spherical.lat) * np.sin(spherical.lon)
    z = spherical.radius * np.sin(spherical.lat)

    return Point(x, y, z)

def rect_to_spherical(rectangular):
    Sphere = namedtuple("Sphere", 'lat, lon, radius')

    rad = magnitude(rectangular)
    if (rad < 1e-15):
      return Sphere(0.0, 0.0, 0.0)
    
    return Sphere(
      np.arcsin(rectangular.z / rad),
      np.arctan2(rectangular.y, rectangular.x),
      rad
    )

def ground_azimuth(ground_pt, sub_pt):
    LatLon = namedtuple("LatLon", 'lat lon')

    if (ground_pt.lat >= 0.0):
      a = (90.0 - sub_pt.lat) * np.pi / 180.0
      b = (90.0 - ground_pt.lat) * np.pi / 180.0
    else:
      a = (90.0 + sub_pt.lat) * np.pi / 180.0
      b = (90.0 + ground_pt.lat) * np.pi / 180.0

    cs = LatLon(0, sub_pt.lon)
    cg = LatLon(0, ground_pt.lon)

    if (cs.lon > cg.lon):
      if ((cs.lon - cg.lon) > 180.0):
        while ((cs.lon - cg.lon) > 180.0): 
           cs = LatLon(0, cs.lon - 360.0)
    if (cg.lon > cs.lon):
      if ((cg.lon-cs.lon) > 180.0):
        while ((cg.lon-cs.lon) > 180.0): 
           cg = LatLon(0, cg.lon - 360.0)
    
    if (sub_pt.lat > ground_pt.lat):
      if (cs.lon < cg.lon):
        quad = 2
      else:
        quad = 1
    elif (sub_pt.lat < ground_pt.lat):
      if (cs.lon < cg.lon):
        quad = 3
      else:
        quad = 4
    else:
      if (cs.lon > cg.lon):
        quad = 1
      elif (cs.lon < cg.lon):
        quad = 2
      else:
        return 0.0

    C = (cg.lon - cs.lon) * np.pi / 180.0
    if (C < 0):
       C = -C

    c = np.arccos(np.cos(a) * np.cos(b) + np.sin(a) * np.sin(b) * np.cos(C))

    azimuth = 0.0

    if (np.sin(b) == 0.0 or np.sin(c) == 0.0):
      return azimuth

    intermediate = (np.cos(a) - np.cos(b) * np.cos(c)) / (np.sin(b) * np.sin(c))
    if (intermediate < -1.0):
      intermediate = -1.0
    elif (intermediate > 1.0):
      intermediate = 1.0

    A = np.arccos(intermediate) * 180.0 / np.pi

    if (ground_pt.lat >= 0.0):
        if (quad == 1 or quad == 4):
            azimuth = A
        elif (quad == 2 or quad == 3):
            azimuth = 360.0 - A
    else: 
      if (quad == 1 or quad == 4):
        azimuth = 180.0 - A
      elif (quad == 2 or quad == 3):
        azimuth = 180.0 + A
    return azimuth

def crossProduct(a_vec, b_vec):
    Point = namedtuple("Point", 'x, y, z')

    x = a_vec.y * b_vec.z - a_vec.z * b_vec.y
    y = a_vec.z * b_vec.x - a_vec.x * b_vec.z
    z = a_vec.x * b_vec.y - a_vec.y * b_vec.x
    return Point(x, y, z)


def unit_vector(vec):
    mag = magnitude(vec)
    return vec / mag

def perpendicular_vector(a_vec, b_vec):
    if (magnitude(a_vec) == 0):
      return b_vec

    a_norm = unit_vector(a_vec)
    b_norm = unit_vector(b_vec)

    angle = a_norm * b_norm
    a_mag = magnitude(a_vec)
    mag_p = angle * a_mag

    p = b_norm * mag_p
    q = a_vec - p

    return q

def scale_vector(vec, scalar):
    Point = namedtuple("Point", 'x, y, z')

    return Point(vec.x * scalar, vec.y * scalar, vec.z * scalar)

def matrixVecProduct(mat, vec):
    Point = namedtuple("Point", 'x, y, z')

    x = mat.a.x * vec.x + mat.a.y * vec.y + mat.a.z * vec.z
    y = mat.b.x * vec.x + mat.b.y * vec.y + mat.b.z * vec.z
    z = mat.c.x * vec.x + mat.c.y * vec.y + mat.c.z * vec.z

    return Point(x, y, z)

def reproject(record, semi_major, semi_minor, source_proj, dest_proj, **kwargs):
    """
    Thin wrapper around PyProj's Transform() function to transform 1 or more three-dimensional
    point from one coordinate system to another. If converting between Cartesian
    body-centered body-fixed (BCBF) coordinates and Longitude/Latitude/Altitude coordinates,
    the values input for semi-major and semi-minor axes determine whether latitudes are
    planetographic or planetocentric and determine the shape of the datum for altitudes.
    If semi_major == semi_minor, then latitudes are interpreted/created as planetocentric
    and altitudes are interpreted/created as referenced to a spherical datum.
    If semi_major != semi_minor, then latitudes are interpreted/created as planetographic
    and altitudes are interpreted/created as referenced to an ellipsoidal datum.

    Parameters
    ----------
    record : object
          Pandas series object

    semi_major : float
              Radius from the center of the body to the equater

    semi_minor : float
              Radius from the pole to the center of mass

    source_proj : str
                      Pyproj string that defines a projection space ie. 'geocent'

    dest_proj : str
                   Pyproj string that defines a project space ie. 'latlon'

    Returns
    -------
    : list
    Transformed coordinates as y, x, z

    """
    transformer = pyproj.Transformer.from_crs(f'+proj={source_proj} +a={semi_major} +b={semi_minor}',
                                              f'+proj={dest_proj} +a={semi_major} +b={semi_minor}',
                                              always_xy=True)
    source_proj = f'+proj={source_proj} +a={semi_major} +b={semi_minor}'
    dest_proj = f'+proj={dest_proj} +a={semi_major} +b={semi_minor}'
    transformer = create_transformer(source_proj, dest_proj)
    x, y, z = transformer.transform(record[0], record[1], record[2], errcheck=True)

    return y, x, z

def create_transformer(source_proj, dest_proj):
    return pyproj.Transformer.from_crs(source_proj,
                                       dest_proj,
                                       always_xy=True)