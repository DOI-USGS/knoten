from knoten import csm, utils, bundle
from csmapi import csmapi
import numpy as np

def phase_angle(image_pt, sensor, shape, illuminator):
    """
    Computes and returns phase angle, in degrees for a given image point.
   
    Phase Angle: The angle between the vector from the intersection point to
    the observer (usually the spacecraft) and the vector from the intersection
    point to the illuminator (usually the sun).

    Parameters
    ----------
    image_pt : tuple
               Pair of x, y (sample, line) coordinates in pixel space

    sensor : object
             A CSM compliant sensor model object

    shape : object
            A shape model object

    illuminator: object
            An illuminator object

    Returns
    -------
     : np.ndarray
            phase angle in degrees
    """
    if not isinstance(image_pt, csmapi.ImageCoord):
        image_pt = csmapi.ImageCoord(*image_pt)

    sensor_state = csm.get_state(sensor, image_pt)
    ground_pt = shape.intersect_surface(sensor_state["sensorPos"], sensor_state["lookVec"])
  
    illum_pos = illuminator.get_position_from_csm_sensor(sensor, ground_pt)

    vec_a = utils.Point(sensor_state["sensorPos"].x - ground_pt.x, 
                        sensor_state["sensorPos"].y - ground_pt.y, 
                        sensor_state["sensorPos"].z - ground_pt.z)

    vec_b = utils.Point(illum_pos.x - ground_pt.x, 
                        illum_pos.y - ground_pt.y, 
                        illum_pos.z - ground_pt.z)
  
    return np.rad2deg(utils.sep_angle(vec_a, vec_b))

def emission_angle(image_pt, sensor, shape):
    """
    Computes and returns emission angle, in degrees, for a given image point.

    Emission Angle: The angle between the surface normal vector at the
    intersection point and the vector from the intersection point to the
    observer (usually the spacecraft). The emission angle varies from 0 degrees
    when the observer is viewing the sub-spacecraft point (nadir viewing) to 90
    degrees when the intercept is tangent to the surface of the target body.
    Thus, higher values of emission angle indicate more oblique viewing of the
    target.

    Parameters
    ----------
    image_pt : tuple
               Pair of x, y (sample, line) coordinates in pixel space

    sensor : object
             A CSM compliant sensor model object

    shape : object
            A shape model object

    Returns
    -------
     : np.ndarray
            emission angle in degrees
    """
    if not isinstance(image_pt, csmapi.ImageCoord):
        image_pt = csmapi.ImageCoord(*image_pt)

    sensor_state = csm.get_state(sensor, image_pt)
    ground_pt = shape.intersect_surface(sensor_state["sensorPos"], sensor_state["lookVec"])

    normal = shape.get_surface_normal(ground_pt)

    sensor_diff = utils.Point(sensor_state["sensorPos"].x - ground_pt.x, 
                              sensor_state["sensorPos"].y - ground_pt.y,
                              sensor_state["sensorPos"].z - ground_pt.z)
    
    return np.rad2deg(utils.sep_angle(normal, sensor_diff))

def slant_distance(image_pt, sensor, shape):
    """
    Computes the slant distance from the sensor to the ground point.

    Parameters
    ----------
    image_pt : tuple
               Pair of x, y (sample, line) coordinates in pixel space

    sensor : object
             A CSM compliant sensor model object

    shape : object
            A shape model object

    Returns
    -------
     : np.ndarray 
        slant distance in meters
    """
    if not isinstance(image_pt, csmapi.ImageCoord):
        image_pt = csmapi.ImageCoord(*image_pt)

    sensor_state = csm.get_state(sensor, image_pt)
    ground_pt = shape.intersect_surface(sensor_state["sensorPos"], sensor_state["lookVec"])

    return utils.distance(sensor_state["sensorPos"], ground_pt)

def target_center_distance(image_pt, sensor):
    """
    Calculates and returns the distance from the spacecraft to the target center.

    Parameters
    ----------
    image_pt : tuple
               Pair of x, y (sample, line) coordinates in pixel space

    sensor : object
             A CSM compliant sensor model object

    Returns
    -------
     : np.ndarray 
        target center distance in meters
    """
    if not isinstance(image_pt, csmapi.ImageCoord):
        image_pt = csmapi.ImageCoord(*image_pt)

    sensor_state = csm.get_state(sensor, image_pt)
    return utils.distance(sensor_state["sensorPos"], utils.Point(0,0,0))

def sub_spacecraft_point(image_pt, sensor):
    """
    Get the latitude and longitude of the sub-spacecraft point.

    Parameters
    ----------
    image_pt : tuple
               Pair of x, y (sample, line) coordinates in pixel space

    sensor : object
             A CSM compliant sensor model object

    Returns
    -------
     : np.ndarray
        sub spacecraft point in degrees
    """
    if not isinstance(image_pt, csmapi.ImageCoord):
        image_pt = csmapi.ImageCoord(*image_pt)

    sensor_state = csm.get_state(sensor, image_pt)
    lat_lon_rad = utils.rect_to_spherical(sensor_state["sensorPos"])

    return utils.radians_to_degrees(lat_lon_rad)

def local_radius(image_pt, sensor, shape):
    """
    Gets the local radius for a given image point.

    Parameters
    ----------
    image_pt : tuple
               Pair of x, y (sample, line) coordinates in pixel space

    sensor : object
             A CSM compliant sensor model object

    shape : object
            A shape model object

    Returns
    -------
     : np.ndarray
        local radius in meters
    """
    if not isinstance(image_pt, csmapi.ImageCoord):
        image_pt = csmapi.ImageCoord(*image_pt)

    sensor_state = csm.get_state(sensor, image_pt)
    ground_pt = shape.intersect_surface(sensor_state["sensorPos"], sensor_state["lookVec"])

    return utils.magnitude(ground_pt)

def right_ascension_declination(image_pt, sensor):
    """
    Computes the right ascension and declination for a given image point.

    Parameters
    ----------
    image_pt : tuple
               Pair of x, y (sample, line) coordinates in pixel space

    sensor : object
             A CSM compliant sensor model object

    Returns
    -------
     : tuple
       in the form (ra, dec) in degrees
    """
    if not isinstance(image_pt, csmapi.ImageCoord):
        image_pt = csmapi.ImageCoord(*image_pt)

    sensor_state = csm.get_state(sensor, image_pt)
    spherical_pt = utils.rect_to_spherical(sensor_state["lookVec"])
    
    ra_dec = utils.radians_to_degrees(spherical_pt)

    return (ra_dec.lon, ra_dec.lat)


def line_resolution(image_pt, sensor, shape):
    """
    Compute the line resolution in meters per pixel for the current set point.

    CSM sensor models do not expose all of the necessary parameters to do the
    same calculation as ISIS sensor models, so this uses a more time consuming but
    more accurate method and thus is equivalent to the oblique line resolution.

    For time dependent sensor models, this may also be the line-to-line resolution
    and not the resolution within a line or framelet. This is determined by the
    CSM model's ground computeGroundPartials method.

    Parameters
    ----------
    image_pt : tuple
               Pair of x, y (sample, line) coordinates in pixel space

    sensor : object
             A CSM compliant sensor model object

    shape : object
            A shape model object

    Returns
    -------
     : np.ndarray
        line resolution in meters/pixel
    """
    if not isinstance(image_pt, csmapi.ImageCoord):
        image_pt = csmapi.ImageCoord(*image_pt)

    sensor_state = csm.get_state(sensor, image_pt)
    ground_pt = shape.intersect_surface(sensor_state["sensorPos"], sensor_state["lookVec"])

    image_partials = bundle.compute_image_partials(sensor, ground_pt)

    return np.sqrt(image_partials[0] * image_partials[0] +
                   image_partials[2] * image_partials[2] +
                   image_partials[4] * image_partials[4])

def sample_resolution(image_pt, sensor, shape):
    """
    Compute the sample resolution in meters per pixel for the current set point.

    CSM sensor models do not expose all of the necessary parameters to do the
    same calculation as ISIS sensor models, so this uses a more time consuming but
    more accurate method and thus is equivalent to the oblique sample resolution.

    Parameters
    ----------
    image_pt : tuple
               Pair of x, y (sample, line) coordinates in pixel space

    sensor : object
             A CSM compliant sensor model object

    shape : object
            A shape model object

    Returns
    -------
     : np.ndarray
        sample resolution in meters/pixel
    """
    if not isinstance(image_pt, csmapi.ImageCoord):
        image_pt = csmapi.ImageCoord(*image_pt)

    sensor_state = csm.get_state(sensor, image_pt)
    ground_pt = shape.intersect_surface(sensor_state["sensorPos"], sensor_state["lookVec"])

    image_partials = bundle.compute_image_partials(sensor, ground_pt)

    return np.sqrt(image_partials[1] * image_partials[1] +
                   image_partials[3] * image_partials[3] +
                   image_partials[5] * image_partials[5])

def pixel_resolution(image_pt, sensor, shape):
    """
    Returns the pixel resolution at the current position in meters/pixel.

    Parameters
    ----------
    image_pt : tuple
               Pair of x, y (sample, line) coordinates in pixel space

    sensor : object
             A CSM compliant sensor model object

    shape : object
            A shape model object

    Returns
    -------
     : np.ndarray
        pixel resolution in meters/pixel
    """
    line_res = line_resolution(image_pt, sensor, shape)
    samp_res = sample_resolution(image_pt, sensor, shape)
    if (line_res < 0.0):
        return None
    if (samp_res < 0.0):
        return None
    return (line_res + samp_res) / 2.0


# def sub_solar_point(image_pt, sensor, illuminator, shape):
#     sensor_state = csm.get_state(sensor, image_pt)

#     illum_pos = illuminator.get_position_from_time(sensor_state)
#     lookVec = utils.Point(-illum_pos.x, -illum_pos.y, -illum_pos.z)

#     pt = shape.intersect_surface(illum_pos, lookVec)
#     return utils.Point(pt.x, pt.y, pt.z)


# def local_solar_time(image_pt, sensor, illuminator, shape):
#     if not isinstance(image_pt, csmapi.ImageCoord):
#         image_pt = csmapi.ImageCoord(*image_pt)

#     sensor_state = csm.get_state(sensor, image_pt)

#     sub_solar_pt = sub_solar_point(image_pt, sensor, illuminator, shape)
    
#     sub_solar_pt_degrees = utils.radians_to_degrees(utils.rect_to_spherical(sub_solar_pt))
  
#     ground_pt = shape.intersect_surface(sensor_state["sensorPos"], sensor_state["lookVec"])
#     spherical_pt = utils.rect_to_spherical(ground_pt)

#     spherical_pt_degrees = utils.radians_to_degrees(spherical_pt)

#     lst = spherical_pt_degrees.lon - sub_solar_pt_degrees.lon + 180.0
#     lst = lst / 15.0
#     if (lst < 0.0):
#         lst += 24.0
#     if (lst > 24.0):
#         lst -= 24.0
#     return lst

# def solar_longitude(image_pt, sensor, illuminator, body):
#     sensor_state = csm.get_state(sensor, image_pt)

#     illum_pos = illuminator.get_position_from_time(sensor_state)
#     illum_vel = illuminator.get_velocity(sensor_state)

#     body_rot = body.rotation(sensor_state)

#     sun_av = utils.unit_vector(utils.cross_product(illum_pos, illum_vel))
    
#     npole = [body_rot[6], body_rot[7], body_rot[8]]

#     z = sun_av
#     x = utils.unit_vector(utils.cross_product(npole, z))
#     y = utils.unit_vector(utils.cross_product(z, x))

#     trans = np.matrix([[x.x, x.y, x.z], [y.x, y.y, y.z], [z.x, z.y, z.z]])

#     pos = np.matmul(trans, illum_pos)
#     spherical_pos = utils.rect_to_spherical(pos)

#     longitude360 = np.rad2deg(spherical_pos.lon)

#     if (longitude360 != 360.0):
#       longitude360 -= 360 * np.floor(longitude360 / 360)

#     return longitude360