import datetime
import json
import os

from csmapi import csmapi
import jinja2
from osgeo import ogr
import numpy as np
import requests
import scipy.stats
from functools import singledispatch

from plio.io.io_gdal import GeoDataset

from knoten import utils
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, datetime.date):
            return obj.isoformat()
        return json.JSONEncoder.default(self, obj)

def get_radii(camera):
    """
    Given a sensor model, get the ellipsoid and return
    the semi_major and semi_minor radii.

    Parameters
    ----------
    camera : object
             A CSM compliant sensor model object

    Returns
    -------
     : tuple
       in the form (semi_major, semi_minor)
    """
    ellipsoid = csmapi.SettableEllipsoid.getEllipsoid(camera)
    semi_major = ellipsoid.getSemiMajorRadius()
    semi_minor = ellipsoid.getSemiMinorRadius()
    return semi_major, semi_minor

def create_camera(label, url='http://pfeffer.wr.usgs.gov/api/1.0/pds/'):
    """
    Given an ALE supported label, create a CSM compliant ISD file. This func
    defaults to supporting PDS labels. The URL kwarg can be used to point
    to the appropriate pfeffernusse endpoint for a given input data type.

    Parameters
    ----------
    label : str
            The image label

    url : str
          The service endpoint what is being acessed to generate an ISD.

    Returns
    -------
    model : object
            A CSM sensor model (or None if no associated model is available.)
    """
    data = json.dumps(label, cls=NumpyEncoder)
    response = requests.post(url, data=data)
    fname = ''
    with open(fname, 'w') as f:
        json.dump(response.json(), f)
    isd = csmapi.Isd(fname)

    plugin = csmapi.Plugin.findPlugin('UsgsAstroPluginCSM')

    model_name = response.json()['name_model']
    if plugin.canModelBeConstructedFromISD(isd, model_name):
        model = plugin.constructModelFromISD(isd, model_name)
        return model

def create_csm(image, verbose=False):
    """
    Given an image file, create a Community Sensor Model.

    Parameters
    ----------
    image : str
            The image filename to create a CSM for
    verbose : bool
              Print information about which plugins and models were attempted

    Returns
    -------
    model : object
            A CSM sensor model (or None if no associated model is available.)
    """
    isd = csmapi.Isd(image)
    plugins = csmapi.Plugin.getList()
    for plugin in plugins:
        if verbose:
            print(f'Trying plugin {plugin.getPluginName()}')
        num_models = plugin.getNumModels()
        for model_index in range(num_models):
            model_name = plugin.getModelName(model_index)
            warnings = csmapi.WarningList()
            if verbose:
                print(f'Trying model {plugin.getModelName(model_index)}')
            if plugin.canModelBeConstructedFromISD(isd, model_name, warnings):
                camera_warnings = csmapi.WarningList()
                camera = plugin.constructModelFromISD(isd, model_name, camera_warnings)
                if verbose:
                    for warning in camera_warnings:
                        print(f'Warning in function {warning.getFunction()}: "{warning.getMessage()}"')
                    print('Success!')
                return plugin.constructModelFromISD(isd, model_name)
            elif verbose:
                for warning in warnings:
                    print(f'Warning in function {warning.getFunction()}: "{warning.getMessage()}"')
                print('Failed!')

@singledispatch
def generate_ground_point(dem, image_pt, camera):
    '''
    Generates a longitude, latitude, and altitude coordinate for a given x, y
    pixel coordinate

    Parameters
    ----------
    dem : float or object
          Either a float that represents the height above the datum or a
          GeoDataset object generated from Plio off of a Digital Elevation
          Model (DEM)
    image_pt : tuple
               Pair of x, y coordinates in pixel space
    camera : object
             USGSCSM camera model object
    max_its : int, optional
              Maximum number of iterations to go through if the height does not
              converge
    tolerance : float, optional
                Number of decimal places to solve to when solving for height

    Returns
    -------
    : object
      CSM EcefCoord containing the newly computed lon, lat, and alt values
      corresponding to the original image_pt coordinates
    '''
    if not isinstance(image_pt, csmapi.ImageCoord):
        # Support a call where image_pt is in the form (x,y)
        image_pt = csmapi.ImageCoord(*image_pt)

    return camera.imageToGround(image_pt, dem)

@generate_ground_point.register(GeoDataset)
def _(dem, image_pt, camera, max_its = 20, tolerance = 0.001):
    if not isinstance(image_pt, csmapi.ImageCoord):
        # Support a call where image_pt is in the form (x,y)
        image_pt = csmapi.ImageCoord(*image_pt)

    intersection = generate_ground_point(0.0, image_pt, camera)
    iterations = 0
    semi_major, semi_minor = get_radii(camera)
    
    source_proj = f'+proj=cart +a={semi_major} +b={semi_minor}'
    dest_proj = f'+proj=lonlat +a={semi_major} +b={semi_minor}'
    transformer = utils.create_transformer(source_proj, dest_proj)

    while iterations != max_its:
        lon, lat, alt = transformer.transform(intersection.x, 
                                              intersection.y, 
                                              intersection.z,
                                              errcheck=True)

        px, py = dem.latlon_to_pixel(lat, lon)
        height = dem.read_array(1, [px, py, 1, 1])[0][0]
        if height == dem.no_data_value:
            raise ValueError(f'No DEM height at {lat}, {lon}')
    
        next_intersection = camera.imageToGround(image_pt, float(height))
        dist = max(abs(intersection.x - next_intersection.x),
            abs(intersection.y - next_intersection.y),
            abs(intersection.z - next_intersection.z))

        intersection = next_intersection
        iterations += 1
        if dist < tolerance:
            break
    return intersection

def generate_boundary(isize, npoints=10):
    '''
    Generates a bounding box given a camera model
    Parameters
    ----------
    isize : list
            image size in the form xsize, ysize
    npoints : int
              Number of points to generate between the corners of the bounding
              box per side.
    Returns
    -------
    boundary : list
               List of full bounding box
    '''
    x = np.linspace(0, isize[0], npoints)
    y = np.linspace(0, isize[1], npoints)
    boundary = [(i,0.) for i in x] + [(isize[0], i) for i in y[1:]] +\
               [(i, isize[1]) for i in x[::-1][1:]] + [(0.,i) for i in y[::-1][1:]]

    return boundary

def generate_latlon_boundary(camera, boundary, dem=0.0, radii=None, **kwargs):
    '''
    Generates a latlon bounding box given a camera model

    Parameters
    ----------
    camera : object
             csmapi generated camera model
    boundary : list
               List of boundary image coordinates
    dem : float or object
          Either a float that represents the height above the datum or a
          GeoDataset object generated from Plio off of a Digital Elevation
          Model (DEM)
    radii : tuple
            in the form (semimajor, semiminor) axes in meters. The default
            None, will attempt to get the radii from the camera object.

    Returns
    -------
    lons : list
           List of longitude values
    lats : list
           List of latitude values
    alts : list
           List of altitude values
    '''

    if radii is None:
        semi_major, semi_minor = get_radii(camera)
    else:
        semi_major, semi_minor = radii

    source_proj = f'+proj=cart +a={semi_major} +b={semi_minor}'
    dest_proj = f'+proj=lonlat +a={semi_major} +b={semi_minor}'
    transformer = utils.create_transformer(source_proj, dest_proj)

    gnds = np.empty((len(boundary), 3))

    for i, b in enumerate(boundary):
        # Could be potential errors or warnings from imageToGround
        try:
            gnd = generate_ground_point(dem, b, camera, **kwargs)
        except:
            pass

        gnds[i] = [gnd.x, gnd.y, gnd.z]

    lons, lats, alts = transformer.transform(gnds[:,0], gnds[:,1], gnds[:,2])
    return lons, lats, alts

def generate_gcps(camera, boundary, radii=None):
    '''
    Generates an area of ground control points formated as:
    <GCP Id="" Info="" Pixel="" Line="" X="" Y="" Z="" /> per record
    Parameters
    ----------
    camera : object
             csmapi generated camera model
    boundary : list
               List of image boundary coordinates
    radii : tuple
            in the form (semimajor, semiminor) axes in meters. The default
            None, will attempt to get the radii from the camera object.

    Returns
    -------
    gcps : list
           List of all gcp records generated
    '''
    if radii is None:
        semi_major, semi_minor = get_radii(camera)
    else:
        semi_major, semi_minor = radii

    lons, lats, alts = generate_latlon_boundary(camera, boundary, radii=radii)

    lla = np.vstack((lons, lats, alts)).T

    tr = zip(boundary, lla)

    gcps = []
    for i, t in enumerate(tr):
        l = '<GCP Id="{}" Info="{}" Pixel="{}" Line="{}" X="{}" Y="{}" Z="{}" />'.format(i, i, t[0][1], t[0][0], t[1][0], t[1][1], t[1][2])
        gcps.append(l)

    return gcps

def generate_latlon_footprint(camera, boundary, dem=0.0, radii=None, **kwargs):
    '''
    Generates a latlon footprint from a csmapi generated camera model
    Parameters
    ----------
    camera : object
             csmapi generated camera model
    boundary : list
               List of boundary image coordinates
    dem : float or object
          Either a float that represents the height above the datum or a
          GeoDataset object generated from Plio off of a Digital Elevation
          Model (DEM)
    radii : tuple
            in the form (semimajor, semiminor) axes in meters. The default
            None, will attempt to get the radii from the camera object.

    Returns
    -------
    : object
      ogr multipolygon containing between one and two polygons
    '''
    if radii is None:
        semi_major, semi_minor = get_radii(camera)
    else:
        semi_major, semi_minor = radii

    lons, lats, _ = generate_latlon_boundary(camera, boundary, dem=dem, radii=radii, **kwargs)

    # Transform coords from -180, 180 to 0, 360
    # Makes crossing the meridian easier to identify
    ll_coords = [*zip(((lons + 180) % 360), lats)]

    ring = ogr.Geometry(ogr.wkbLinearRing)
    wrap_ring = ogr.Geometry(ogr.wkbLinearRing)

    poly = ogr.Geometry(ogr.wkbPolygon)
    wrap_poly = ogr.Geometry(ogr.wkbPolygon)

    multipoly = ogr.Geometry(ogr.wkbMultiPolygon)

    current_ring = ring
    switch_point = (None, None)
    previous_point = ll_coords[0]

    for coord in ll_coords:

        coord_diff = previous_point[0] - coord[0]

        if coord_diff > 0 and np.isclose(previous_point[0], 360, rtol = 1e-03) and \
                              np.isclose(coord[0], 0, atol=1e0, rtol=1e-01):
            regression = scipy.stats.linregress([previous_point[0], coord[0]], [previous_point[1], coord[1]])
            slope, b = regression.slope, regression.intercept
            current_ring.AddPoint(360 - 180, (slope*360 + b))
            current_ring = wrap_ring
            switch_point = 0 - 180, (slope*0 + b)
            current_ring.AddPoint(*switch_point)

        elif coord_diff < 0 and np.isclose(previous_point[0], 0, atol=1e0, rtol=1e-01) and \
                                np.isclose(coord[0], 360, rtol = 1e-03):
            regression = scipy.stats.linregress([previous_point[0], coord[0]], [previous_point[1], coord[1]])
            slope, b = regression.slope, regression.intercept
            current_ring.AddPoint(0 - 180, (slope*0 + b))
            current_ring.AddPoint(*switch_point)
            current_ring = ring
            current_ring.AddPoint(360 - 180, (slope*360 + b))

        lat, lon = coord
        current_ring.AddPoint(lat - 180, lon)
        previous_point = coord

    poly.AddGeometry(ring)
    wrap_poly.AddGeometry(wrap_ring)

    if not wrap_poly.IsEmpty():
        multipoly.AddGeometry(wrap_poly)

    if not poly.IsEmpty():
        multipoly.AddGeometry(poly)

    return multipoly

def generate_bodyfixed_footprint(camera, boundary, radii=None):
    '''
    Generates a bodyfixed footprint from a csmapi generated camera model
    Parameters
    ----------
    camera : object
             csmapi generated camera model
    boundary : list
               List of boundary image coordinates
    radii : tuple
            in the form (semimajor, semiminor) axes in meters. The default
            None, will attempt to get the radii from the camera object.

    Returns
    -------
    : object
      ogr polygon
    '''
    if radii is None:
        semi_major, semi_minor = get_radii(camera)
    else:
        semi_major, semi_minor = radii

    latlon_fp = generate_latlon_footprint(camera, boundary, radii=radii)

    source_proj = f'+proj=lonlat +a={semi_major} +b={semi_minor}'
    dest_proj = f'+proj=cart +a={semi_major} +b={semi_minor}'
    transformer = utils.create_transformer(source_proj, dest_proj)
    # Step over all geometry objects in the latlon footprint
    for i in range(latlon_fp.GetGeometryCount()):
        latlon_coords = np.array(latlon_fp.GetGeometryRef(i).GetGeometryRef(0).GetPoints())

        # Check if the geometry object is populated with points
        if len(latlon_coords) > 0:
            x, y, z = transformer.transform(latlon_coords[:,0], latlon_coords[:,1], latlon_coords[:,2])

            # Step over all coordinate points in a geometry object and update said point
            for j, _ in enumerate(latlon_coords):
                latlon_fp.GetGeometryRef(i).GetGeometryRef(0).SetPoint(j, x[j], y[j], 0)

    return latlon_fp

def generate_vrt(raster_size, gcps, fpath,
                 no_data_value=0,
                 proj='+proj=longlat +a=3396190 +b=3376200 +no_defs'):
    """
    Create a GDAl VRT string from a list of ground control points and
    a projection.

    Parameters
    ----------
    raster_size : iterable
                  in the form xsize, ysize

    gcps : list
           of GCPs (likely created by the generate_gcp function)

    fpath : str
            path to the source file that the VRT points to

    no_data_value : numeric
                    the no data value for the VRT (default=0)

    proj : str
           A proj4 projection string for the VRT.

    Returns
    -------
    template : str
               The rendered VRT string

    Example
    --------
    >>> camera = create_camera(label)  # Get a camera
    >>> boundary = generate_boundary((100,100), 10)  # Compute the boundary points in image space
    >>> gcps = generate_gcps(camera, boundary)  # Generate GCPs using the boundary in lat/lon
    >>> vrt = generate_vrt((100,100), gcps, 'my_original_image.img')  # Create the vrt
    >>> # Then optionally, warp the VRT to render to disk
    >>> warp_options = gdal.WarpOptions(format='VRT', dstNodata=0)
    >>> gdal.Warp(some_output.tif, vrt, options=warp_options)
    """
    xsize, ysize = raster_size
    vrt = r'''<VRTDataset rasterXSize="{{ xsize }}" rasterYSize="{{ ysize }}">
     <Metadata/>
     <GCPList Projection="{{ proj }}">
     {% for gcp in gcps -%}
       {{gcp}}
     {% endfor -%}
    </GCPList>
     <VRTRasterBand dataType="Float32" band="1">
       <NoDataValue>{{ no_data_value }}</NoDataValue>
       <Metadata/>
       <ColorInterp>Gray</ColorInterp>
       <SimpleSource>
         <SourceFilename relativeToVRT="0">{{ fpath }}</SourceFilename>
         <SourceBand>1</SourceBand>
         <SourceProperties rasterXSize="{{ xsize }}" rasterYSize="{{ ysize }}"
    DataType="Float32" BlockXSize="512" BlockYSize="512"/>
         <SrcRect xOff="0" yOff="0" xSize="{{ xsize }}" ySize="{{ ysize }}"/>
         <DstRect xOff="0" yOff="0" xSize="{{ xsize }}" ySize="{{ ysize }}"/>
       </SimpleSource>
     </VRTRasterBand>
    </VRTDataset>'''

    context = {'xsize':xsize, 'ysize':ysize,
               'gcps':gcps,
               'proj':proj,
               'fpath':fpath,
               'no_data_value':no_data_value}
    template = jinja2.Template(vrt)
    return template.render(context)

def triangulate_ground_pt(cameras, image_pts):
    """
    Given a set of cameras and image points, find the ground point closest
    to the image rays for the image points.

    This function minimizes the sum of the squared distances from the ground
    point to the image rays.

    Parameters
    ----------
    cameras : list
              A list of CSM compliant sensor model objects
    image_pts : list
                A list of x, y image point tuples

    Returns
    -------
     : tuple
       The ground point as an (x, y, z) tuple
    """
    if len(cameras) != len(image_pts):
        raise ValueError("Lengths of cameras ({}) and image_pts ({}) must be the "
                         "same".format(len(cameras), len(image_pts)))

    M = np.zeros((3,3))
    b = np.zeros(3)
    unit_x = np.array([1, 0, 0])
    unit_y = np.array([0, 1, 0])
    unit_z = np.array([0, 0, 1])
    for camera, image_pt in zip(cameras, image_pts):
        if not isinstance(image_pt, csmapi.ImageCoord):
            image_pt = csmapi.ImageCoord(*image_pt)
        locus = camera.imageToRemoteImagingLocus(image_pt)
        look = np.array([locus.direction.x, locus.direction.y, locus.direction.z])
        pos = np.array([locus.point.x, locus.point.y, locus.point.z])
        look_squared = np.dot(look, look)
        M[0] += look[0] * look - look_squared * unit_x
        M[1] += look[1] * look - look_squared * unit_y
        M[2] += look[2] * look - look_squared * unit_z
        b += np.dot(pos, look) * look - look_squared * pos
    return tuple(np.dot(np.linalg.inv(M), b))
