import datetime
from functools import singledispatch, update_wrapper
import json
import os

from csmapi import csmapi

from gdal import ogr
import numpy as np
from plio.io.io_gdal import GeoDataSet
import pvl
import pyproj
import requests
import scipy.stats

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, datetime.date):
            return obj.isoformat()
        return json.JSONEncoder.default(self, obj)


def methodispatch(func):
    dispatcher = singledispatch(func)
    def wrapper(*args, **kw):
        return dispatcher.dispatch(args[1].__class__)(*args, **kw)
    wrapper.register = dispatcher.register
    update_wrapper(wrapper, func)
    return wrapper


class Sensor(object):
    def __init__(self, image):
        isd = csmapi.Isd(image)
        plugins = csmapi.Plugin.getList()
        for plugin in plugins:
            num_models = plugin.getNumModels()
            for model_index in range(num_models):
                model_name = plugin.getModelName(model_index)
                if plugin.canModelBeConstructedFromISD(isd, model_name):
                    self.model = plugin.constructModelFromISD(isd, model_name)

    @classmethod
    def request_isd(cls, image, url='http://pfeffer.wr.usgs.gov/api/1.0/pds/')
        """
        Instantiate a CSM sensor model by first requesting an ISD from a remote
        service and writing the ISD to a json file adjacent to the passed image.

        The image is opened and the method attempts to parse the label using the
        PVL library. The label is then encoded and passed to the service.

        Parameters
        ----------
        image : str
                path to the image file

        url : str
              The URL of thge remote service used to generate an ISD.
        """
        label = pvl.dumps(image).decode()
        image = json.dumps(label, cls=NumpyEncoder)
        response = requests.post(url, data=image)
        fname = os.path.splitext(image)[0] + '.json'
        with open(fname, 'w') as f:
            json.dump(response.json(), f)
        
        return cls(image)

    @property
    def semimajor_radii(self):
        ellipsoid = csmapi.SettableEllipsoid.getEllipsoid(self.model)
        return ellipsoid.getSemiMajorRadius()

    @property
    def semiminor_radii(self):
        ellipsoid = csmapi.SettableEllipsoid.getEllipsoid(self.model)
        return ellipsoid.getSemiMinorRadius()

    @methodispatch
    def intersect(self, height, image_pt):
        if not isinstance(image_pt, csmapi.ImageCoord):
            image_pt = csmapi.ImageCoord(*image_pt)
        return self.model.imageToGround(image_pt, height)

    @intersect.register(GeoDataset)
    def _(self, dem, image_pt, tolerance=1, max_its=20):
        if not isinstance(image_pt, csmapi.ImageCoord):
            # Support a call where image_pt is in the form (x,y)
            image_pt = csmapi.ImageCoord(*image_pt)

        intersection = self.camera.imageToGround(image_pt, 0.0)
        iterations = 0
        while iterations < max_its:
            height = 
            next_intersection = camera.imageToGround(image_pt, height)
            dist = max(abs(intersection.x - next_intersection.x),
                    abs(intersection.y - next_intersection.y),
                    abs(intersextion.z - next_intersection.z))

            intersection = next_intersection
            if dist < tolerance:
                return intersection
            iteration += 1


class Footprint():
    def __init__(self, sensor, geodata, dem=0):
        self.sensor = sensor
        self.dem = dem


    def footprint_bcbf(self):
        '''
        Generates a latlon bounding box given a camera model

        Parameters
        ----------
        camera : object
                csmapi generated camera model
        boundary : list
                of boundary image coordinates
        radii : tuple
                in the form (semimajor, semiminor) axes in meters. The default
                None, will attempt to get the radii from the camera object.
    
        Returns
        -------
        : ndarray
        of ground coordinates
        '''

        gnds = []

        for i, b in enumerate(self.geodata.boundary):
            # Could be potential errors or warnings from imageToGround
            try:
                gnd = self.sensor.intersect(b, self.dem)
                gnds.append([gnd.x, gnd.y, gnd.z])
            except: pass

        return np.array(gnds)

def footprint_latlon(ground_points, radii=None):
    '''
    Generates a latlon footprint from a csmapi generated camera model
    Parameters
    ----------
    camera : object
             csmapi generated camera model
    nnodes : int
             Not sure
    semi_major : int
                 Semimajor axis of the target body
    semi_minor : int
                 Semiminor axis of the target body
    n_points : int
               Number of points to generate between the corners of the bounding
               box per side.
    Returns
    -------
    : object
      ogr multipolygon containing between one and two polygons
    '''
    if radii is None:
        semi_major, semi_minor = get_radii(camera)
    else:
        semi_major, semi_minor = radii

    lons, lats, _ = generate_latlon_boundary(camera, boundary, radii=radii)

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
    nnodes : int
             Not sure
    n_points : int
               Number of points to generate between the corners of the bounding
               box per side.
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

    ecef = pyproj.Proj(proj='geocent', a=semi_major, b=semi_minor)
    lla = pyproj.Proj(proj='latlon', a=semi_major, b=semi_minor)

    # Step over all geometry objects in the latlon footprint
    for i in range(latlon_fp.GetGeometryCount()):
        latlon_coords = np.array(latlon_fp.GetGeometryRef(i).GetGeometryRef(0).GetPoints())

        # Check if the geometry object is populated with points
        if len(latlon_coords) > 0:
            x, y, z = pyproj.transform(lla, ecef,  latlon_coords[:,0], latlon_coords[:,1], latlon_coords[:,2])

            # Step over all coordinate points in a geometry object and update said point
            for j, _ in enumerate(latlon_coords):
                latlon_fp.GetGeometryRef(i).GetGeometryRef(0).SetPoint(j, x[j], y[j], 0)

    return latlon_fp

def generate_gcps(camera, boundary, radii=None):
    '''
    Generates an area of ground control points formated as:
    <GCP Id="" Info="" Pixel="" Line="" X="" Y="" Z="" /> per record
    Parameters
    ----------
    camera : object
             csmapi generated camera model
    boundary : list
               of image boundary coordinates
    nnodes : int
             Not sure
    semi_major : int
                 Semimajor axis of the target body
    semi_minor : int
                 Semiminor axis of the target body
    n_points : int
               Number of points to generate between the corners of the bounding
               box per side.
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




