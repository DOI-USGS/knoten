import datetime
import json
import os

from csmapi import csmapi
import jinja2

import numpy as np
import scipy.stats
import pyproj
import gdal
from gdal import ogr
import pvl


def generate_boundary(isize, nnodes=5, n_points=10):
    '''
    Generates a bounding box given a camera model

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
    boundary : list
               List of full bounding box
    '''
    x = np.linspace(0, isize[0], n_points)
    y = np.linspace(0, isize[1], n_points)
    boundary = [(i,0.) for i in x] + [(isize[0], i) for i in y[1:]] +\
               [(i, isize[1]) for i in x[::-1][1:]] + [(0.,i) for i in y[::-1][1:]]

    return boundary

def generate_latlon_boundary(camera, nnodes=5, semi_major=3396190, semi_minor=3376200, n_points=10):
    '''
    Generates a latlon bounding box given a camera model

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
    lons : list
           List of longitude values

    lats : list
           List of latitude values

    alts : list
           List of altitude values
    '''
    isize = camera.getImageSize()
    isize = (isize.line, isize.samp)

    boundary = generate_boundary(isize, nnodes=nnodes, n_points = n_points)

    ecef = pyproj.Proj(proj='geocent', a=semi_major, b=semi_minor)
    lla = pyproj.Proj(proj='latlon', a=semi_major, b=semi_minor)

    gnds = np.empty((len(boundary), 3))

    for i, b in enumerate(boundary):
        # Could be potential errors or warnings from imageToGround
        try:
            gnd = camera.imageToGround(csmapi.ImageCoord(*b), 0)
        except:
            pass

        gnds[i] = [gnd.x, gnd.y, gnd.z]

    lons, lats, alts = pyproj.transform(ecef, lla, gnds[:,0], gnds[:,1], gnds[:,2])
    return lons, lats, alts

def generate_gcps(camera, nnodes=5, semi_major=3396190, semi_minor=3376200, n_points=10):
    '''
    Generates an area of ground control points formated as:
    <GCP Id="" Info="" Pixel="" Line="" X="" Y="" Z="" /> per record

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
    gcps : list
           List of all gcp records generated
    '''
    lons, lats, alts = generate_latlon_boundary(camera, nnodes=nnodes,
                                                semi_major=semi_major,
                                                semi_minor=semi_minor,
                                                n_points=n_points)

    lla = np.vstack((lons, lats, alts)).T

    tr = zip(boundary, lla)

    gcps = []
    for i, t in enumerate(tr):
        l = '<GCP Id="{}" Info="{}" Pixel="{}" Line="{}" X="{}" Y="{}" Z="{}" />'.format(i, i, t[0][1], t[0][0], t[1][0], t[1][1], t[1][2])
        gcps.append(l)

    return gcps

def generate_latlon_footprint(camera, nnodes=5, semi_major=3396190, semi_minor=3376200, n_points=10):
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
    lons, lats, _ = generate_latlon_boundary(camera, nnodes=nnodes,
                                                semi_major=semi_major,
                                                semi_minor=semi_minor,
                                                n_points=n_points)

    # Transform coords from -180, 180 to 0, 360
    # Makes crossing the maridian easier to identify
    ll_coords = [*zip(((lons + 180) % 360), lats)]

    ring = ogr.Geometry(ogr.wkbLinearRing)
    wrap_ring = ogr.Geometry(ogr.wkbLinearRing)

    poly = ogr.Geometry(ogr.wkbPolygon)
    wrap_poly = ogr.Geometry(ogr.wkbPolygon)

    multipoly = ogr.Geometry(ogr.wkbMultiPolygon)

    current_ring = ring
    switch_point = None
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

def generate_bodyfixed_footprint(camera, nnodes=5, semi_major=3396190, semi_minor=3376200, n_points=10):
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
    latlon_fp = generate_latlon_footprint(camera, nnodes=5, semi_major = semi_major, semi_minor = semi_minor, n_points = n_points)

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

def generate_vrt(camera, raster_size, fpath, outpath=None, no_data_value=0):
    gcps = generate_gcps(camera)
    xsize, ysize = raster_size

    if outpath is None:
        outpath = os.path.dirname(fpath)
    outname = os.path.splitext(os.path.basename(fpath))[0] + '.vrt'
    outname = os.path.join(outpath, outname)

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
               'proj':'+proj=longlat +a=3396190 +b=3376200 +no_defs',
               'fpath':fpath,
               'no_data_value':no_data_value}
    template = jinja2.Template(vrt)
    tmp = template.render(context)
    warp_options = gdal.WarpOptions(format='VRT', dstNodata=0)
    gdal.Warp(outname, tmp, options=warp_options)
