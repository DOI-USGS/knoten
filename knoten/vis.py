import tempfile

import json
import pvl
import pyproj
import csmapi

from knoten import csm

from numbers import Number

import numpy as np
import pandas as pd

from pysis import isis
from pysis.exceptions import ProcessError

import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.figure_factory as ff

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
    source_pyproj = pyproj.Proj(proj = source_proj, a = semi_major, b = semi_minor)
    dest_pyproj = pyproj.Proj(proj = dest_proj, a = semi_major, b = semi_minor)

    y, x, z = pyproj.transform(source_pyproj, dest_pyproj, record[0], record[1], record[2], **kwargs)

    return y, x, z


def point_info(cube_path, x, y, point_type, allow_outside=False):
    """
    Use Isis's campt to get image/ground point info from an image

    Parameters
    ----------
    cube_path : str
                path to the input cube

    x : float
        point in the x direction. Either a sample or a longitude value
        depending on the point_type flag

    y : float
        point in the y direction. Either a line or a latitude value
        depending on the point_type flag

    point_type : str
                 Options: {"image", "ground"}
                 Pass "image" if  x,y are in image space (sample, line) or
                 "ground" if in ground space (longitude, lattiude)

    Returns
    -------
    : PvlObject
      Pvl object containing campt returns
    """
    point_type = point_type.lower()

    if point_type not in {"image", "ground"}:
        raise Exception(f'{point_type} is not a valid point type, valid types are ["image", "ground"]')


    if isinstance(x, Number) and isinstance(y, Number):
        x, y = [x], [y]

    with tempfile.NamedTemporaryFile("w+") as f:
        # ISIS wants points in a file, so write to a temp file
        if point_type == "ground":
            # campt uses lat, lon for ground but sample, line for image.
            # So swap x,y for ground-to-image calls
            x,y = y,x
        elif point_type == "image":
            # convert to ISIS pixels
            x = np.add(x, .5)
            y = np.add(y, .5)

        f.write("\n".join(["{}, {}".format(xval,yval) for xval,yval in zip(x, y)]))
        f.flush()

        with tempfile.NamedTemporaryFile("r+") as campt_output:
            try:
                isis.campt(from_=cube_path, coordlist=f.name, allowoutside=allow_outside, usecoordlist=True, coordtype=point_type, to=campt_output.name)
            except ProcessError as e:
                warn(f"CAMPT call failed, image: {cube_path}\n{e.stderr}")
                return

            pvlres = pvl.load(campt_output.name)

        if len(x) > 1 and len(y) > 1:
            for r in pvlres:
                # convert all pixels to PLIO pixels from ISIS
                r[1]["Sample"] -= .5
                r[1]["Line"] -= .5
        else:
            pvlres["GroundPoint"]["Sample"] -= .5
            pvlres["GroundPoint"]["Line"] -= .5
    return pvlres


def plot_diff(data, title='diff plot', colx='x', coly='y', coldx='diffx', coldy='diffy', colmag='magnitude', width=500, height=500):
    import matplotlib.cm as cm
    from matplotlib.colors import Normalize

    fig = make_subplots(rows=2, cols=2, column_widths=[0.9, .1], row_width=[.9, .1],
                        shared_xaxes=True, shared_yaxes=True, horizontal_spacing = 0.01, vertical_spacing = 0.01)

    quiver_plot = ff.create_quiver(data[colx],
                           data[coly],
                           data[coldx],
                           data[coldy],
                           scale=1,
                           line_color='#e3838e',
                           name='offset direction',
                           arrow_scale=0.1)

    for i in range(len(quiver_plot.data)):
        quiver_plot.data[i].xaxis='x1'
        quiver_plot.data[i].yaxis='y1'

    quiver_plot.layout.xaxis1.update({'anchor': 'y1'})
    quiver_plot.layout.yaxis1.update({'anchor': 'x1', 'domain': [.55, 1]})

    fig.add_trace(quiver_plot.data[0], row=2, col=1)

    text = [f'{coldx}: {r[coldx]}<br>{coldy}: {r[coldy]}<br>{colmag}: {r[colmag]}' for i,r in data.iterrows()]
    fig.add_trace(go.Scatter(x=data[colx], y=data[coly],
                customdata=data,
                mode='markers',
                name=f'{colx},{coly}',
                hovertext=text,
                marker=dict(
                        color=data[colmag],
                        colorbar=dict(
                            thickness = 5,
                            outlinewidth = 0,
                            ypad=0,
                            title=f'{colmag}'
                        ),
                        colorscale="viridis",
                        reversescale=True
                )), row=2, col=1)

    xavg = data.groupby(colx).apply(np.mean)
    fig.add_trace(go.Scatter(x=xavg.index, y=xavg[colmag],
        customdata=xavg,
        name=f'{colx} mean error',
        mode='lines+markers',
        line_shape='spline',
        line_color='seagreen',
        marker=dict(
                color=xavg[colmag],
                colorscale="viridis",
                reversescale=True
        )), row=1, col=1)

    yavg = data.groupby(coly).apply(np.mean)
    fig.add_trace(go.Scatter(x=yavg[colmag],y=yavg.index,
        customdata=yavg,
        name=f'{coly} mean error',
        mode='lines+markers',
        line_shape='spline',
        line_color='purple',
        marker=dict(
                color=yavg[colmag],
                colorscale="viridis",
                reversescale=True
        )), row=2, col=2)

    fig.update_layout(width=width, height=height, showlegend=True,legend_orientation="h", title_text=title)
    fig.update_yaxes(autorange="reversed", title_text=coly, row=2, col=1)
    fig.update_xaxes(title_text=colx, row=2, col=1)

    fig.update_xaxes(title_text=f'mean error', row=2, col=2)

    return fig


def plot_diff_3d(data, title='3D diff plot', colx='x', coly='y', colz='z', coldx='diffx', coldy='diffy', coldz='diffz', colmag='magnitude', width=500, height=500):

    text = [f'{coldx}: {r[coldx]}<br>{coldy}: {r[coldy]}<br>{coldz}: {r[coldz]}<br>{colmag}: {r[colmag]}' for i,r in data.iterrows()]

    plot_data = {
    "type": "scatter3d",
    "mode": "markers",
    "name": "original",
    "text": text,
    "x": data[colx],
    "y": data[coly],
    "z": data[colz],
    "marker": { "colorscale": 'viridis',
                "opacity": .8,
                "size": 5,
                "color": data[colmag],
                "colorbar":{
                            "thickness": 5,
                            "outlinewidth": 0,
                            "ypad": 0,
                            "title":f'{colmag}'
                }
              }
    }


    layout = {
        "title": title,
        "width": width,
        "height": height,
        "scene": {
            "aspectratio": {"x": 1, "y": 1, "z": 0.8},
        }
    }

    fig = go.Figure(data=[plot_data], layout=layout)

    return fig

def plot_diff_3d_cone(data, title='3D diff plot', colx='x', coly='y', colz='z',
                                                  colu='u', colv='v',colw='w',
                                                  coldx='diffx', coldy='diffy', coldz='diffz',
                                                  coldv = 'diffu', coldu='diffv', coldw='diffw',
                                                  colxyz_mag='xyz_magnitude', coluvw_mag='uvw_magnitude',
                                                  width=500, height=500):
    text = [f'{coldx}: {r[coldx]}<br>\
              {coldy}: {r[coldy]}<br>\
              {coldz}: {r[coldz]}<br>\
              {coldu}: {r[coldu]}<br>\
              {coldv}: {r[coldv]}<br>\
              {coldw}: {r[coldw]}<br>\
              {colxyz_mag}: {r[colxyz_mag]}<br>\
              {coluvw_mag}: {r[coluvw_mag]}' for i,r in data.iterrows()]

    plot_data = {
        "type": "cone",
        "text":text,
        "x": data[colx],
        "y": data[coly],
        "z": data[colz],
        "u": data[colu],
        "v": data[colv],
        "w": data[colw],
        "sizeref": 5,
        "colorscale": 'viridis',
        "colorbar":{
                    "thickness": 5,
                    "outlinewidth": 0,
                    "ypad": 0,
                    "title": f'mean({colxyz_mag},{coluvw_mag})'
        }
    }

    layout = {
        "title": title,
        "width": width,
        "height": height,
        "scene": {
            "aspectratio": {"x": 1, "y": 1, "z": 0.8},
        }
    }

    fig = go.Figure(data=[plot_data], layout=layout)

    return fig

def reprojection_diff(isd, cube, nx=10, ny=50, width=500, height=500):
    """
    """

    isdjson = json.load(open(isd))

    nlines = isdjson['image_lines']
    nsamples = isdjson['image_samples']

    # generate meshgrid
    xs, ys = np.mgrid[0:nsamples:nsamples/nx, 0:nlines:nlines/ny]
    xs, ys = xs.flatten(), ys.flatten()

    csmcam = csm.create_csm(isd)

    # get data for isis image to ground, csm ground to image
    isis_pts = point_info(cube, xs, ys, 'image')
    isisgnds = np.asarray([np.asarray(g[1]['BodyFixedCoordinate'].value)*1000 for g in isis_pts])
    csm_pts = np.asarray([[p.samp, p.line] for p in [csmcam.groundToImage(csmapi.EcefCoord(*bf)) for bf in isisgnds]])

    isis2csm_diff = np.asarray([xs,ys]).T - csm_pts
    isis2csm_diffmag = np.linalg.norm(isis2csm_diff, axis=1)
    isis2csm_angles = np.arctan2(*isis2csm_diff.T[::-1])

    isis2csm_data = np.asarray([csm_pts.T[0], csm_pts.T[1], xs, ys,  isis2csm_diff.T[0], isis2csm_diff.T[1], isis2csm_diffmag, isis2csm_angles]).T
    isis2csm_data = pd.DataFrame(isis2csm_data, columns=['csm sample', 'csm line', 'isis sample','isis line', 'diff sample', 'diff line', 'magnitude', 'angles'])

    isis2csm_plot = plot_diff(isis2csm_data, colx='isis sample', coly='isis line',
                     coldx='diff sample', coldy='diff line',
                     title="ISIS2Ground->CSM2Image", width=width, height=height)


    # get data for csm image to ground, isis ground to image
    csmgnds = np.asarray([[p.x, p.y, p.z] for p in [csmcam.imageToGround(csmapi.ImageCoord(y,x), 0) for x,y in zip(xs,ys)]])
    csmlon, csmlat, _ = reproject(csmgnds.T, isdjson['radii']['semimajor'], isdjson['radii']['semimajor'], 'geocent', 'latlong')

    isis_imgpts = point_info(cube, (csmlon+360)%360, csmlat, 'ground')
    isis_imgpts = np.asarray([(p[1]['Sample'], p[1]['Line']) for p in isis_imgpts])

    csm2isis_diff = np.asarray([xs,ys]).T - isis_imgpts
    csm2isis_diffmag = np.linalg.norm(csm2isis_diff, axis=1)
    csm2isis_angles = np.arctan2(*(csm2isis_diff/csm2isis_diffmag[:,np.newaxis]).T[::-1])
    csm2isis_data = np.asarray([xs, ys, isis_imgpts.T[0], isis_imgpts.T[1], csm2isis_diff.T[0], csm2isis_diff.T[1], csm2isis_diffmag, csm2isis_angles]).T
    csm2isis_data = pd.DataFrame(csm2isis_data, columns=['csm sample', 'csm line', 'isis sample','isis line', 'diff sample', 'diff line', 'magnitude', 'angles'])

    csm2isis_plot = plot_diff(csm2isis_data, colx='csm sample', coly='csm line',
                 coldx='diff sample', coldy='diff line',
                 title="CSM2Ground->ISIS2Image", width=width, height=height)

    # get data for footprint comparison
    isis_lonlat = np.asarray([[p[1]['PositiveEast360Longitude'].value, p[1]['PlanetocentricLatitude'].value] for p in isis_pts])
    csm_lonlat = np.asarray([(csmlon+360)%360, csmlat]).T

    isiscsm_difflatlon = csm_lonlat - isis_lonlat
    isiscsm_difflatlonmag = np.linalg.norm(isiscsm_difflatlon, axis=1)
    isiscsm_angleslatlon = np.arctan2(*isiscsm_difflatlon.T[::-1])
    isiscsm_latlondata = np.asarray([isis_lonlat.T[0], isis_lonlat.T[1], csm_lonlat.T[0], csm_lonlat.T[1], isiscsm_difflatlon.T[0], isiscsm_difflatlon.T[1], isiscsm_difflatlonmag, isiscsm_angleslatlon]).T
    isiscsm_latlondata = pd.DataFrame(isiscsm_latlondata, columns=['isis lon', 'isis lat', 'csm lon','csm lat', 'diff lon', 'diff lat', 'magnitude', 'angles'])

    isiscsm_latlonplot = plot_diff(isiscsm_latlondata, colx='isis lon', coly='isis lat',
                 coldx='diff lon', coldy='diff lat',
                 title="ISIS Lat/Lon vs CSM Lat/Lon", width=width, height=height)


    isiscsm_diffbf = isisgnds - csmgnds
    isiscsm_diffbfmag = np.linalg.norm(isiscsm_diffbf, axis=1)
    isiscsm_anglesbf = np.arctan2(*isiscsm_diffbf.T[::-1])

    isiscsm_bfdata = np.asarray([isisgnds.T[0], isisgnds.T[1], isisgnds.T[2], csmgnds.T[0], csmgnds.T[1], csmgnds.T[2], isiscsm_diffbf.T[0], isiscsm_diffbf.T[1], isiscsm_diffbf.T[2], isiscsm_diffbfmag, isiscsm_anglesbf]).T
    isiscsm_bfdata = pd.DataFrame(isiscsm_bfdata, columns=['isisx', 'isisy', 'isisz', 'csmx','csmy', 'csmz', 'diffx', 'diffy', 'diffz', 'magnitude', 'angles'])
    isiscsm_bfplot = plot_diff_3d(isiscsm_bfdata, colx='isisx', coly='isisy', colz='isisz',
                                  title='ISIS Body-Fixed vs CSM Body-Fixed',  width=width, height=height)

    return isis2csm_plot, csm2isis_plot, isiscsm_latlonplot, isiscsm_bfplot, isis2csm_data, csm2isis_data, isiscsm_latlondata, isiscsm_bfdata


def external_orientation_diff(isd, cube, nx=4, ny=4, width=500, height=500):
    csmcam = csm.create_csm(isd)
    isdjson = json.load(open(isd))
    nlines, nsamples = isdjson['image_lines'], isdjson['image_samples']

    xs, ys = np.mgrid[0:nsamples:nsamples/nx, 0:nlines:nlines/ny]
    xs, ys = xs.flatten(), ys.flatten()

    isis_pts = point_info(cube, xs, ys, "image")
    isis_lv_bf = np.asarray([p[1]['LookDirectionBodyFixed'] for p in isis_pts])
    isis_pos_bf = np.asarray([p[1]['SpacecraftPosition'].value for p in isis_pts])*1000
    isis_ephem_times = np.asarray([p[1]['EphemerisTime'].value for p in isis_pts])

    csm_locus = [csmcam.imageToRemoteImagingLocus(csmapi.ImageCoord(y, x)) for x,y in zip(xs,ys)]
    csm_lv_bf = np.asarray([[lv.direction.x, lv.direction.y, lv.direction.z] for lv in csm_locus])
    csm_pos_bf = np.asarray([[lv.point.x, lv.point.y, lv.point.z] for lv in csm_locus])
    csm_ephem_times = np.asarray([csmcam.getImageTime(csmapi.ImageCoord(y, x)) for x,y in zip(xs,ys)])
    csm_ephem_times += isdjson['center_ephemeris_time']

    csmisis_diff_pos = csm_pos_bf - isis_pos_bf
    csmisis_diff_lv = csm_lv_bf - isis_lv_bf
    csmisis_diff_ephem = csm_ephem_times - isis_ephem_times
    csmisis_diff_pos_mag = np.linalg.norm(csmisis_diff_pos, axis=1)
    csmisis_diff_lv_mag = np.linalg.norm(csmisis_diff_lv, axis=1)

    csmisis_diff_pos_data = np.asarray([csm_pos_bf.T[0], csm_pos_bf.T[1], csm_pos_bf.T[2],
                                      isis_pos_bf.T[0], isis_pos_bf.T[1], isis_pos_bf.T[2],
                                      csmisis_diff_pos.T[0], csmisis_diff_pos.T[1], csmisis_diff_pos.T[2],
                                      csmisis_diff_pos_mag])
    csmisis_diff_pos_data = pd.DataFrame(csmisis_diff_pos_data.T, columns=['csm pos x', 'csm pos y', 'csm pos z',
                                                                       'isis pos x', 'isis pos y', 'isis pos z',
                                                                       'diffx', 'diffy', 'diffz',
                                                                        'magnitude'])

    csmisis_diff_pos_plot = plot_diff_3d(csmisis_diff_pos_data, colx='isis pos x', coly='isis pos y', colz='isis pos z',
                                                title='ISIS CSM Position Difference')

    csmisis_diff_lv_data = np.asarray([isis_pos_bf.T[0], isis_pos_bf.T[1], isis_pos_bf.T[2],
                                      csm_lv_bf.T[0], csm_lv_bf.T[1], csm_lv_bf.T[2],
                                      isis_lv_bf.T[0], isis_lv_bf.T[1], isis_lv_bf.T[2],
                                      csmisis_diff_pos.T[0], csmisis_diff_pos.T[1], csmisis_diff_pos.T[2],
                                      csmisis_diff_lv.T[0], csmisis_diff_lv.T[1], csmisis_diff_lv.T[2],
                                      csmisis_diff_pos_mag, csmisis_diff_lv_mag, isis_ephem_times, csm_ephem_times, csmisis_diff_ephem])
    csmisis_diff_lv_data = pd.DataFrame(csmisis_diff_lv_data.T, columns=['isis pos x', 'isis pos y', 'isis pos z',
                                                                        'csm lv x', 'csm lv y', 'csm lv z',
                                                                       'isis lv x', 'isis lv y', 'isis lv z',
                                                                       'diffx', 'diffy', 'diffz',
                                                                       'diffu', 'diffv', 'diffw',
                                                                        'xyz_magnitude', 'uvw_magnitude', 'isis ephem time', 'csm ephem time', 'diff ephem'])

    csmisis_diff_lv_plot = plot_diff_3d_cone(csmisis_diff_lv_data, colx='isis pos x', coly='isis pos y', colz='isis pos z',
                                                                   colu='isis lv x', colv='isis lv y', colw='isis lv z',
                                                                    title='ISIS CSM Position and Look Vector Difference', width=width, height=height)

    csmisis_diff_ephem_plot = go.Figure(go.Scatter(x=np.linspace(0, nlines, ny), y=csmisis_diff_ephem, line_shape='spline')).update_layout(title='ISIS CSM Ephem Time Difference', width=width, height=height/2).update_xaxes(title_text='Line').update_yaxes(title_text='Time Delta Seconds')

    return csmisis_diff_lv_plot, csmisis_diff_ephem_plot, csmisis_diff_lv_data
