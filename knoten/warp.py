import jinja2

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