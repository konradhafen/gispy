from osgeo import gdal, ogr, osr
import os
import raster
import vector


def clipByFeature(inputdir, outputdir, rasterfiles, shapefile, fieldname, nodata=-9999, xres=None, yres=None):
    """
    Clip multiple rasters by features in a shapefile. Creates a new directory for each feature and saves clipped rasters in the directory

    Note:
        All input rasters and shapefiles should have the same spatial reference

    Args:
        inputdir: directory where rasters to be clipped are located
        outputdir: directory where a new directory for each clipped raster will be created
        rasterfiles: list of raster file names and extensions to be clipped. e.g. ['raster1.tif', 'raster2.tif]
        shapefile: shapefile containing features to clip with
        fieldname: name of field to select feature with
        nodata: nodata value (default: -9999)
        xres: x resolution of output raster (default: None), with default resolution is taken from the input rater
        yres: y resolution of output rater (defaul: None), with default resolution is taken from the input rater

    Returns:
        None

    """
    fieldValues, fids = vector.getFieldValues(shapefile, fieldname)
    fieldValues_unique = list(set(fieldValues)) #get unique field values (otherwise the same operation may be done twice)
    for value in fieldValues_unique:
        dirvalue = outputdir + "/" + str(value)
        if not os.path.isdir(dirvalue):
            os.mkdir(dirvalue)
        for raster in rasterfiles:
            infile = inputdir + "/" + str(raster)
            if os.path.exists(infile):
                clipRasterWithPolygon(infile, shapefile, dirvalue + "/" + str(raster), nodata=nodata, xres=xres, yres=yres,
                                      fieldValue=value, field=fieldname)
    return None

def clipRasterWithPolygon(rasterpath, polygonpath, outputpath, nodata=-9999, xres=None, yres=None, field=None, fieldValue=None):
    """

    Args:
        rasterpath: raster to clip
        polygonpath: shapefile containing features to clip with
        outputpath: path of output clipped raster
        nodata: nodata value (default: -9999)
        xres: x resolution of output raster (default: None, resolution of input raster)
        yres: y resolution of output raster (default: None, resolution of input raster)
        field: name of shapefile field to select features and name directories
        fieldValue: list of unique values for input field

    Returns:

    """
    if xres is None or yres is None: xres, yres = raster.getXYResolution(rasterpath)
    if field is not None and fieldValue is not None: wexp = str(field) + " = \'" + str(fieldValue) + "\'"
    else: wexp = None

    warpOptions = gdal.WarpOptions(format='GTiff', cutlineDSName=polygonpath, cropToCutline=True, cutlineWhere=wexp, xRes=xres, yRes=abs(yres), dstNodata=nodata)
    gdal.WarpOptions()
    gdal.Warp(outputpath, rasterpath, options=warpOptions)
    return None

def polygonToRaster(rasterpath, vectorpath, fieldname, rows, cols, geot, prj, allcells=False, nodata=-9999, datatype = gdal.GDT_Float32):
    """
    Convert polygon shapefile to raster dataset.

    Args:
        rasterpath: Path of raster to be created.
        vectorpath: Path of polygon shapefile to rasterize.
        fieldname: Name of shapefile field to use as values in new raster.
        rows: Number of rows in new raster.
        cols: Number of columns in new raster.
        geot: Affine geotransform of new raster.
        prj: Spatial reference of new raster.
        allcells: If all cells intersected by polygons should be rasterized, or just when polygon includes cell center (defaul: False).
        nodata: No data value.
        datatype: GDAL datatype of new raster (default: gdal.GDT_Float32).

    Returns:

    """
    inds = ogr.Open(vectorpath)
    lyr = inds.GetLayer()
    outds = gdal.GetDriverByName('GTiff').Create(rasterpath, cols, rows, 1, datatype)
    outds.SetProjection(prj)
    outds.SetGeoTransform(geot)
    band = outds.GetRasterBand(1).ReadAsArray()
    band.fill(nodata)
    outds.GetRasterBand(1).WriteArray(band)
    outds.GetRasterBand(1).SetNoDataValue(nodata)

    ALL_TOUCHED = 'FALSE'
    if allcells: ALL_TOUCHED = 'TRUE'

    if vector.fieldExists(lyr, fieldname):
        status = gdal.RasterizeLayer(outds, [1], lyr, options=['ALL_TOUCHED='+ALL_TOUCHED, 'ATTRIBUTE='+fieldname, 'NODATA='+str(nodata)])
        if status is not 0:
            print "Rasterize not successful"
    else:
        print "Rasterize field does not exist"

    outds = None
    return None

def zonalStatistics(vectorpath, rasterpath, write=True):
    rasterds = raster.openGDALRaster(rasterpath)
    vectords = vector.openOGRDataSource(vectorpath)
    lyr = vectords.GetLayer()
    geot = rasterds.GetGeoTransform()
    feat = lyr.GetNextFeature()
    while feat:
        tmpds = vector.createOGRDataSource('temp', 'Memory')
        tmplyr = tmpds.CreateLayer('polygons', None, ogr.wkbPolygon)
        tmplyr.CreateFeature(feat.Clone())
    return None