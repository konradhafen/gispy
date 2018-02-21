from osgeo import gdal, ogr, osr
import numpy as np
import math
import struct
import os
import raster
import vector


def bboxToOffsets(bbox, geot):
    col1 = int((bbox[0] - geot[0]) / geot[1])
    col2 = int((bbox[1] - geot[0]) / geot[1]) + 1
    row1 = int((bbox[3] - geot[3]) / geot[5])
    row2 = int((bbox[2] - geot[3]) / geot[5]) + 1
    return (row1, row2, col1, col2)

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

def polygonToRaster(rasterpath, vectorpath, fieldname, rows, cols, geot, prj=None, drivername='GTiff', allcells=False, nodata=-9999, datatype = gdal.GDT_Float32, islayer=False):
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
        drivername: Name of GDAL driver to use to create raster (default: 'GTiff')
        allcells: If all cells intersected by polygons should be rasterized, or just when polygon includes cell center (defaul: False).
        nodata: No data value.
        datatype: GDAL datatype of new raster (default: gdal.GDT_Float32).
        islayer (bool): True if vector path is an OGRLayer (default: False)

    Returns:

    """
    if islayer:
        lyr = vectorpath
    else:
        inds = ogr.Open(vectorpath)
        lyr = inds.GetLayer()
    driver = gdal.GetDriverByName(drivername)
    outds = driver.Create(rasterpath, cols, rows, 1, datatype)
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

    if inds: inds = None
    outds = None
    return None

def rasterValueAtPoints(pointshapefile, rasterpath, fieldname, datatype=ogr.OFTReal, idxfield=None):
    """
    Get the value of a raster at point locations.

    Args:
        pointshapefile: Path to point shapefile.
        rasterpath: Path to raster dataset.
        fieldname: Name of field to create and write raster values.
        datatype: OGR datatype of created field (default: ogr.OFTReal)
        idxfield: Shapefile field containing band number to use on multipband rasters (default: None, get value from band 1)

    Returns:

    """
    ras = raster.openGDALRaster(rasterpath)
    geot = ras.GetGeoTransform()
    shp = vector.openOGRDataSource(pointshapefile, 1)
    lyr = shp.GetLayer()
    if lyr.GetGeomType() is not ogr.wkbPoint:
        print "incorrect geometry type, should be points", lyr.GetGeomType()
        return None
    vector.createFields(lyr, [fieldname], datatype)
    feat = lyr.GetNextFeature()
    while feat:
        nband = 1
        if idxfield is not None:
            nband = feat.GetField(idxfield)
        geom = feat.GetGeometryRef()
        row, col = raster.getCellAddressOfPoint(geom.GetX(), geom.GetY(), geot)
        if nband <= 0 or nband > ras.RasterCount:
            value = -9999.0
        else:
            value = struct.unpack('f'*1, ras.GetRasterBand(int(nband)).ReadRaster(xoff=col, yoff=row, xsize=1, ysize=1,
                                                                      buf_xsize=1, buf_ysize=1,
                                                                      buf_type=gdal.GDT_Float32))[0]
        feat.SetField(fieldname, value)
        lyr.SetFeature(feat)
        feat = lyr.GetNextFeature()
    lyr = None
    shp.Destroy()

def zonalStatistics(vectorpath, rasterpath, write=['min', 'max', 'sd', 'mean'], prepend=None, idxfield=None):
    rasterds = raster.openGDALRaster(rasterpath)
    vectords = vector.openOGRDataSource(vectorpath)
    lyr = vectords.GetLayer()
    geot = rasterds.GetGeoTransform()
    array = rasterds.ReadAsArray()
    nodata = rasterds.GetRasterBand(1).GetNoDataValue()
    stats=[]
    feat = lyr.GetNextFeature()
    while feat:
        tmpds = vector.createOGRDataSource('temp', 'Memory')
        tmplyr = tmpds.CreateLayer('polygons', None, ogr.wkbPolygon)
        tmplyr.CreateFeature(feat.Clone())
        offsets = bboxToOffsets(feat.GetGeometryRef().GetEnvelope(), geot)
        newgeot= raster.getOffsetGeot(offsets[0], offsets[2], geot)
        tmpras = raster.createGDALRaster('', offsets[1]-offsets[0], offsets[3]-offsets[2], datatype=gdal.GDT_Byte, drivername='MEM', geot=newgeot)
        gdal.RasterizeLayer(tmpras, [1], tmplyr, burn_values=[1])
        tmparray = tmpras.ReadAsArray()
        maskarray = np.ma.MaskedArray(array[offsets[0]:offsets[1], offsets[2]:offsets[3]],
                                      mask=np.logical_or(array[offsets[0]:offsets[1], offsets[2]:offsets[3]]==nodata, np.logical_not(tmparray)))
        featstats = {
            'min' : float(maskarray.min()),
            'mean': float(maskarray.mean()),
            'max': float(maskarray.max()),
            'sd': float(maskarray.std()),
            'sum': float(maskarray.sum()),
            'count': float(maskarray.count()),
            'fid': float(feat.GetFID())
        }
        stats.append(featstats)
        tmpras = None
        tmpds = None
        feat = lyr.GetNextFeature()
    return stats