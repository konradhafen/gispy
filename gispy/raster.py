from osgeo import gdal, ogr, osr
import numpy as np
import vector
from scipy import stats
import os

def raster_test():
    return "this is the raster module of the gispy package"

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
    if xres is None or yres is None: xres, yres = getXYResolution(rasterpath)
    if field is not None and fieldValue is not None: wexp = str(field) + " = \'" + str(fieldValue) + "\'"
    else: wexp = None

    warpOptions = gdal.WarpOptions(format='GTiff', cutlineDSName=polygonpath, cropToCutline=True, cutlineWhere=wexp, xRes=xres, yRes=abs(yres), dstNodata=nodata)
    gdal.WarpOptions()
    gdal.Warp(outputpath, rasterpath, options=warpOptions)
    return None

def createBandIndex(rasterPath, minValue, maxValue):
    """

    Args:
        rasterPath:
        minValue:
        maxValue:

    Returns:

    """
    band = getRasterBandAsArray(rasterPath)
    if band is not None:
        array = band - minValue
        array = np.where((array < 0) | (array > (maxValue-minValue)), 0, array)
        return array
    else:
        return None

def createMask(rasterPath, minValue, maxValue, band=1):
    """
    Create a mask from a raster band

    Args:
        rasterPath: path to raster
        minValue: minimum data value
        maxValue: maximum data value
        band: raster band to use (default: 1)

    Returns:
        Integer array with a value of 0 where the input raster band is < minValue or > maxValue

    """
    array = getRasterBandAsArray(rasterPath, band)
    return np.where((array < minValue) | (array > maxValue), 0, 1)

def getGeoTransform(rasterPath):
    """
    Get the affine geotransformation information for a raster dataset

    Args:
        rasterPath: path to rater dataset

    Returns: 6 element list if successful, None if not successful

    """
    ds = gdal.Open(rasterPath)
    if ds is not None:
        geot = ds.GetGeoTransform()
        ds = None
        return geot
    else:
        return None

def getGeoTransformAndSize(rasterPath):
    """
    Get affine transformation information and the number of rows and columns for a raster dataset.

    Args:
        rasterPath: Path to raster dataset.

    Returns:
        int, int, list: number of rows, number of columns, geotransform (six item list)

    """
    geot = getGeoTransform(rasterPath)
    if geot is not None:
        ds = gdal.Open(rasterPath)
        if ds is not None:
            rows = ds.RasterYSize
            cols = ds.RasterXSize
            return rows, cols, geot
        else:
            return None, None, geot
    else:
        return None, None, None

def getProjection(rasterPath):
    ds = openGDALRaster(rasterPath)
    return ds.GetProjection()

def getRasterAsArray(rasterPath):
    """
    Returns raster as numpy array

    Args:
        rasterPath: Path to input rater

    Returns:
        numpy array if file exists or None if it does not

    """
    if os.path.isfile(rasterPath):
        ds = gdal.Open(rasterPath)
        array = ds.ReadAsArray()
        ds = None
        return array
    else:
        return None

def getRasterBandAsArray(rasterPath, band = 1):
    """
    Returns raster band as 2d numpy array

    Args:
        rasterPath: Path to input raster

    Returns:
        numpy array if band exists or None if band does not exist

    """
    if os.path.isfile(rasterPath):
        ds = gdal.Open(rasterPath)
        bands = ds.RasterCount
        if band > 0 and band <= bands:
            array = ds.GetRasterBand(band).ReadAsArray()
            ds = None
            return array
        else:
            return None
    else:
        return None

def getXYResolution(rasterPath):
    """
    Get X and Y pixel resolution from a rater dataset

    Args:
        rasterPath: path to raster dataset

    Returns:
        X resolution (positive), Y resolution (negative) on success or None on failure

    """
    geot = getGeoTransform(rasterPath)
    if geot is not None: return geot[1], geot[5]
    else: return None

def linearTake(values, indices):
    """
    Get 2d array of band values from a multiband raster of shape (bands, rows, columns)
    Args:
        values: input 3d array containing data values
        indices: 2d array containing indices to bands in the value raster

    Returns:
        2d array where the value in the array corresponds the band value from indices

    """
    _, nR, nC = values.shape
    idx = nC*nR*indices + nC*np.arange(nR)[:, None] + np.arange(nC) #convert 2d indices to linear indices
    return np.take(values, idx), idx

def maskArray(array, mask, nodata=-9999):
    """
    Replace array values with a no data value where a mask is false (0)

    Args:
        array: input array to be masked
        mask: boolean array or integer array with values of 0 (false) and 1 (true)
        nodata: value to write where mask is false (default: -9999)

    Returns:
        array with nodata where mask is false

    """
    if array.shape != mask.shape:
        print "error masking array: array shapes are different", array.shape, mask.shape
    return np.where(mask, array, nodata)

def maskRaster(rasterPath, array, nodata=-9999, band=1):
    """
    Replace values in a raster band with no data where another array is equal to no data

    Args:
        rasterPath: Path of raster to mask.
        array: Array to mask with (must be same shape as array from raster).
        nodata: No data value of array (default: -9999).
        band: Band of raster to mask (default: 1).

    Returns:

    """
    mask = np.where(array==nodata, 0, 1)
    ds = openGDALRaster(rasterPath, gdal.GA_Update)
    band = ds.GetRasterBand(band).ReadAsArray()
    ds.GetRasterBand(1).WriteArray(maskArray(band, mask, nodata))
    ds.GetRasterBand.SetNoDataValue(nodata)
    ds = None
    return None

def maskRasterWithRaster(inputraster, maskraster, inputband=1, maskband=1):
    """
    Mask raster with another raster according to no data extent.

    Args:
        inputraster: Path of raster to be masked.
        maskraster: Path of raster to use as mask.
        inputband: Band of input raster to be masked (default: 1)
        maskband: Band of mask raster to use as mask (default: 1)

    Returns:

    """
    ds = openGDALRaster(inputraster, gdal.GA_Update)
    dsmask = openGDALRaster(maskraster)
    band = ds.GetRasterBand(inputband).ReadAsArray()
    maskarray = dsmask.GetRasterBand(maskband).ReadAsArray()
    nodata = dsmask.GetRasterBand(maskband).GetNoDataValue()
    ds.GetRasterBand(inputband).WriteArray(maskArray(band, np.where(maskarray==nodata, 0, 1), nodata))
    ds.GetRasterBand(inputband).SetNoDataValue(nodata)
    ds = None
    dsmask = None

def openGDALRaster(rasterPath, access=gdal.GA_ReadOnly):
    ds = gdal.Open(rasterPath, access)
    if ds is not None:
        return ds

def percentileMultiband(multi, index):
    """
    Calculate the percentile of a specified value at a position in a multiband raster

    Args:
        multi: Path to multiband raster
        index: Numpy array containing the band index

    Returns:
        2d arrays result (percentile of value based on all bands), and score (value of requested band)

    """
    mean = np.mean(multi, axis=0) #mean of all bands at each row,col
    sd = np.std(multi, axis=0) #standard deviation of all bands at each row,col
    score, idx = linearTake(multi, index) #value of index band at each row,col
    result = stats.norm.cdf(score, loc=mean, scale=sd)*100 #percentile
    return result, score

def percentileOfMultibandIndex(datapath, index, percentilepath, scorepath=None, mask = None):
    """
    For a multiband raster, calculates the percentile of a band value at each row,col

    Args:
        datapath: input path to multiband raster
        index: 2d array with each row,col containing an index to a band in data path
        percentilepath: output path for 2d raster of percentile values
        scorepath: ouput path for 2d raster of the scores of each band index (optional)
        maskpath: boolean array to mask outputs (optional)

    Returns:
        None

    """
    multi = getRasterAsArray(datapath)
    if index is not None and index.shape == multi.shape[1:]:
        result, score = percentileMultiband(multi, index)
        if mask is None or mask.shape != result.shape:
            mask = np.ones((result.shape))

        ds = gdal.Open(datapath)
        writeArrayAsGTiff(percentilepath, maskArray(result, mask), rows=result.shape[0], cols=result.shape[1],
                          geot=ds.GetGeoTransform(), srs=ds.GetProjection())
        if scorepath is not None:
            writeArrayAsGTiff(scorepath, maskArray(score, mask), rows=result.shape[0], cols=result.shape[1],
                              geot=ds.GetGeoTransform(), srs=ds.GetProjection())
        ds = None
    else:
        print "problem with input index array", index.shape, multi.shape[1:]

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
        allcells: If all cells intersected by polygons should be rasterized, or just when polygon includes cell center (defaul: False)
        datatype: GDAL datatype of new raster (default: gdal.GDT_Float32)

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

def writeArrayAsGTiff(path, array, rows, cols, geot, srs, nodata=-9999, nan=-9999, type=gdal.GDT_Float32):
    """
    Write array to a GeoTiff raster

    Args:
        path: output file for raster
        array: array containing data
        rows: number of rows in array
        cols: number of columns in array
        geot: affine geotransformation for the output raster
        srs: spatial reference for the output raster
        nodata: no data value for the output raster
        nan: value in array that should be written as nodata
        type: gdal data type of output raster (default: GDT_Float32)

    Returns:
        None

    """
    ds = gdal.GetDriverByName("GTiff").Create(path, xsize=cols, ysize=rows, bands=1, eType=type)
    ds.SetProjection(srs)
    ds.SetGeoTransform(geot)
    array = np.where((array==np.nan) | (array==nan), nodata, array)
    ds.GetRasterBand(1).WriteArray(array)
    ds.GetRasterBand(1).SetNoDataValue(nodata)
    ds.GetRasterBand(1).FlushCache()
    ds = None
    return None