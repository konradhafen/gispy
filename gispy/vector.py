from osgeo import ogr, osr
import os
import geopandas as gpd

def vector_test():
    return "this is the vector module of the gispy package"

def createFields(lyr, fieldnames, fieldtype = ogr.OFTReal):
    """
    Add fields to a shapefile layer

    Args:
        lyr: ogr layer to add fields to
        fieldnames: names of fields to add
        fieldtype: data type of fields to add

    Returns:
        ogr.Layer: input layer with fields added

    """
    for name in fieldnames:
        field = ogr.FieldDefn(name, fieldtype)
        index = lyr.FindFieldIndex(name, 1)
        if index is not -1: lyr.DeleteField(index)
        lyr.CreateField(field)
    return lyr

def copyFields(lyr, lyrdefn):
    for i in range(0, lyrdefn.GetFieldCount()):
        fldefn = lyrdefn.GetFieldDefn(i)
        lyr.CreateField(fldefn)
    return lyr

def copyFieldValues(newfeat, copyfeat, defn):
    for i in range(0, defn.GetFieldCount()):
        newfeat.SetField(defn.GetFieldDefn(i).GetNameRef(), copyfeat.GetField(i))

def createOGRDataSource(filename, driver='ESRI Shapefile', layertype=ogr.wkbPoint, srs=None ):
    driver = ogr.GetDriverByName(driver)
    if os.path.exists(filename):
        driver.DeleteDataSource(filename)
    ds = driver.CreateDataSource(filename)
    lyr = ds.CreateLayer(os.path.split(filename)[0], srs=srs, geom_type=layertype)
    return ds, lyr

def getFieldValues(dataset, fieldName):
    """
    Get list of all values for shapefile field

    Args:
        dataset: path to shapefile
        fieldName: name of field

    Returns:
        list of field values and FIDs if successful, None if not successful

    """
    values = []
    fids = []
    ds = ogr.Open(dataset)
    lyr = ds.GetLayer()
    if lyr is not None:
        for feat in lyr:
            values.append(feat.GetField(fieldName))
            fids.append(feat.GetFID())
        return values, fids
    else:
        return None

def joinZonalStatsToSHP(inshp, zsresult, id, stats, fieldnames, stattype=ogr.OFTReal):
    """
    Join zonal stats results to shapefile

    Args:
        inshp: shapefile to add zonal stats to
        zsresult: result of zonal stats (from rasterstats module)
        id: FID of feature corresponding to stats
        stats: names of statistics to join
        fieldnames: names of fields to create (corresponding to stats)
        stattype: ogr OFT data type of field to create (default: ogr.OFTReal)

    Returns:
        None

    """
    ds = ogr.Open(inshp, 1)
    lyr = ds.GetLayer(0)
    lyr = createFields(lyr, fieldnames, stattype)
    for result in zsresult:
        feat = lyr.GetFeature(int(result[id]))
        for i in range(len(stats)):
            feat.SetField(fieldnames[i], result["properties"][stats[i]])
        lyr.SetFeature(feat)
    ds.Destroy()
    return None

def fieldExists(lyr, fieldname):
    """
    Determine if field exists in a layer.

    Args:
        lyr: OGR Layer.
        fieldname: Name of field.

    Returns:
        bool: True if field exists, False if it does not.

    """
    defn = lyr.GetLayerDefn()
    if defn.GetFieldIndex(fieldname) is not -1:
        return True
    else:
        return False

def openOGRDataSource(file, access=0):
    ds = ogr.Open(file, access)
    if ds is not None:
        return ds
    else:
        print 'Problem opening shapefile'

def reprojectShapefileLayer(shapefile, newshapefile, newsrs):
    ds = openOGRDataSource(shapefile)
    layer = ds.GetLayer()
    srs = layer.GetSpatialRef()
    transform = osr.CoordinateTransformation(srs, newsrs)

    outds, outlyr = createOGRDataSource(newshapefile, srs=newsrs, layertype=layer.GetGeomType())
    copyFields(outlyr, layer.GetLayerDefn())
    defn = outlyr.GetLayerDefn()
    feat = layer.GetNextFeature()
    while feat:
        geom = feat.GetGeometryRef()
        geom.Transform(transform)
        outfeat = ogr.Feature(defn)
        outfeat.SetGeometry(geom)
        copyFieldValues(outfeat, feat, defn)
        outlyr.CreateFeature(outfeat)
        outfeat = None
        feat = layer.GetNextFeature()
    ds = None
    outds = None