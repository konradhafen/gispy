from osgeo import gdal, ogr, osr
import os

def vector_test():
    return "this is the vector module of the gispy package"

def copyFields(lyr, lyrdefn):
    """
    Create fields in a new layer based on fields from an existing OGRLayerDefn.

    Args:
        lyr: OGRLayer to create fields for.
        lyrdefn: OGRLayerDefn to copy from.

    Returns:
        OGRLayer with fields from lyrdefn added.

    """
    for i in range(0, lyrdefn.GetFieldCount()):
        fldefn = lyrdefn.GetFieldDefn(i)
        lyr.CreateField(fldefn)
    return lyr

def copyFieldValues(newfeat, copyfeat, defn):
    """
    Copy field values from one feature to another.

    Args:
        newfeat: OGRFeature to copy to
        copyfeat: OGRFeature containing field values to copy
        defn: OGRFieldDefinition of copyfeat

    Returns:

    """
    for i in range(0, defn.GetFieldCount()):
        newfeat.SetField(defn.GetFieldDefn(i).GetNameRef(), copyfeat.GetField(i))

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
        if index is not -1:
            lyr.DeleteField(index)
        lyr.CreateField(field)
    return lyr

def createOGRDataSource(filename, driver='ESRI Shapefile'):
    """
    Create OGRDataSource
    Args:
        filename: Path of file to create
        driver: Name of OGR Driver to use (default: 'ESRI Shapefile')

    Returns:
        OGRDataSource

    """
    driver = ogr.GetDriverByName(driver)
    if os.path.exists(filename):
        driver.DeleteDataSource(filename)
    ds = driver.CreateDataSource(filename)
    return ds

def deleteFields(filename, fields):
    """
    Delete fields from shapefile

    Args:
        filename: Path to shapefile.
        fields: List of fields to delete

    Returns:

    """
    ds = gdal.OpenEx(filename, gdal.OF_VECTOR | gdal.OF_UPDATE)
    layername = getFilenameWithoutExtenstion(filename)
    for field in fields:
        sql = 'ALTER TABLE ' + layername + ' DROP COLUMN ' + field
        ds.ExecuteSQL(sql)
    ds = None

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

def getFilenameWithoutExtenstion(filepath):
    """
    Return the name of a file without the director or extension

    Args:
        filepath: Path to file.

    Returns:
        string: name of file without extension

    """
    base = os.path.basename(filepath)
    return os.path.splitext(base)[0]

def joinZonalStatsToSHP(inshp, zsresult, id, stats, fieldnames, stattype=ogr.OFTReal):
    """
    Join zonal stats results to shapefile

    Args:
        inshp: shapefile to add zonal stats to
        zsresult: result of zonal stats (from rasterstats package)
        id: name of zonal stat with FID of feature corresponding to stats
        stats: names of statistics to join
        fieldnames: names of fields to create (corresponding to stats)
        stattype: ogr OFT data type of field to create (default: ogr.OFTReal)

    Returns:
        None

    """
    print "starting join"
    ds = ogr.Open(inshp, 1)
    lyr = ds.GetLayer(0)
    lyr = createFields(lyr, fieldnames, stattype)
    print "fields created"
    for result in zsresult:
        feat = lyr.GetFeature(int(result[id]))
        for i in range(len(stats)):
            feat.SetField(fieldnames[i], result[stats[i]]*1.0)
        lyr.SetFeature(feat)
    ds.Destroy()
    return None

def joinZonalStatsToSHP_rasterstats(inshp, zsresult, id, stats, fieldnames, stattype=ogr.OFTReal):
    """
    Join zonal stats results to shapefile

    Args:
        inshp: shapefile to add zonal stats to
        zsresult: result of zonal stats (from rasterstats package)
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
    """
    Open and existing OGRDataSource.

    Args:
        file: Path of file to open.
        access: Access type, 0 read only, 1 update (default: 0)

    Returns:
        OGRDataSource

    """
    ds = ogr.Open(file, access)
    if ds is not None:
        return ds
    else:
        print 'Problem opening shapefile'

def reprojectShapefileLayer(shapefile, newshapefile, newsrs):
    """
    Create a new shapefile with a different projection.

    Args:
        shapefile: Path of shapefile to reproject.
        newshapefile: Path of shapefile to create with new projection.
        newsrs: OSRSpatialReference of the new shapefile.

    Returns:

    """
    ds = openOGRDataSource(shapefile)
    layer = ds.GetLayer()
    srs = layer.GetSpatialRef()
    transform = osr.CoordinateTransformation(srs, newsrs)

    outds = createOGRDataSource(newshapefile)
    outlyr = outds.CreateLayer(getFilenameWithoutExtenstion(newshapefile), srs=newsrs, geom_type=layer.GetGeomType())
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

def saveFilteredLayerAsShapefile(filename, layer, filter, driver='ESRI Shapefile', close=False):
    """
    Filter OGRLayer with SQL statement and save as a shapefile
    Args:
        filename: Path to save new shapefile.
        layer: Layer to apply SQL filter to.
        filter: Filter to apply to layer.
        driver: Name of OGR driver used to create the new shapefile (default: ESRI Shapefile)
        close: If the new OGRDataSource and OGRLayer should be closed. False: data source and layer returned, True: data source and layer closed (default: False)

    Returns:
        OGRDataSource and OGRLayer if 'close' is False.

    """
    layer.SetAttributeFilter(filter)
    return saveLayerAsShapefile(filename, layer, driver, close)

def saveLayerAsShapefile(filename, layer, driver='ESRI Shapefile', close=False):
    """
    Save a vector layer to a new OGR DataSource.

    Args:
        filename: Name of new data source/shapefile
        layer: Layer to copy.
        driver: Name of driver to create DataSource (default: ESRI Shapefile)
        close: Close created DataSource and layer, if True new ds and new layer will be deleted, if False they will be returned (default: False)

    Returns:
        OGRDataSource and OGRLayer, if 'close' is False.

    """
    ds = createOGRDataSource(filename, driver)
    newlyr = ds.CopyLayer(layer, getFilenameWithoutExtenstion(filename))
    if close:
        del ds, newlyr
        return None
    else:
        return ds, newlyr