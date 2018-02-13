from osgeo import ogr

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