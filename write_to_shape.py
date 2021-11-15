import glob
import os
import tempfile
import gzip
import json
import random
import shutil
import tarfile
import sys
import signal
import datetime
import logging
import pickle
import warnings
import itertools
from collections import OrderedDict
from functools import partial, wraps
from time import gmtime, strftime
import fnmatch
import platform
import struct
import importlib
import re as regexp

# External libs
import pandas as pd
import numpy as np
from scipy import stats
import xarray as xr
import shapely.geometry as shpg
from shapely.ops import transform as shp_trafo
import netCDF4

# Optional libs
try:
    import geopandas as gpd
except ImportError:
    pass
try:
    import salem
except ImportError:
    pass
try:
    from salem import wgs84
    from salem.gis import transform_proj
except ImportError:
    pass
try:
    import pyproj
except ImportError:
    pass


# Locals
from oggm import __version__
from oggm.utils._funcs import (calendardate_to_hydrodate, date_to_floatyear,
                               tolist, filter_rgi_name, parse_rgi_meta,
                               haversine, multipolygon_to_polygon, clip_scalar)
from oggm.utils._downloads import (get_demo_file, get_wgms_files,
                                   get_rgi_glacier_entities)

from oggm.utils._workflow import (_write_shape_to_disk,_get_centerline_lonlat)
from oggm import cfg
from oggm.exceptions import InvalidParamsError, InvalidWorkflowError

log = logging.getLogger(__name__)

def write_centerlines_to_shape_bis(gdirs, *, path=True, to_tar=False,
                               filesuffix='', flowlines_output=False,
                               geometrical_widths_output=False,
                               corrected_widths_output=False):
    """Write the centerlines in a shapefile.
    Parameters
    ----------
    gdirs: the list of GlacierDir to process.
    path:
        Set to "True" in order  to store the info in the working directory
        Set to a path to store the file to your chosen location
    to_tar : bool
        put the files in a .tar file. If cfg.PARAMS['use_compression'],
        also compress to .gz
    filesuffix : str
        add suffix to output file
    flowlines_output : bool
        output the flowlines instead of the centerlines
    geometrical_widths_output : bool
        output the geometrical widths instead of the centerlines
    corrected_widths_output : bool
        output the corrected widths instead of the centerlines
    """
    from oggm.workflow import execute_entity_task

    if path is True:
        path = os.path.join(cfg.PATHS['working_dir'],
                            'glacier_centerlines' + filesuffix + '.shp')

    log.workflow('write_centerlines_to_shape on {} ...'.format(path))

    olist = execute_entity_task(_get_centerline_lonlat, gdirs,
                                flowlines_output=flowlines_output,
                                geometrical_widths_output=geometrical_widths_output,
                                corrected_widths_output=corrected_widths_output)
    # filter for none
    olist = [o for o in olist if o is not None]
    odf = gpd.GeoDataFrame(itertools.chain.from_iterable(olist))
    odf = odf.sort_values(by='SEGMENT_ID')
    odf.crs = {'init': 'epsg:4326'}
    _write_shape_to_disk(odf, path, to_tar=to_tar)
