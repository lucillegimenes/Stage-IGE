# Built ins
import logging
import warnings

# External libs
import numpy as np
from scipy.interpolate import griddata
from scipy import optimize

# Locals
from oggm import utils, cfg
from oggm import entity_task
from oggm.core.gis import gaussian_blur
from oggm.exceptions import InvalidParamsError, InvalidWorkflowError

# Module logger
log = logging.getLogger(__name__)

try:
    import salem
except ImportError:
    pass


# Module logger
log = logging.getLogger(__name__)


def add_data_thickness(gdir,input_file,name_var,long_name_var):

    dsb = salem.GeoTiff(input_file)
    thick = utils.clip_min(dsb.get_vardata(), 0)

    in_volume = thick.sum() * dsb.grid.dx ** 2
    thick = gdir.grid.map_gridded_data(thick, dsb.grid, interp='linear')


	# Correct for volume
    thick = utils.clip_min(thick.filled(0), 0)
    out_volume = thick.sum() * gdir.grid.dx ** 2
    if out_volume > 0:
        thick *= in_volume / out_volume
    
    # We mask zero ice as nodata
    thick = np.where(thick == 0, np.NaN, thick)

	# Write
    with utils.ncDataset(gdir.get_filepath('gridded_data'), 'a') as nc:

        vn = name_var
        if vn in nc.variables:
            v = nc.variables[vn]
        else:
            v = nc.createVariable(vn, 'f4', ('y', 'x', ), zlib=True)
        v.units = 'm'
        ln = long_name_var
        v.long_name = ln
		#v.base_url = base_url
        v[:] = thick
        
def add_data_elev(gdir,input_file,name_var,long_name_var):
    #same as add_data_thickness but allow the data to be negative
    dsb = salem.GeoTiff(input_file)
    elev = dsb.get_vardata()
    in_volume = elev.sum() * dsb.grid.dx ** 2
    elev = gdir.grid.map_gridded_data(elev, dsb.grid, interp='linear')

	# Correct for volume
    #thick = utils.clip_min(thick.filled(0), 0)
    out_volume = elev.sum() * gdir.grid.dx ** 2
    if out_volume > 0:
        elev *= in_volume / out_volume
    
    # We mask zero ice as nodata
    #thick = np.where(thick == 0, np.NaN, thick)

	# Write
    with utils.ncDataset(gdir.get_filepath('gridded_data'), 'a') as nc:

        vn = name_var
        if vn in nc.variables:
            v = nc.variables[vn]
        else:
            v = nc.createVariable(vn, 'f4', ('y', 'x', ), zlib=True)
        v.units = 'm'
        ln = long_name_var
        v.long_name = ln
		#v.base_url = base_url
        v[:] = elev


def mapped_extracted_thickness(gdir):
    grids_file = gdir.get_filepath('gridded_data')
   # See if we have the masks, else compute them
    with utils.ncDataset(grids_file) as nc:
        has_masks = 'glacier_ext_erosion' in nc.variables
    if not has_masks:
        from oggm.core.gis import gridded_attributes
        gridded_attributes(gdir)

    with utils.ncDataset(grids_file) as nc:
        topo_smoothed = nc.variables['topo_smoothed'][:]
        glacier_mask = nc.variables['glacier_mask'][:]
        dis_from_border = nc.variables['dis_from_border'][:]
    
    # Thickness 
    thick = glacier_mask * np.NaN
    thick_hist= glacier_mask * np.NaN
    thick[glacier_mask == 1] = 0
    
   # Along the lines
    cls = gdir.read_pickle('inversion_output')
    fls = gdir.read_pickle('inversion_flowlines')
    vs = []
    for cl, fl in zip(cls, fls):
        vs.extend(cl['volume'])
        x, y = utils.tuple2int(fl.line.xy)
        thick[y, x] = cl['thick']
        thick_hist[y, x] = cl['thick']
           
     # Re-mask
    thick[glacier_mask == 0] = np.NaN
    thick_hist[glacier_mask == 0] = np.NaN
    
    # Write
    with utils.ncDataset(grids_file, 'a') as nc:
        vn = 'extracted_thickness'
        if vn in nc.variables:
            v = nc.variables[vn]
        else:
            v = nc.createVariable(vn, 'f4', ('y', 'x', ), zlib=True)
        v.units = '-'
        v.long_name = 'extracted_thickness'
        v[:] = thick
        
        vo = 'extracted_thickness_hist'
        if vo in nc.variables:
            w = nc.variables[vo]
        else:
            w = nc.createVariable(vo, 'f4', ('y', 'x', ), zlib=True)
        w.units = '-'
        w.long_name = 'extracted_thickness_for_histogram'
        w[:] = thick_hist