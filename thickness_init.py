import logging
import warnings
import pickle


# Module logger
log = logging.getLogger(__name__)

import numpy as np
from scipy.interpolate import griddata
from scipy import optimize

from oggm import cfg, utils, workflow, tasks, graphics, core
from oggm.core.inversion import _vol_below_water 

def ice_thickness_from_data(gdir, glen_a=None, fs=None, write=True,
                                filesuffix='', water_level=None,
                                t_lambda=None,data_file=''):
    """ Download from data the glacier thickness along the flowlines
    
    Parameters
    ----------
    gdir : :py:class:`oggm.GlacierDirectory`
        the glacier directory to process
    glen_a : float
        glen's creep parameter A. Defaults to cfg.PARAMS.
    fs : float
        sliding parameter. Defaults to cfg.PARAMS.
    write: bool
        default behavior is to compute the thickness and write the
        results in the pickle. Set to False in order to spare time
        during calibration.
    filesuffix : str
        add a suffix to the output file
    water_level : float
        to compute volume below water level - adds an entry to the output dict
    t_lambda : float
        defining the angle of the trapezoid walls (see documentation). Defaults
        to cfg.PARAMS.
    data_file : the file containing the thicknesses along the flowlines
    """

    # Defaults
    if glen_a is None:
        glen_a = cfg.PARAMS['inversion_glen_a']
    if fs is None:
        fs = cfg.PARAMS['inversion_fs']
    if t_lambda is None:
        t_lambda = cfg.PARAMS['trapezoid_lambdas']


    # Ice flow params
    fd = 2. / (cfg.PARAMS['glen_n']+2) * glen_a
    a3 = fs / fd
    rho = cfg.PARAMS['ice_density']

    # Inversion with shape factors?
    sf_func = None
    use_sf = cfg.PARAMS.get('use_shape_factor_for_inversion', None)
    if use_sf == 'Adhikari' or use_sf == 'Nye':
        sf_func = utils.shape_factor_adhikari
    elif use_sf == 'Huss':
        sf_func = utils.shape_factor_huss

    # Clip the slope, in rad
    min_slope = 'min_slope_ice_caps' if gdir.is_icecap else 'min_slope'
    min_slope = np.deg2rad(cfg.PARAMS[min_slope])

    out_volume = 0.

    cls = gdir.read_pickle('inversion_input')
    
    with open(data_file, 'rb') as db_file:
        thick_data = pickle.load(db_file)

    k=-1 #future loop index for thick_data
    for cl in cls:
        k+=1
        
        # Clip slope to avoid negative and small slopes
        slope = cl['slope_angle']
        slope = utils.clip_array(slope, min_slope, np.pi/2.)

        # Glacier width
        w = cl['width']
        a0s = - cl['flux_a0'] / ((rho*cfg.G*slope)**3*fd)
        sf = np.ones(slope.shape)  # Default shape factor is 1

        #Glacier thickness
                #Glacier thickness
        #temporaire (supprimez dernier élément si un en trop)
        while (len(thick_data[k])!=len(cl['width'])):
            thick_data[k]=np.delete(thick_data[k],len(thick_data[k])-1)
        
        out_thick = thick_data[k]

        # volume
        is_rect = cl['is_rectangular']
        fac = np.where(is_rect, 1, 2./3.)
        volume = fac * out_thick * w * cl['dx']

        # Now recompute thickness where parabola is too flat
        is_trap = cl['is_trapezoid']
        if cl['invert_with_trapezoid']:
            min_shape = cfg.PARAMS['mixed_min_shape']
            bed_shape = 4 * out_thick / w ** 2
            is_trap = ((bed_shape < min_shape) & ~ cl['is_rectangular'] & (cl['flux'] > 0)) | is_trap
            for i in np.where(is_trap)[0]:
                try:
                    out_thick[i] = thick_data[k][i]
                    sect = (2*w[i] - t_lambda * out_thick[i]) / 2 * out_thick[i]
                    volume[i] = sect * cl['dx']
                except ValueError:
                    # no solution error - we do with rect
                    out_thick[i] = thick_data[k][i]
                    is_rect[i] = True
                    is_trap[i] = False
                    volume[i] = out_thick[i] * w[i] * cl['dx']

        # Sanity check
        if np.any(out_thick <= 0):
            log.warning("Found zero or negative thickness: "
                        "this should not happen.")

        if write:
            cl['is_trapezoid'] = is_trap
            cl['is_rectangular'] = is_rect
            cl['thick'] = out_thick
            cl['volume'] = volume

            # volume below sl
            try:
                bed_h = cl['hgt'] - out_thick
                bed_shape = 4 * out_thick / w ** 2
                if np.any(bed_h < 0):
                    cl['volume_bsl'] = _vol_below_water(cl['hgt'], bed_h,
                                                        bed_shape, out_thick,
                                                        w,
                                                        cl['is_rectangular'],
                                                        cl['is_trapezoid'],
                                                        fac, t_lambda,
                                                        cl['dx'], 0)
                if water_level is not None and np.any(bed_h < water_level):
                    cl['volume_bwl'] = _vol_below_water(cl['hgt'], bed_h,
                                                        bed_shape, out_thick,
                                                        w,
                                                        cl['is_rectangular'],
                                                        cl['is_trapezoid'],
                                                        fac, t_lambda,
                                                        cl['dx'],
                                                        water_level)
            except KeyError:
                # cl['hgt'] is not available on old prepro dirs
                pass

        out_volume += np.sum(volume)

    if write:
        gdir.write_pickle(cls, 'inversion_output', filesuffix=filesuffix)
        gdir.add_to_diagnostics('inversion_glen_a', glen_a)
        gdir.add_to_diagnostics('inversion_fs', fs)

    return out_volume
