import logging
import warnings
import pickle
import os
import shutil
from collections.abc import Sequence

import numpy as np
import pandas as pd
from scipy.interpolate import griddata
from scipy import optimize as optimization

import oggm
from oggm import cfg, utils, workflow, tasks, graphics, core

from oggm.core import centerlines, flowline, climate
from oggm.exceptions import InvalidParamsError, InvalidWorkflowError
from oggm.utils import global_task
from oggm.core.inversion import _vol_below_water 

import multiprocessing

#MPI
try:
    import oggm.mpi as ogmpi
    _have_ogmpi = True
except ImportError:
    _have_ogmpi = False

# Module logger
log = logging.getLogger(__name__)

# Multiprocessing Pool
_mp_manager = None
_mp_pool = None

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
    data_file : the .pkl file containing the thicknesses along the flowlines
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


def calibrate_inversion_from_millan(gdirs,path_to_hdf, vol_bias=1, ignore_missing=True,
                                       fs=0, a_bounds=(0.1, 10),
                                       apply_fs_on_mismatch=False,
                                       error_on_mismatch=True,
                                       filter_inversion_output=True):
    """Fit the total volume of the glaciers to the Millan estimate.
    This method finds the "best Glen A" to match all glaciers in gdirs with
    a valid inverted volume.
    Parameters
    ----------
    gdirs : list of :py:class:`oggm.GlacierDirectory` objects
        the glacier directories to process
    ignore_missing : bool
        set this to true to silence the error if some glaciers could not be
        found in the consensus estimate.
    fs : float
        invert with sliding (default: no)
    a_bounds: tuple
        factor to apply to default A
    apply_fs_on_mismatch: false
        on mismatch, try to apply an arbitrary value of fs (fs = 5.7e-20 from
        Oerlemans) and try to optimize A again.
    error_on_mismatch: bool
        sometimes the given bounds do not allow to find a zero mismatch:
        this will normally raise an error, but you can switch this off,
        use the closest value instead and move on.
    filter_inversion_output : bool
        whether or not to apply terminus thickness filtering on the inversion
        output (needs the downstream lines to work).
    Returns
    -------
    a dataframe with the individual glacier volumes
    """

    gdirs = utils.tolist(gdirs)

    # Get the ref data for the glaciers we have
    df = pd.read_hdf(path_to_hdf)
    rids = [gdir.rgi_id for gdir in gdirs]

    found_ids = df.index.intersection(rids)
    if not ignore_missing and (len(found_ids) != len(rids)):
        raise InvalidWorkflowError('Could not find matching indices in the '
                                   'millan estimate for all provided '
                                   'glaciers. Set ignore_missing=True to '
                                   'ignore this error.')

    df = df.reindex(rids)
    df= df*vol_bias 

    # Optimize the diff to ref
    def_a = cfg.PARAMS['inversion_glen_a']

    def compute_vol(x):
        workflow.inversion_tasks(gdirs, glen_a=x*def_a, fs=fs,
                        filter_inversion_output=filter_inversion_output)
        odf = df.copy()
        odf['oggm'] = workflow.execute_entity_task(tasks.get_inversion_volume, gdirs)
        return odf.dropna()

    def to_minimize(x):
        log.workflow('Consensus estimate optimisation with '
                     'A factor: {} and fs: {}'.format(x, fs))
        odf = compute_vol(x)
        return odf.volume_m3.sum() - odf.oggm.sum()

    try:
        out_fac, r = optimization.brentq(to_minimize, *a_bounds, rtol=1e-2,
                                         full_output=True)
        if r.converged:
            log.workflow('calibrate_inversion_from_consensus '
                         'converged after {} iterations and fs={}. The '
                         'resulting Glen A factor is {}.'
                         ''.format(r.iterations, fs, out_fac))
        else:
            raise ValueError('Unexpected error in optimization.brentq')
    except ValueError:
        # Ok can't find an A. Log for debug:
        odf1 = compute_vol(a_bounds[0]).sum() * 1e-9
        odf2 = compute_vol(a_bounds[1]).sum() * 1e-9
        msg = ('calibration from millan estimate CAN\'T converge with fs={}.\n'
               'Bound values (km3):\nRef={:.3f} OGGM={:.3f} for A factor {}\n'
               'Ref={:.3f} OGGM={:.3f} for A factor {}'
               ''.format(fs,
                         odf1.volume_m3, odf1.oggm, a_bounds[0],
                         odf2.volume_m3, odf2.oggm, a_bounds[1]))
        if apply_fs_on_mismatch and fs == 0 and odf2.oggm > odf2.volume_m3:
            return calibrate_inversion_from_millan(gdirs,path_to_hdf,vol_bias=vol_bias,
                                                      ignore_missing=ignore_missing,
                                                      fs=5.7e-20, a_bounds=a_bounds,
                                                      apply_fs_on_mismatch=False,
                                                      error_on_mismatch=error_on_mismatch)
        if error_on_mismatch:
            raise ValueError(msg)

        out_fac = a_bounds[int(abs(odf1.volume_m3 - odf1.oggm) >
                               abs(odf2.volume_m3 - odf2.oggm))]
        log.workflow(msg)
        log.workflow('We use A factor = {} and fs = {} and move on.'
                     ''.format(out_fac, fs))

    # Compute the final volume with the correct A
    workflow.inversion_tasks(gdirs, glen_a=out_fac*def_a, fs=fs,
                    filter_inversion_output=filter_inversion_output)
    df['vol_oggm_m3'] = workflow.execute_entity_task(tasks.get_inversion_volume, gdirs)
    return df

def calibrate_inversion_from_consensus(gdirs, vol_bias=1, ignore_missing=True,
                                       fs=0, a_bounds=(0.1, 10),
                                       apply_fs_on_mismatch=False,
                                       error_on_mismatch=True,
                                       filter_inversion_output=True):
    """Fit the total volume of the glaciers to the 2019 consensus estimate.
    This method finds the "best Glen A" to match all glaciers in gdirs with
    a valid inverted volume.
    Parameters
    ----------
    gdirs : list of :py:class:`oggm.GlacierDirectory` objects
        the glacier directories to process
    ignore_missing : bool
        set this to true to silence the error if some glaciers could not be
        found in the consensus estimate.
    fs : float
        invert with sliding (default: no)
    a_bounds: tuple
        factor to apply to default A
    apply_fs_on_mismatch: false
        on mismatch, try to apply an arbitrary value of fs (fs = 5.7e-20 from
        Oerlemans) and try to optimize A again.
    error_on_mismatch: bool
        sometimes the given bounds do not allow to find a zero mismatch:
        this will normally raise an error, but you can switch this off,
        use the closest value instead and move on.
    filter_inversion_output : bool
        whether or not to apply terminus thickness filtering on the inversion
        output (needs the downstream lines to work).
    Returns
    -------
    a dataframe with the individual glacier volumes
    """

    gdirs = utils.tolist(gdirs)

    # Get the ref data for the glaciers we have
    df = pd.read_hdf(utils.get_demo_file('rgi62_itmix_df.h5'))
    rids = [gdir.rgi_id for gdir in gdirs]

    found_ids = df.index.intersection(rids)
    if not ignore_missing and (len(found_ids) != len(rids)):
        raise InvalidWorkflowError('Could not find matching indices in the '
                                   'consensus estimate for all provided '
                                   'glaciers. Set ignore_missing=True to '
                                   'ignore this error.')

    df = df.reindex(rids)
    df= df*vol_bias

    # Optimize the diff to ref
    def_a = cfg.PARAMS['inversion_glen_a']

    def compute_vol(x):
        workflow.inversion_tasks(gdirs, glen_a=x*def_a, fs=fs,
                        filter_inversion_output=filter_inversion_output)
        odf = df.copy()
        odf['oggm'] = workflow.execute_entity_task(tasks.get_inversion_volume, gdirs)
        return odf.dropna()

    def to_minimize(x):
        log.workflow('Consensus estimate optimisation with '
                     'A factor: {} and fs: {}'.format(x, fs))
        odf = compute_vol(x)
        return odf.vol_itmix_m3.sum() - odf.oggm.sum()

    try:
        out_fac, r = optimization.brentq(to_minimize, *a_bounds, rtol=1e-2,
                                         full_output=True)
        if r.converged:
            log.workflow('calibrate_inversion_from_consensus '
                         'converged after {} iterations and fs={}. The '
                         'resulting Glen A factor is {}.'
                         ''.format(r.iterations, fs, out_fac))
        else:
            raise ValueError('Unexpected error in optimization.brentq')
    except ValueError:
        # Ok can't find an A. Log for debug:
        odf1 = compute_vol(a_bounds[0]).sum() * 1e-9
        odf2 = compute_vol(a_bounds[1]).sum() * 1e-9
        msg = ('calibration from consensus estimate CAN\'T converge with fs={}.\n'
               'Bound values (km3):\nRef={:.3f} OGGM={:.3f} for A factor {}\n'
               'Ref={:.3f} OGGM={:.3f} for A factor {}'
               ''.format(fs,
                         odf1.vol_itmix_m3, odf1.oggm, a_bounds[0],
                         odf2.vol_itmix_m3, odf2.oggm, a_bounds[1]))
        if apply_fs_on_mismatch and fs == 0 and odf2.oggm > odf2.vol_itmix_m3:
            return calibrate_inversion_from_consensus(gdirs,
                                                      ignore_missing=ignore_missing,
                                                      fs=5.7e-20, a_bounds=a_bounds,
                                                      apply_fs_on_mismatch=False,
                                                      error_on_mismatch=error_on_mismatch)
        if error_on_mismatch:
            raise ValueError(msg)

        out_fac = a_bounds[int(abs(odf1.vol_itmix_m3 - odf1.oggm) >
                               abs(odf2.vol_itmix_m3 - odf2.oggm))]
        log.workflow(msg)
        log.workflow('We use A factor = {} and fs = {} and move on.'
                     ''.format(out_fac, fs))

    # Compute the final volume with the correct A
    workflow.inversion_tasks(gdirs, glen_a=out_fac*def_a, fs=fs,
                    filter_inversion_output=filter_inversion_output)
    df['vol_oggm_m3'] = workflow.execute_entity_task(tasks.get_inversion_volume, gdirs)
    return df
