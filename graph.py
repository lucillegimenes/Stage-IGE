"""Useful plotting functions"""
import os
import functools
import logging
from collections import OrderedDict
import itertools
import textwrap
import xarray as xr

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import shapely.geometry as shpg
from matplotlib import cm as colormap

try:
    import salem
except ImportError:
    pass

OGGM_CMAPS = dict()

from oggm.core.flowline import FileModel
from oggm import cfg, utils
from oggm.core import gis

# Module logger
log = logging.getLogger(__name__)


def set_oggm_cmaps():
    # Set global colormaps
    global OGGM_CMAPS
    OGGM_CMAPS['terrain'] = colormap.terrain
    OGGM_CMAPS['section_thickness'] = plt.cm.get_cmap('YlGnBu')
    OGGM_CMAPS['glacier_thickness'] = plt.get_cmap('viridis')


set_oggm_cmaps()

#Not useful anymore
def plot_superpo(gdir,ds,data_source,ax=None,smap=None, linewidth=3, vmax=None):

    with utils.ncDataset(gdir.get_filepath('gridded_data')) as nc:
            topo = nc.variables['topo'][:]

        # Dirty optim
    smap = ds.salem.get_map(countries=False)
    
    try:
        smap.set_topography(topo)
    except ValueError:
        pass

    toplot_th = np.array([])
    toplot_lines = []
    toplot_crs = []
    vol = []
    crs = gdir.grid.center_grid
    geom = gdir.read_pickle('geometries')
    inv = gdir.read_pickle('inversion_output')
    # Plot boundaries
    poly_pix = geom['polygon_pix']
    smap.set_geometry(poly_pix, crs=crs, fc='none', zorder=2,linewidth=.2)
    for l in poly_pix.interiors:
        smap.set_geometry(l, crs=crs, color='black', linewidth=0.5)

        # plot Centerlines
    cls = gdir.read_pickle('inversion_flowlines')
    for l, c in zip(cls, inv):

        smap.set_geometry(l.line, crs=crs, color='gray',linewidth=1.2, zorder=50)
        toplot_th = np.append(toplot_th, c['thick'])
        for wi, cur, (n1, n2) in zip(l.widths, l.line.coords, l.normals):
            line = shpg.LineString([shpg.Point(cur + wi / 2. * n1),shpg.Point(cur + wi / 2. * n2)])
            toplot_lines.append(line)
            toplot_crs.append(crs)
            vol.extend(c['volume'])
    
    dl = salem.DataLevels(cmap='Blues',data=toplot_th, vmin=0, vmax=vmax)
    colors = dl.to_rgb()

    for l, c, crs in zip(toplot_lines, colors, toplot_crs):
        if (data_source=='millan'):
            smap.set_data(ds.millan_thickness)
        else:
            smap.set_data(ds.consensus_ice_thickness)
        smap.set_geometry(l, crs=crs, color=c,
                          linewidth=linewidth, zorder=50)
    
    
    #smap.set_data(ds.millan_thickness)
    #smap.set_cmap('Blues')
    smap.plot(ax)
    smap.set_cmap('Blues')
    smap.visualize(ax=ax,title='original map and thickness along flowlines')

#function for colormmap
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

def plot_inversion_diff(gdirs,pkl_diff, ax=None, smap=None, linewidth=3, vmax=None):
    """Plots the result of the inversion out of a glacier directory."""
    
    gdir = gdirs[0]
    with xr.open_dataset(gdir.get_filepath('gridded_data')) as ds:
        ds = ds.load()

    smap = ds.salem.get_map(countries=False)

    with utils.ncDataset(gdir.get_filepath('gridded_data')) as nc:
        topo = nc.variables['topo'][:]

    # Dirty optim
    try:
        smap.set_topography(topo)
    except ValueError:
        pass

    toplot_th = np.array([])
    toplot_lines = []
    toplot_crs = []
    vol = []
    for gdir in gdirs:
        crs = gdir.grid.center_grid
        geom = gdir.read_pickle('geometries')
        inv = gdir.read_pickle('inversion_output')
        # Plot boundaries
        poly_pix = geom['polygon_pix']
        smap.set_geometry(poly_pix, crs=crs, fc='none', zorder=2,
                          linewidth=.2)
        for l in poly_pix.interiors:
            smap.set_geometry(l, crs=crs, color='black', linewidth=0.5)

        # plot Centerlines
        cls = gdir.read_pickle('inversion_flowlines')
        for l, c in zip(cls, pkl_diff):

            smap.set_geometry(l.line, crs=crs, color='gray',
                              linewidth=1.2, zorder=50)
            toplot_th = np.append(toplot_th, c['thick'])
            for wi, cur, (n1, n2) in zip(l.widths, l.line.coords, l.normals):
                line = shpg.LineString([shpg.Point(cur + wi / 2. * n1),
                                        shpg.Point(cur + wi / 2. * n2)])
                toplot_lines.append(line)
                toplot_crs.append(crs)
            vol.extend(c['volume'])

    
    orig_cmap = matplotlib.cm.RdBu
    shifted_cmap = shiftedColorMap(orig_cmap, midpoint=0.5, name='shifted')
    
    dl = salem.DataLevels(cmap=shifted_cmap,data=toplot_th, vmin=-20, vmax=20)
    colors_ = dl.to_rgb()
    
    for l, c, crs in zip(toplot_lines, colors_, toplot_crs):
        smap.set_geometry(l, crs=crs, color=c,
                          linewidth=linewidth, zorder=50)

    smap.plot(ax)
    smap.set_cmap(shifted_cmap)
    dl.append_colorbar(ax,label='ice thickness (m)')
    #smap.append_colorbar(ax, label='ice thickness (m)')
    smap.visualize(ax=ax,title='Difference between thickness before and after calibration')
    
