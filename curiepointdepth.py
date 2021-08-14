# -*- coding: utf-8 -*-
"""
A set of functions to calculate curie point depths from an aeromagnetic\
data grid using the Tanaka method (Tanaka 1999, Tectonophysics).

@author: craigm
Created on Wed May 04 11:09:00 2016

"""

from osgeo import gdal
import os
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from matplotlib.patches import Rectangle
from skimage.exposure import equalize_hist
import pandas as pd
import statsmodels.api as sm
import cm_utils
#import pygmt
from cmcrameri import cm as sci_cm
import geosoft.gxpy as gx

# SET GLOABL FONT PARAMETERS
font = {'family': 'Arial',
        'weight': 'normal',
        'size': 8}
plt.rc('font', **font)


def splitgrid_geosoft(working_dir, gridfilename, window_width, overlap_factor):
    """
    Split the main grid into smaller grids.

    Expected input format is geosoft grid (*.grd)

    This function takes in the large initial grid and splits it into windows
    of specified width and overlap.

    It creates various output directories and makes a plot of each grid subset
    which can be used for checking for empty cells before running the fft.

    Parameters
    ----------
    working_dir: main directory where the main data grid lives and where
                 processing subdirectories will be created.

    gridfilename: main data grid to be split.  Expected format is Surfer 7.

    window_width:  sub grid width in metres (ie 100 km = 100000.)

    overlap_factor: overlap of adjacent grids
    use factor of 2. with 100000. window_width to get 50 km spaced points
    use factor of 4. with 200000. window_width to get 50 km spaced points
    use factor of 6. with 300000. window_width to get 50 km spaced points

    Returns
    -------
    Many subgrids in a sub directory "subsetted grids" under the working_dir.

    """
    gx.gxc = gx.gx.GXpy(log=print)

    dsep = '/'

    main_dir = working_dir

    grid = (main_dir + dsep + gridfilename)
    # grid = 'D:/dropbox/TVZ_Curie_point/Window_64_geosoft.grd'

    # prefix for output grids
    out_dir = main_dir + dsep + 'subsetted_grids'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    if not os.path.exists(out_dir + dsep + 'window_' + str(window_width/1000)):
        os.makedirs(out_dir + dsep + 'window_' + str(window_width/1000),
                    exist_ok=True)

    if not os.path.exists(out_dir + dsep + 'fft_output' +
                          str(window_width/1000)):
        os.makedirs(main_dir + dsep + 'fft_output'
                    + dsep + 'window_' + str(window_width/1000), exist_ok=True)

    if not os.path.exists(out_dir + dsep + 'spectra_plots' +
                          str(window_width/1000)):
        os.makedirs(main_dir + dsep + 'spectra_plots'
                    + dsep + 'window_' + str(window_width/1000), exist_ok=True)

    if not os.path.exists(out_dir + dsep + 'results' + str(window_width/1000)):
        os.makedirs(main_dir + dsep + 'results'
                    + dsep + 'window_' + str(window_width/1000), exist_ok=True)

    outfile = out_dir + dsep + 'window_' + str(window_width/1000) + dsep + \
        'window'

    os.chdir(out_dir + dsep + 'window_' + str(window_width/1000))

    # open grid
    datagrid = gx.grid.Grid.open(grid)

    # setup windowing of the main grid
    pixelwidth = datagrid.dx
    window_step = window_width/pixelwidth/overlap_factor
    window_pixels_x = window_width / datagrid.dx
    window_pixels_y = window_width / datagrid.dy

    print("number of columns : " + str(datagrid.nx))
    print("number of rows : " + str(datagrid.ny))
    print("window step (pixels): " + str(window_step))
    print("\n")

    # set subset counter
    subset = 1
    window = 1

    # process subsets loop moving from left to right and top to bottom
    # while window_uly > grid_lowerrightY:

    center = []
    window_num = []

    for ix in range(0, datagrid.nx, int(window_step)):
        for iy in range(0, datagrid.ny, int(window_step)):
            if (ix + window_pixels_x > datagrid.nx) or (iy + window_pixels_y >
                                                        datagrid.ny):
                print('skipping, window grid outside of main grid')
            else:
                print(f'processing subset # {subset} : X0={ix}, Y0={iy}')

                window_grid = gx.grid.Grid.index_window(datagrid,
                                                        name=outfile + '_' + str(subset) + '.grd',
                                                        x0=ix,
                                                        y0=iy,
                                                        nx=window_pixels_x,
                                                        ny=window_pixels_y,
                                                        overwrite=True)
                # time.sleep(1)
                center.append(window_grid.centroid_xy)
                window_num.append(window)
                subset += 1
                window += 1

    window_num = np.array(window_num)
    center_coords = np.c_[window_num, center]
    np.savetxt('window_' + str(window_width/1000) + '_km_center_coords.csv',
               center_coords, fmt='%.0f %.1f %.1f', delimiter=',')
    time.sleep(2)

    # plot centroid locations
    # convert grid to gxf
    # https://gist.github.com/jesserobertson/59050674e2871ea03acd3fa312ac9c02
    # convert to gxf format then rasterio/GDAL to load and matplotlib to plot

    gxf_file = main_dir + dsep + datagrid.name + '.gxf'
    gx.grid.Grid.copy(datagrid, gxf_file + '(GXF)', overwrite=True)

    gxf_grid = gdal.Open(gxf_file)
    gridinfo = gxf_grid.GetGeoTransform()
    grid_topleftX = gridinfo[0]  # xmin
    grid_topleftY = gridinfo[3]  # ymax
    pixelwidth = gridinfo[1]
    pixelheight = gridinfo[5]
    cols = gxf_grid.RasterXSize
    rows = gxf_grid.RasterYSize
    grid_lowerrightX = grid_topleftX + (cols * pixelwidth)  # xmax
    grid_lowerrightY = grid_topleftY + (rows * pixelheight)  # ymin

    xmin = grid_topleftX
    xmax = grid_lowerrightX
    ymin = grid_lowerrightY
    ymax = grid_topleftY

    # create a numpy array from grid
    blank = -1.e+12
    mag = np.array(gxf_grid.GetRasterBand(1).ReadAsArray())
    mag[mag == blank] = np.nan

    # plot
    fig, ax = plt.subplots(figsize=(10, 10))
    ax = plt.subplot(aspect='equal')

    norm_eq = cm_utils.MidpointNormalize(vmin=np.nanpercentile(mag, 5),
                                         vmax=np.nanpercentile(mag, 95),
                                         midpoint=0)

    cmap = sci_cm.roma

    ax.imshow(mag, extent=(xmin, xmax, ymin, ymax), norm=norm_eq, alpha=0.8)

    plt.xlabel('Easting')
    plt.ylabel('Northing')

    colors = np.random.rand(len(center))
    cmap = plt.cm.tab20
    c = cmap(colors)

    for win, x, y, cl in zip(np.array(center_coords)[:, 0],
                             np.array(center_coords)[:, 1],
                             np.array(center_coords)[:, 2],
                             np.array(c)):
        ax.add_patch(Rectangle(xy=(x-window_width/2, y-window_width/2),
                               width=window_width-1000,
                               height=window_width-1000,
                               linewidth=2, color=cl, fill=False, alpha=0.5))
        ax.scatter(x, y, color=cl)
        ax.annotate(str(int(win)), (x, y))

    plt.savefig(out_dir + dsep + 'window_' + str(window_width/1000) + dsep +
                'windows_centers.png', dpi=300)


def splitgrid(working_dir, gridfilename, window_width, overlap_factor):
    """
    Split the main grid into smaller grids.

    Expected input format is "surfer 7"

    This function takes in the large initial grid and splits it into windows
    of specified width and overlap.

    It creates various output directories and makes a plot of each grid subset
    which can be used for checking for empty cells before running the fft.

    update the locations of gdalTranslate and gdalinfo if required

    update grid blank (nan) value if required.

    Parameters
    ----------
    working_dir: main directory where the main data grid lives and where
                 processing subdirectories will be created.

    gridfilename: main data grid to be split.  Expected format is Surfer 7.

    window_width:  sub grid width in metres (ie 100 km = 100000.)

    overlap_factor: overlap of adjacent grids
    use factor of 2. with 100000. window_width to get 50 km spaced points
    use factor of 4. with 200000. window_width to get 50 km spaced points
    use factor of 6. with 300000. window_width to get 50 km spaced points

    Returns
    -------
    Many subgrids in a sub directory "subsetted grids" under the working_dir.

    """
    plt.close('all')
    # set no data value
    blank = 1.701410009187828e+038

    # define window width
    window_width = window_width  # in grid units (degrees or metres)

    window_step = window_width/overlap_factor

    dsep = '/'

    main_dir = working_dir

    grid = (main_dir + dsep + gridfilename)

    # prefix for output grids
    out_dir = main_dir + dsep + 'subsetted_grids'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    if not os.path.exists(out_dir + dsep + 'window_' + str(window_width/1000)):
        os.makedirs(out_dir + dsep + 'window_' + str(window_width/1000), exist_ok=True)

    if not os.path.exists(out_dir + dsep + 'fft_output' +
                          str(window_width/1000)):
        os.makedirs(main_dir + dsep + 'fft_output'
                    + dsep + 'window_' + str(window_width/1000), exist_ok=True)

    if not os.path.exists(out_dir + dsep + 'spectra_plots' +
                          str(window_width/1000)):
        os.makedirs(main_dir + dsep + 'spectra_plots'
                    + dsep + 'window_' + str(window_width/1000), exist_ok=True)

    if not os.path.exists(out_dir + dsep + 'results' + str(window_width/1000)):
        os.makedirs(main_dir + dsep + 'results'
                    + dsep + 'window_' + str(window_width/1000), exist_ok=True)

    outfile = out_dir + dsep + 'window_' + str(window_width/1000) + dsep + \
        'window'

    gdal.Info(str(grid))

    # get grid extents
    datagrid = gdal.Open(grid)
    gridinfo = datagrid.GetGeoTransform()
    grid_topleftX = gridinfo[0]  # xmin
    grid_topleftY = gridinfo[3]  # ymax
    pixelwidth = gridinfo[1]
    pixelheight = gridinfo[5]

    cols = datagrid.RasterXSize
    rows = datagrid.RasterYSize

    grid_lowerrightX = grid_topleftX + (cols * pixelwidth)  # xmax
    grid_lowerrightY = grid_topleftY + (rows * pixelheight)  # ymin

    xmin = grid_topleftX
    xmax = grid_lowerrightX
    ymin = grid_lowerrightY
    ymax = grid_topleftY

    print("Grid top left X: " + str(grid_topleftX))
    print("Grid top left Y: " + str(grid_topleftY))
    print("Grid lower right X: " + str(grid_lowerrightX))
    print("Grid lower right Y: " + str(grid_lowerrightY))
    print("Grid width : " + str(grid_lowerrightX-grid_topleftX) + " m")
    print("Grid height : " + str(grid_topleftY-grid_lowerrightY) + " m")
    print("xmin : " + str(xmin) + "m")
    print("xmax : " + str(xmax) + "m")
    print("ymin : " + str(ymin) + "m")
    print("ymax : " + str(ymax) + "m")
    print("\n")

    # setup windowing of the main grid
    # define initial window corners in the top left corner of the grid
    window_ulx = grid_topleftX
    window_uly = grid_topleftY
    window_lrx = window_ulx + window_width
    window_lry = window_uly - window_width

    # set subset counter
    subset = 1

    # process subsets loop moving from left to right and top to bottom
    while window_uly > grid_lowerrightY:
        while window_ulx < grid_lowerrightX:
            print('processing subset #' + str(subset))
            gdal.Translate(outfile + '_' + str(subset) + '.grd', str(grid),
                           format='GS7BG', projWin=[window_ulx, window_uly,
                                                    window_lrx, window_lry]
                           )

            # update subset coords - move in +X direction
            window_ulx = window_ulx + window_step
            window_lrx = window_lrx + window_step
            time.sleep(1)
            subset += 1

        # reset x coords to left hand side of grid
        window_ulx = grid_topleftX
        window_lrx = window_ulx + window_width
        # increment y coords down a step
        window_uly = window_uly - window_step
        window_lry = window_lry - window_step

    # loop through created subgrids and make a file of the grid centers
    os.chdir(out_dir + dsep + 'window_' + str(window_width/1000))
    centerx = []
    centery = []
    window_num = []

    for fn in os.listdir(out_dir + dsep + 'window_' + str(window_width/1000)):
        if fn.endswith(".grd"):
            print("file:", fn)
            datagrid = gdal.Open(fn)
            gridinfo = datagrid.GetGeoTransform()
            centerx.append(gridinfo[0] + window_width/2)
            centery.append(gridinfo[3] - window_width/2)
            window_num.append(os.path.splitext(fn)[0])
            time.sleep(1)

            # add making a plot of each grid here
            # get data grid extents

            grid_topleftX = gridinfo[0]  # xmin
            grid_topleftY = gridinfo[3]  # ymax
            pixelwidth = gridinfo[1]
            pixelheight = gridinfo[5]
            cols = datagrid.RasterXSize
            rows = datagrid.RasterYSize
            grid_lowerrightX = grid_topleftX + (cols * pixelwidth)  # xmax
            grid_lowerrightY = grid_topleftY + (rows * pixelheight)  # ymin

            xmin = grid_topleftX
            xmax = grid_lowerrightX
            ymin = grid_lowerrightY
            ymax = grid_topleftY

            # create a numpy array from grid
            mag = np.array(datagrid.GetRasterBand(1).ReadAsArray())
            mag[mag == blank] = np.nan

            # test if array contains all nans (or no data)
            test = mag[np.logical_not(np.isnan(mag))]

            if len(test) != len(mag.flatten()):
                print("skipping, grid contains nans")

            else:
                # plot
                fig, ax = plt.subplots(figsize=(5, 5))
                mag_equalised = equalize_hist(mag)
                ls = LightSource(azdeg=45, altdeg=45)
                norm_eq = cm_utils.MidpointNormalize(vmin=mag_equalised.min(),
                                                     vmax=mag_equalised.max(),
                                                     midpoint=0)

                cmap = sci_cm.vik
                relief = ls.shade(mag_equalised, cmap=cmap,
                                  blend_mode='overlay',
                                  vert_exag=50, vmin=0, vmax=1)

                ax.imshow(relief, extent=(xmin, xmax, ymin, ymax),
                          norm=norm_eq, alpha=0.8)

                plt.xlabel('Easting')
                plt.ylabel('Northing')
                plt.savefig(os.path.splitext(fn)[0] + '.png', dpi=100)
                plt.close()

    centerx = np.array(centerx)
    centery = np.array(centery)
    window_num = np.array(window_num)
    center_coords = np.c_[window_num, centerx, centery]
    np.savetxt('window_' + str(window_width/1000) + '_km_center_coords.csv',
               center_coords, fmt='%s', delimiter=',')

    # Plot center points and subset outlines overlain on maggrid
    # read in grid
    datagrid = gdal.Open(grid)
    gridinfo = datagrid.GetGeoTransform()
    grid_topleftX = gridinfo[0]  # xmin
    grid_topleftY = gridinfo[3]  # ymax
    pixelwidth = gridinfo[1]
    pixelheight = gridinfo[5]
    cols = datagrid.RasterXSize
    rows = datagrid.RasterYSize
    grid_lowerrightX = grid_topleftX + (cols * pixelwidth)  # xmax
    grid_lowerrightY = grid_topleftY + (rows * pixelheight)  # ymin

    xmin = grid_topleftX
    xmax = grid_lowerrightX
    ymin = grid_lowerrightY
    ymax = grid_topleftY

    # create a numpy array from grid
    mag = np.array(datagrid.GetRasterBand(1).ReadAsArray())
    mag[mag == blank] = np.nan

    # plot
    fig, ax = plt.subplots(figsize=(10, 10))
    ax = plt.subplot(aspect='equal')

    mag_equalised = equalize_hist(mag)
    ls = LightSource(azdeg=45, altdeg=45)
    norm_eq = cm_utils.MidpointNormalize(vmin=mag_equalised.min(),
                                         vmax=mag_equalised.max(),
                                         midpoint=0)

    cmap = sci_cm.vik
    relief = ls.shade(mag_equalised, cmap=cmap, blend_mode='overlay',
                      vert_exag=50, vmin=0, vmax=1)

    ax.imshow(relief, extent=(xmin, xmax, ymin, ymax), norm=norm_eq, alpha=0.8)

    plt.xlabel('Easting')
    plt.ylabel('Northing')

    colors = np.random.rand(len(centerx))
    cmap = plt.cm.tab20
    c = cmap(colors)

    for x, y, cl in zip(centerx, centery, c):
        ax.add_patch(Rectangle(xy=(x-window_width/2, y-window_width/2),
                               width=window_width-1000,
                               height=window_width-1000,
                               linewidth=2, color=cl, fill=False))
        ax.scatter(x, y, color=cl)

    plt.savefig(out_dir + dsep + 'window_' + str(window_width/1000) + dsep +
                'windows_centers.png', dpi=300)


def runfft(working_dir, window_dir):
    """
    Run the FFT on each grid.

    This function runs the fft on each of the grids created from splitgrid.
    It does not produce a file if the grid contains ANY nans.

    Uses pygmt to call grdfft

    Parameters
    ----------
    working_dir: head directory where the main grid file is located
    window_dir: name of the window_dir created when the grids were split

    Returns
    -------
    file: *.txt file with FFT output
    """
    main_dir = working_dir
    dsep = '/'
    griddir = (main_dir + dsep + 'subsetted_grids' + dsep + window_dir)
    outdir = main_dir + dsep + 'fft_output' + dsep + window_dir
    os.chdir(griddir)
    print(griddir)

    blank = 1.701410009187828e+038

    # run fft loop through all the grd files in the dir
    for fn in os.listdir(griddir):
        if fn.endswith(".grd"):
            print("Doing fft on file:", fn)
            datagrid = gdal.Open(fn)
            mag = np.array(datagrid.GetRasterBand(1).ReadAsArray())
            mag[mag == blank] = np.nan

            # test if array contains all nans (or no data)
            test = mag[np.logical_not(np.isnan(mag))]
            if len(test) != len(mag.flatten()):
                print("skipping, grid contains nans")
            else:
                with pygmt.clib.Session() as session:
                    fout = outdir + '/' + os.path.splitext(fn)[0] + str('.txt')
                    args = f"{fn} -Nq+n+a -Ern -G{fout}" # test with -Ern to normalise
                    session.call_module('grdfft', args)

                    # Plot raw spectrum
                    data = pd.read_csv(fout, names=['Wavenumber', 'Power',
                                                    'Stdev'],
                                       delim_whitespace=True)

                    fig, ax = plt.subplots(figsize=(4, 4))

                    ax.scatter(data.Wavenumber*1000, data.Power, marker='o',
                               s=2, label="Data")
                    ax.set_xlim(0, data.Wavenumber.max()*1000)
                    ax.set_xlabel(r'$|k| [km^{-1}$]')
                    ax.set_ylabel('Power')
                    plt.title(fn.split('.')[0])
                    ax.semilogy()
                    ax.set_ylim(data.Power.min(), data.Power.max())
                    plt.tight_layout()
                    plt.savefig(outdir + '/' + os.path.splitext(fn)[0] +
                                str('_spectrum.png'))


def runfft_geosoft(working_dir, window_dir):
    """
    Run the FFT on each grid using the geosoft fft. This requires a geosoft
    licence.

    This function runs the fft on each of the grids created from splitgrid.

    Parameters
    ----------
    working_dir: head directory where the main grid file is located
    window_dir: name of the window_dir created when the grids were split

    Returns
    -------
    file: *.txt file with FFT output
    """
    main_dir = working_dir
    dsep = '/'
    griddir = (main_dir + dsep + 'subsetted_grids' + dsep + window_dir)
    outdir = main_dir + dsep + 'fft_output' + dsep + window_dir
    os.chdir(griddir)
    print(griddir)
    gx.gxc = gx.gx.GXpy(log=print)

    for fn in os.listdir(griddir):
        if fn.endswith(".grd"):
            print("\n Doing fft on file:", fn)

            # open geosoft grid
            grd = gx.grid.Grid.open(fn)

            # test for empty grids and skip if encountered.
            if grd.statistics().get('min') is not None:
                fft = gx.grid_fft.GridFFT(grd)

                # this returns a numpy array shaped (n, 5), where each element
                # contains the wavenumber, number_samples, log(average_power),
                # depth_3 and depth_5)
                # normalise source spectrum by subtracting average power density
                source_spectrum = fft.radially_averaged_spectrum()
                source_spectrum[:, 2] = source_spectrum[:, 2] - fft.log_average_spectral_density()

                # the power is relative to the average density
                print(f'average spectral density: {fft.log_average_spectral_density()}')

                # save spectrum to file
                fout = outdir + '/' + os.path.splitext(fn)[0] + str('.txt')
                np.savetxt(fout, source_spectrum)

                fig, ax = plt.subplots(figsize=(8, 5))
                ax.plot(source_spectrum[:, gx.grid_fft.I_WAVENUMBER],
                        source_spectrum[:, gx.grid_fft.I_LOG_POWER],
                        label='source', marker='o', ms=1, lw=1)
                ax.set_xlabel('wavenumber [cycles/km]')
                ax.set_ylabel('Ln Power')
                plt.title(fn.split('.')[0])
                plt.legend()
                plt.grid()
                plt.savefig(outdir + '/' + os.path.splitext(fn)[0] +
                            str('_spectrum.png'))
                plt.close()
            else:
                print('skipping fft, grid is blank')


def calccpd(working_dir, window_dir, Zo_start, Zo_end, Zt_start, Zt_end, Beta,
            K=2.5, CT=570):
    """
    Calculate the Curie Point Depth.

    Also calculates an error for each using the stdev value from the fft
    Uses the Tanaka method to calculate CPD depths (no fractal correction).

    Read in the output from the GMT grdfft.  Use the output in wavenumbers not
    wavelengths.

    Loops through each file in directory and then incorporates the centre
    coords of each window to the results.

    Incorporates the fractal correction "Beta" from Bansal et al 2011,
    Geophysics vol 76

    Also calculates heatflow in mW/m^2

    Parameters
    ----------
    working_dir: main directory where the main data grid lives and where
                 processing subdirectories will be created
    window_dir:
    Zo_start: starting wavenumnber for centroid depth
    Zo_end: end wavenumnber for centroid depth
    Zt_start: starting wavenumnber for top depth
    Zt_end: end wavenumnber for top depth

    Z* parameters are the wavenumber values used to fit the straight line to.
    May take some experimenting to find appropriate values.

    Beta: fractal parameter.  Beta <= 2. Beta = 0 is the same as no correction

    #Thermal parameters
    K = 2.5 # W/m*K
    CT = 570 # 580C Curie Temp - 10C surface temp

    Returns
    -------
    file: "cpd_results_beta_*.csv' in the /results directory
    plots: Spectra plots are written in the /spectra_plots directory
    """
    plt.close('all')

    dsep = '/'

    main_dir = working_dir

    griddir = main_dir + dsep + 'subsetted_grids/' + window_dir
    datadir = main_dir + dsep + 'fft_output/' + window_dir
    plotdir = main_dir + dsep + 'spectra_plots/' + window_dir
    resultdir = main_dir + dsep + 'results/' + window_dir

    gridcenters = pd.read_csv(griddir + dsep + window_dir +
                              '_km_center_coords.csv', comment='#',
                              names=['window', 'east', 'north'],
                              index_col='window')

    # Tables to store results in
    centroid_depth_table = []
    centroid_depth_err_table = []
    centroid_rsq_table = []
    top_depth_table = []
    top_depth_err_table = []
    top_rsq_table = []
    base_depth_table = []
    base_depth_err_table = []
    window_table = []
    heatflow_table = []
    gradient_table = []

    # This loops through all the files spectra files in the directory
    # and writes the results to file

    for fn in os.listdir(datadir):
        if fn.endswith(".txt"):
            print("\nMaking plot of file:", fn)
            plotfile = plotdir + dsep + fn.split('.')[0] + '.png'
            window_table.append(fn.split('.')[0])
            data = pd.read_csv(datadir + dsep + fn, delim_whitespace=True,
                               names=['Wavenumber', 'Power', 'Stdev'],
                               comment='#')

            print(data.columns)
            # Calculate centroid depth, Zo
            # centroid fit points (k*1000/2pi), may need to adjust values
            # slightly if index errors, they need to match actual values
            Zo_start = Zo_start
            Zo_end = Zo_end

            assert (Zo_start < Zo_end), 'Zo_start value must be smaller than Zo_end value'

            # calculate centroid
            # this is the y axis value
            data['y_centroid'] = np.log(data.Wavenumber**(Beta) *
                                        (data.Power / data.Wavenumber**2)) / (2*np.pi)

            # this is the x axis value (wavenumber/2pi) in cycles per km
            data['wavenumber_2pi'] = np.abs(data.Wavenumber * 1000) / (2*np.pi)

            # select the data to fit
            centroid_fit_start = data[data.wavenumber_2pi.round(3) == Zo_start].index[0]
            centroid_fit_end = data[data.wavenumber_2pi.round(3) == Zo_end].index[0]
            centroid_fit_pts = slice(centroid_fit_start, centroid_fit_end)

            centroid_fit_wavenumber = np.abs(data.wavenumber_2pi[centroid_fit_pts])
            centroid_fit_power = data.y_centroid[centroid_fit_pts]

            # fit the data
            centroid_fit_wavenumber = sm.add_constant(centroid_fit_wavenumber)
            cent_model = sm.OLS(centroid_fit_power, centroid_fit_wavenumber)
            cent_results = cent_model.fit()

            centroid_depth = cent_results.params[1]
            centroid_depth_table.append(centroid_depth)
            centroid_depth_err = cent_results.bse[1]
            centroid_depth_err_table.append(centroid_depth_err)
            centroid_rsq = cent_results.rsquared
            centroid_rsq_table.append(centroid_rsq)

            print(r"centroid depth = %0.1f \u00B1 %0.1f km" % (centroid_depth,
                                                               centroid_depth_err))

            # Calculate top depth Zt
            Zt_start = Zt_start
            Zt_end = Zt_end
            assert (Zt_start < Zt_end), 'Zt_start value must be smaller than Zt_end value'

            # calculate top depth
            data['y_top'] = np.log(data.Wavenumber**(Beta) * data.Power) / (2*np.pi)
            data['y_top_std'] = np.log(data.Wavenumber**(Beta) * data.Stdev) / (2*np.pi)

            # select  the data to fit
            top_fit_start = data[data.wavenumber_2pi.round(3) == Zt_start].index[0]
            top_fit_end = data[data.wavenumber_2pi.round(3) == Zt_end].index[0]
            top_fit_pts = slice(top_fit_start, top_fit_end)

            top_fit_wavenumber = np.abs(data.wavenumber_2pi[top_fit_pts])
            top_fit_power = data.y_top[top_fit_pts]

            # fit the data
            top_fit_wavenumber = sm.add_constant(top_fit_wavenumber)
            top_model = sm.OLS(top_fit_power, top_fit_wavenumber)
            top_results = top_model.fit()

            top_depth = top_results.params[1]
            top_depth_table.append(top_depth)
            top_depth_err = top_results.bse[1]
            top_depth_err_table.append(top_depth_err)
            top_rsq = top_results.rsquared
            top_rsq_table.append(top_rsq)

            print(r'top depth = %0.1f $\pm$ %0.1f km' % (top_depth,
                                                         top_depth_err))
            print(r'centroid depth = %0.1f $\pm$ %0.1f km' % (centroid_depth,
                                                              centroid_depth_err))

            # CALCULATE BASE DEPTH
            base_depth = 2 * centroid_depth - top_depth
            base_depth_table.append(base_depth)
            base_depth_err = np.sqrt(np.sum(top_depth_err**2 + centroid_depth_err**2))
            base_depth_err_table.append(base_depth_err)
            print('base depth = %0.1f ' % base_depth)

            # CALCULATE GEOTHERMAL GRADIENT
            gradient = CT/base_depth * -1
            gradient_table.append(gradient)

            # CALCULATE HEATFLOW
            heatflow = gradient * K
            heatflow_table.append(heatflow)

            # MAKE PLOTS
            plt.figure(figsize=(12, 4))

            # Plot raw spectrum
            plt.subplot(131)
            plt.scatter(data.Wavenumber*1000, data.Power, marker='o', s=2,
                        label="Data")
            plt.xlim(0, data.Wavenumber.max()*1000)
            plt.xlabel(r'$|k| [km^{-1}$]')
            plt.ylabel('Power')
            plt.title(fn.split('.')[0])
            plt.semilogy()
            plt.ylim(data.Power.min(), data.Power.max())

            # plot centroid
            plt.subplot(132)
            plt.scatter(data.wavenumber_2pi, data.y_centroid, marker='o', s=2,
                        label="Data")

            plt.plot(centroid_fit_wavenumber, cent_results.fittedvalues, 'r-',
                     label="Linear fit")
            plt.xlim(0, data.wavenumber_2pi.max())
            plt.xlabel(r'$|k|/(2 \pi)$ [km$^{-1}$]')
            plt.ylabel(r'$\ln(\Phi_{\Delta T}(|k|)^{1/2}/|k|)/(2\pi)$')
            plt.legend(loc="best")
            plt.title('Centroid')
            plt.annotate(r'centroid depth %0.1f $\pm$ %0.1f km'
                         % (centroid_depth, centroid_depth_err),
                         xy=(0.01, 1.5), xytext=(0.01, 1.5))

            # plot top
            plt.subplot(133)
            plt.scatter(data.wavenumber_2pi, data.y_top, marker='o', s=2,
                        label="Data")
            plt.plot(top_fit_wavenumber, top_results.fittedvalues, 'r-',
                     label="Linear fit")
            plt.xlim(0, data.wavenumber_2pi.max())
            plt.xlabel(r'$|k|/(2 \pi)$ [km$^{-1}$]')
            plt.ylabel(r'$\ln(\Phi_{\Delta T}(|k|)^{1/2})/(2\pi)$')
            plt.legend(loc="best")
            plt.title('Top')
            plt.annotate(r'top depth %0.1f $\pm$ %0.1f km'
                         % (top_depth, top_depth_err),
                         xy=(0.01, 0.2), xytext=(0.01, 0.2))
            plt.annotate(r'base depth %0.1f $\pm$ %0.1f km'
                         % (base_depth, base_depth_err),
                         xy=(0.01, 0.3), xytext=(0.01, 0.3))
            plt.savefig(plotfile, dpi=300, bbox_inches='tight')

    # create dataframe of results tables
    results = pd.DataFrame(centroid_depth_table, index=window_table,
                           columns=['centroid_depth_table'])
    results['centroid_depth_err_table'] = centroid_depth_err_table
    results['centroid_rsq_table'] = centroid_rsq_table
    results['top_depth_table'] = top_depth_table
    results['top_depth_err_table'] = top_depth_err_table
    results['top_rsq_table'] = top_rsq_table
    results['base_depth_table'] = base_depth_table
    results['base_depth_err_table'] = base_depth_err_table
    results['gradient_table'] = gradient_table
    results['heatflow_table'] = heatflow_table

    # combine cell center file with results table
    result = pd.concat([gridcenters, results], axis=1).dropna()
    result['window'] = result.index
    result.to_csv(resultdir + dsep + 'cpd_results_beta_%.0f' % Beta + '.csv',
                  header=['window', 'east', 'north', 'centroid_depth',
                          'cent_depth_err', 'cent_rsq', 'top_depth',
                          'top_depth_err', 'top_rsq', 'base_depth',
                          'base_depth_err', 'gradient', 'heatflow'],
                  columns=['window', 'east', 'north', 'centroid_depth_table',
                           'centroid_depth_err_table', 'centroid_rsq_table',
                           'top_depth_table', 'top_depth_err_table',
                           'top_rsq_table', 'base_depth_table',
                           'base_depth_err_table', 'gradient_table',
                           'heatflow_table'],
                  float_format='%.2f', index=False)

    # Plot curie depth vs heatflow
    plt.figure()
    plt.scatter(results.heatflow_table, results.base_depth_table)
    plt.xlabel('heatflow (mW/m$^2$)')
    plt.ylabel('Curie Pt depth (km)')
    plt.savefig(resultdir + dsep + 'heatflow_CPD.png', dpi=300)


def calccpd_geosoft(working_dir, window_dir, Zo_start, Zo_end, Zt_start,
                    Zt_end, Beta, K=2.5, CT=570):
    """
    Calculate the Curie Point Depth.

    Read in the output from the geosoft raidal powerspectrum.
    Use the output in wavenumbers not wavelengths.

    Loops through each file in directory and then incorporates the centre
    coords of each window to the results.

    Incorporates the fractal correction "Beta" from Bansal et al 2011,
    Geophysics vol 76

    Also calculates heatflow in mW/m^2

    Parameters
    ----------
    working_dir: main directory where the main data grid lives and where
                 processing subdirectories will be created
    window_dir:
    Zo_start (int): starting wavenumber index for centroid depth
    Zo_end (int): end wavenumber index for centroid depth
    Zt_start (int): starting wavenumber index for top depth
    Zt_end (int): end wavenumber index for top depth

    Z* parameters are the wavenumber index values used to fit the straight line.
    They start from 0 at lowest wavenumber.

    Beta: fractal parameter.  Beta <= 2. Beta = 0 is the same as no correction

    #Thermal parameters
    K = 2.5 # W/m*K
    CT = 570 # 580C Curie Temp - 10C surface temp

    Returns
    -------
    file: "cpd_results_beta_*.csv' in the /results directory
    plots: Spectra plots are written in the /spectra_plots directory
    """
    plt.close('all')

    dsep = '/'

    main_dir = working_dir

    griddir = main_dir + dsep + 'subsetted_grids/' + window_dir
    datadir = main_dir + dsep + 'fft_output/' + window_dir
    plotdir = main_dir + dsep + 'spectra_plots/' + window_dir
    resultdir = main_dir + dsep + 'results/' + window_dir

    gridcenters = pd.read_csv(griddir + dsep + window_dir +
                              '_km_center_coords.csv', comment='#',
                              names=['window', 'east', 'north'],
                              index_col='window', delim_whitespace=True)
    gridcenters.index = str('window_') + gridcenters.index.astype(str)

    # Tables to store results in
    centroid_depth_table = []
    centroid_depth_err_table = []
    centroid_rsq_table = []
    top_depth_table = []
    top_depth_err_table = []
    top_rsq_table = []
    base_depth_table = []
    base_depth_err_table = []
    window_table = []
    heatflow_table = []
    gradient_table = []

    # This loops through all the files spectra files in the directory
    # and writes the results to file

    for fn in os.listdir(datadir):
        if fn.endswith(".txt"):
            print("\nMaking plot of file:", fn)
            plotfile = plotdir + dsep + fn.split('.')[0] + '.png'
            window_table.append(fn.split('.')[0])

            data = pd.read_csv(datadir + dsep + fn, delim_whitespace=True,
                               names=['Wavenumber', 'Num_samples', 'Ln_P',
                                      'Depth_3', 'Depth_5'], comment='/')

            # extract the log average power density value and add back. Geosoft
            # subtracts the log average power density from the spectrum.
            '''
            ave_power = pd.read_csv(datadir + dsep + fn,
                                    nrows=1, header=None, skiprows=[0, 1],
                                    delim_whitespace=True)
            print(ave_power)
            data['Ln_raw'] = data.Ln_P + ave_power.values
            '''
            data['Power'] = np.exp(data.Ln_P)

            # Calculate centroid depth, Zo
            Zo_start = int(Zo_start)
            Zo_end = int(Zo_end)

            assert (Zo_start < Zo_end), 'Zo_start value must be smaller than Zo_end value'

            # Calculate centroid. This is the y axis value.

            # this is the y axis value Eqn7 Bansal
            # data['y_centroid'] = np.log(data.Wavenumber**(Beta) * (data.Power / data.Wavenumber**2))/(2*np.pi)

            # Eqn 6 # Andres 2017 JGR
            # Eqn5 Wang and Li 2015
            # Beta -1 converts 3D to 2D
            data['y_centroid'] = np.log(data.Wavenumber**(Beta-1)/2 * data.Power / data.Wavenumber) / (2*np.pi)

            # this is the x axis value (wavenumber/2pi) in cycles per km
            data['wavenumber_2pi'] = np.abs(data.Wavenumber)/(2*np.pi)

            # select the data to fit
            centroid_fit_pts = slice(Zo_start, Zo_end)

            centroid_fit_wavenumber = np.abs(data.wavenumber_2pi[centroid_fit_pts])
            centroid_fit_power = data.y_centroid[centroid_fit_pts]

            # fit the data
            centroid_fit_wavenumber = sm.add_constant(centroid_fit_wavenumber)
            cent_model = sm.OLS(centroid_fit_power, centroid_fit_wavenumber)
            cent_results = cent_model.fit()

            centroid_depth = cent_results.params[1]
            centroid_depth_table.append(centroid_depth)
            centroid_depth_err = cent_results.bse[1]
            centroid_depth_err_table.append(centroid_depth_err)
            centroid_rsq = cent_results.rsquared
            centroid_rsq_table.append(centroid_rsq)

            print(r"centroid depth = %0.1f +\- %0.1f km" % (centroid_depth,
                                                            centroid_depth_err))

            # Calculate top depth Zt
            Zt_start = int(Zt_start)
            Zt_end = int(Zt_end)
            assert (Zt_start < Zt_end), 'Zt_start value must be smaller than Zt_end value'

            # calculate top depth
            # data['y_top'] = np.log(data.Wavenumber**Beta * data.Power) / (2*np.pi)

            # Eqn5 Andres 2017 JGR and Eqn4 Wang and Li 2015
            data['y_top'] = np.log(data.Wavenumber**(Beta-1)/2 * data.Power) / (2*np.pi)

            # select  the data to fit
            top_fit_pts = slice(Zt_start, Zt_end)

            top_fit_wavenumber = np.abs(data.wavenumber_2pi[top_fit_pts])
            top_fit_power = data.y_top[top_fit_pts]

            # fit the data
            top_fit_wavenumber = sm.add_constant(top_fit_wavenumber)
            top_model = sm.OLS(top_fit_power, top_fit_wavenumber)
            top_results = top_model.fit()

            top_depth = top_results.params[1]
            top_depth_table.append(top_depth)
            top_depth_err = top_results.bse[1]
            top_depth_err_table.append(top_depth_err)
            top_rsq = top_results.rsquared
            top_rsq_table.append(top_rsq)

            print(r'top depth = %0.1f +\- %0.1f km' % (top_depth,
                                                       top_depth_err))
            print(r'centroid depth = %0.1f +\- %0.1f km' % (centroid_depth,
                                                            centroid_depth_err))

            # CALCULATE BASE DEPTH
            base_depth = 2 * centroid_depth - top_depth
            base_depth_table.append(base_depth)
            base_depth_err = np.sqrt(np.sum(top_depth_err**2 + centroid_depth_err**2))
            base_depth_err_table.append(base_depth_err)
            print('base depth = %0.1f ' % base_depth)

            # CALCULATE GEOTHERMAL GRADIENT
            gradient = CT/base_depth * -1
            gradient_table.append(gradient)

            # CALCULATE HEATFLOW
            heatflow = gradient * K
            heatflow_table.append(heatflow)

            # MAKE PLOTS
            plt.figure(figsize=(12, 4))

            # Plot raw spectrum
            plt.subplot(131)
            plt.plot(data.Wavenumber, data.Ln_P, marker='o', ms=2, lw=1,
                     label='data')
            plt.xlim(-0.5, data.Wavenumber.max())
            plt.xlabel(r'$|k| [km^{-1}$]')
            plt.ylabel('Ln Power')
            plt.title(fn.split('.')[0])
            plt.ylim(data.Ln_P.min(), data.Ln_P.max())

            # plot centroid
            plt.subplot(132)
            plt.plot(data.wavenumber_2pi, data.y_centroid, marker='o', ms=2,
                     lw=1, label="Data")
            plt.plot(centroid_fit_wavenumber, cent_results.fittedvalues, 'r-',
                     lw=1, label="Linear fit")
            plt.xlim(-0.05, data.wavenumber_2pi.max())
            plt.xlabel(r'$|k|/(2 \pi)$ [km$^{-1}$]')
            plt.ylabel(r'$\ln(\Phi_{\Delta T}(|k|)^{(\beta -1)/2}/|k|)/(2\pi)$')
            plt.legend(loc="best")
            plt.title('Centroid')
            plt.annotate(r'centroid depth %0.1f $\pm$ %0.1f km'
                         % (centroid_depth, centroid_depth_err),
                         xy=(0.01, 0), xytext=(0.01, 0))

            # plot top
            plt.subplot(133)
            plt.plot(data.wavenumber_2pi, data.y_top, marker='o', ms=2, lw=1,
                     label="Data")
            plt.plot(top_fit_wavenumber, top_results.fittedvalues, 'r-', lw=1,
                     label="Linear fit")
            plt.xlim(-0.05, data.wavenumber_2pi.max())
            plt.xlabel(r'$|k|/(2 \pi)$ [km$^{-1}$]')
            plt.ylabel(r'$\ln(\Phi_{\Delta T}(|k|)^{(\beta -1)/2})/(2\pi)$')
            plt.legend(loc="best")
            plt.title('Top')
            plt.annotate(r'top depth %0.1f $\pm$ %0.1f km'
                         % (top_depth, top_depth_err),
                         xy=(0.1, 0.1), xytext=(0.1, 0.1))
            plt.annotate(r'base depth %0.1f $\pm$ %0.1f km'
                         % (base_depth, base_depth_err),
                         xy=(0.1, 0), xytext=(0.1, 0))
            # plt.tight_layout()
            plt.savefig(plotfile, dpi=300, bbox_inches='tight')
            plt.close()

    # create dataframe of results tables
    results = pd.DataFrame(centroid_depth_table, index=window_table,
                           columns=['centroid_depth_table'])
    results['centroid_depth_err_table'] = centroid_depth_err_table
    results['centroid_rsq_table'] = centroid_rsq_table
    results['top_depth_table'] = top_depth_table
    results['top_depth_err_table'] = top_depth_err_table
    results['top_rsq_table'] = top_rsq_table
    results['base_depth_table'] = base_depth_table
    results['base_depth_err_table'] = base_depth_err_table
    results['gradient_table'] = gradient_table
    results['heatflow_table'] = heatflow_table

    # combine cell center file with results table
    result = pd.concat([gridcenters, results], axis=1).dropna()
    result['window'] = result.index

    result.to_csv(resultdir + dsep + 'cpd_results_beta_%.0f' % Beta + '.csv',
                  header=['window', 'east', 'north', 'centroid_depth',
                          'cent_depth_err', 'cent_rsq', 'top_depth',
                          'top_depth_err', 'top_rsq', 'base_depth',
                          'base_depth_err', 'gradient', 'heatflow'],
                  columns=['window', 'east', 'north', 'centroid_depth_table',
                           'centroid_depth_err_table', 'centroid_rsq_table',
                           'top_depth_table', 'top_depth_err_table',
                           'top_rsq_table', 'base_depth_table',
                           'base_depth_err_table', 'gradient_table',
                           'heatflow_table'],
                  float_format='%.2f', index=False)

    # Plot curie depth vs heatflow
    plt.figure()
    plt.scatter(results.heatflow_table, results.base_depth_table)
    plt.xlabel('heatflow (mW/m$^2$)')
    plt.ylabel('Curie Pt depth (km)')
    plt.savefig(resultdir + dsep + 'heatflow_CPD.png', dpi=300)

