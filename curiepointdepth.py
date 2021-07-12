# -*- coding: utf-8 -*-
"""
A set of functions to calculate curie point depths from an aeromagnetic\
data grid using the Tanaka method (Tanaka 1999, Tectonophysics).

@author: craigm
Created on Wed May 04 11:09:00 2016

"""

import subprocess as sp
from osgeo import gdal
import os
import time
import numpy as np
import matplotlib.pyplot as plt
import graphics
import pandas as pd
import statsmodels.api as sm
import cm_utils


def splitgrid(working_dir, gridfilename, window_width, overlap_factor):
    """
    Split the main grid into smaller grids.
    
    Expected input format is "surfer 7"

    Requires "graphics" from:
    https://github.com/jobar8/graphics

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

    window_width:  sub grid width in metres (ie 100 km = 100000)

    overlap_factor: overlap of adjacent grids
    use factor of 2 with 100000 window_width to get 50 km spaced points
    use factor of 4 with 200000 window_width to get 50 km spaced points
    use factor of 6 with 300000 window_width to get 50 km spaced points

    Returns
    -------
    Many subgrids in a sub directory "subsetted grids" under the working_dir.

    """
    plt.close('all')
    # set no data value
    blank = 1.701410009187828e+038
    # set gdal locations
    gdalTranslate = r'C:/Program Files/GDAL/gdal_translate.exe'
    gdalinfo = r'C:/Program Files/GDAL/gdalinfo.exe'

    # define window width
    window_width = window_width  # in grid units (degrees or metres)

    window_step = window_width/overlap_factor

    dsep = '/'

    main_dir = working_dir

    grid = (main_dir + dsep + gridfilename)

    # prefix for output grids
    out_dir = main_dir + dsep + 'subsetted_grids'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if not os.path.exists(out_dir + dsep + 'window_' + str(window_width/1000)):
        os.makedirs(out_dir + dsep + 'window_' + str(window_width/1000))
        os.makedirs(main_dir + dsep + 'fft_output'
                    + dsep + 'window_' + str(window_width/1000))
        os.makedirs(main_dir + dsep + 'spectra_plots'
                    + dsep + 'window_' + str(window_width/1000))
        os.makedirs(main_dir + dsep + 'results'
                    + dsep + 'window_' + str(window_width/1000))
    outfile = out_dir + dsep + 'window_' + str(window_width/1000) + dsep + 'window'

    # load the input grid file
    # these 2 lines prevent cmd window opening suinfo = startupinfo
    suinfo = sp.STARTUPINFO()
    suinfo.dwFlags |= sp.STARTF_USESHOWWINDOW

    # run gdalinfo on input grid
    info_cmd = str(gdalinfo) + ' ' + grid
    sp.Popen(info_cmd, -1, stdout=sp.PIPE, stderr=sp.PIPE,
             startupinfo=suinfo, creationflags=0x08000000)
    # print info_process.communicate()[0]
    time.sleep(1)
    # sp.Popen.terminate(info_process)

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

    # os.system('taskkill /F /FI "WindowTitle eq C:\PROGRAM FILES\GDAL\gdalinfo.exe" /T')

    # setup windowing of the main grid

    # define initial window corners in the top left corner of the grid
    window_ulx = grid_topleftX
    window_uly = grid_topleftY
    window_lrx = window_ulx + window_width
    window_lry = window_uly - window_width

    # set subset counter
    subset = 1

    # process subsets loop moving from left to right and top to bottom along grid
    while window_uly > grid_lowerrightY:

        while window_ulx < grid_lowerrightX:
            print('processing subset #' + str(subset))
            cmd = (str(gdalTranslate) + ' -of GS7BG -projwin ' + str(window_ulx) +
                   ' ' + str(window_uly) + ' ' + str(window_lrx) + ' ' +
                   str(window_lry) + ' ' + grid + ' ' + outfile + '_' +
                   str(subset) + '.grd')
            sp.Popen(cmd, -1, stdout=sp.PIPE, stderr=sp.PIPE,
                     startupinfo=suinfo, creationflags=0x08000000)
            # print subset_process.communicate()[0]
            # update subset coords - move in +X direction
            window_ulx = window_ulx + window_step
            window_lrx = window_lrx + window_step
            time.sleep(1)
            # os.system('taskkill /F /FI "WindowTitle eq gdal_translate.exe" /T')
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
            info_cmd = str(gdalinfo) + ' ' + fn
            sp.Popen(info_cmd, -1, stdout=sp.PIPE, stderr=sp.PIPE,
                     startupinfo=suinfo, creationflags=0x08000000)
            # print info_process.communicate()[0]
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
            if len(test) == 0:
                print("empty")
            else:
                # plot using graphics.imshow library
                fig, ax = plt.subplots(figsize=(5, 5))
                graphics.imshow_hs(mag, ax, cmap_norm='equalize', hs=True,
                                   colorbar=True, cb_ticks='stats', nSigma=2,
                                   azdeg=45, altdeg=45, blend_mode='alpha',
                                   alpha=0.7, extent=(xmin, xmax, ymin, ymax),
                                   origin='upper', contours=False, labels=False,
                                   fmt='%.0f',
                                   levels=np.arange(round(np.nanmin(mag), 0),
                                                    round(np.nanmax(mag), 0),
                                                    2))
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

    os.system('taskkill /F /FI "WindowTitle eq C:\PROGRAM FILES\GDAL\gdalinfo.exe" /T')

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

    # plot using graphics.imshow library
    fig, ax = plt.subplots(figsize=(10, 10))
    ax = plt.subplot(aspect='equal')

    graphics.imshow_hs(mag, ax, cmap_norm='equalize', hs=True, colorbar=True,
                       cb_ticks='stats', nSigma=2, azdeg=45, altdeg=45,
                       blend_mode='alpha', alpha=0.7,
                       extent=(xmin, xmax, ymin, ymax),
                       origin='upper', contours=False, labels=False,
                       fmt='%.0f',
                       levels=np.arange(round(np.nanmin(mag), 0),
                                        round(np.nanmax(mag), 0), 2))

    ax.scatter(centerx, centery)

    cm_utils.rectangles(centerx, centery, window_width, window_width,
                        facecolor='None')

    plt.savefig(out_dir + dsep + 'window_' + str(window_width/1000) + dsep +
                'windows_centers.png', dpi=300)


def runfft(working_dir, window_dir):
    """
    Run the FFT on each grid.

    This function runs the fft on each of the grids created from splitgrid.
    It does not produce a file if the grid contains ANY nans.

    Requires GMT installation for grdfft
    On some machines grdfft is install here.
    grdfft = 'C:/programs/gmt5/bin/grdfft.exe'
    adjust as needed if you get a "cannot find the file specified" error

    Parameters
    ----------
    working_dir: main directory where the main grid file is located
    window_dir: name of the window_dir created when the grids were split

    Returns
    -------
    file: *.txt file with FFT output
    """
    grdfft = 'C:/Program Files/GMT/gmt5/bin/grdfft.exe'

    main_dir = working_dir
    griddir = main_dir + '/subsetted_grids/' + window_dir
    outdir = main_dir + '/fft_output/' + window_dir

    os.chdir(griddir)

    # run fft loop through all the grd files in the dir

    suinfo = sp.STARTUPINFO()
    suinfo.dwFlags |= sp.STARTF_USESHOWWINDOW

    for fn in os.listdir(griddir):
        if fn.endswith(".grd"):
            print("Doing fft on file:", fn)
            gmt_cmd = (grdfft + ' ' + fn + ' -G' + outdir + '/' +
                       os.path.splitext(fn)[0] + str('.txt ') + '-Nq -Er -V')
            sp.Popen(gmt_cmd, -1, stdout=sp.PIPE, stderr=sp.PIPE,
                     startupinfo=suinfo, creationflags=0x08000000)
            time.sleep(0.5)


def calccpd(working_dir, window_dir, Zo_start, Zo_end, Zt_start, Zt_end, Beta,
            K, CT):
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
            data = pd.read_csv(datadir + dsep + fn, sep='\t',
                               names=['Wavenumber', 'Power', 'Stdev'],
                               comment='#')

            # Calculate centroid depth, Zo
            # centroid fit points (k/2pi), may need to adjust values slightly if errors
            Zo_start = Zo_start
            Zo_end = Zo_end

            assert (Zo_start < Zo_end), 'Zo_start value must be smaller than Zo_end value'

            # calculate centroid
            data['y_centroid'] = np.log(data.Wavenumber**(Beta) * (data.Power / data.Wavenumber**2)) / (2*np.pi)  # this is the y axis value
            data['y_centroid_std'] = np.log(data.Wavenumber**(Beta) * (data.Stdev / data.Wavenumber**2)) / (2*np.pi)
            # data['y_centroid'] = np.log(data.Power**0.5 / np.abs(data.Wavenumber)) / (2*np.pi)  # this is the y axis value
            # data['y_centroid_std'] = np.log(data.Stdev**0.5 / np.abs(data.Wavenumber)) / (2*np.pi)
            data['wavenumber_2pi'] = np.abs(data.Wavenumber * 1000) / (2*np.pi)  # this is the x axis value (wavenumber/2pi) in cycles per km

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
            # data['y_top'] = np.log(data.Power**0.5) / (2*np.pi)
            # data['y_top_std'] = np.log(data.Stdev**0.5) / (2*np.pi)

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

            print(r'top depth = %0.1f \u00B1 %0.1f km' % (top_depth,
                                                          top_depth_err))
            print(r'centroid depth = %0.1f \u00B1 %0.1f km' % (centroid_depth,
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
            plt.tight_layout()

            # plot centroid
            plt.subplot(132)
            plt.scatter(data.wavenumber_2pi, data.y_centroid, marker='o', s=2,
                        label="Data")
            # plt.errorbar(data.wavenumber_2pi, data.y_centroid, yerr=data.y_centroid_std, marker='o', markersize=2, errorevery=10, label="Data")

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
            plt.tight_layout()

            # plot top
            plt.subplot(133)
            plt.scatter(data.wavenumber_2pi, data.y_top, marker='o', s=2,
                        label="Data")
            # plt.errorbar(data.wavenumber_2pi, data.y_top, yerr=data.y_top_std, marker='o', markersize=2, errorevery=10, label="Data")
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
            plt.tight_layout()
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
    plt.scatter(results.heatflow_table, results.base_depth_table)
    plt.xlabel('heatflow (mW/m$^2$)')
    plt.ylabel('Curie Pt depth (km)')
    plt.savefig(resultdir + 'heatflow_CPD.png', dpi=600)
    plt.show()
