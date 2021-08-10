# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 11:56:47 2021

@author: craigm
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

def calccpd_geosoft(working_dir, window_dir, Zo_start, Zo_end, Zt_start, Zt_end, Beta,
            K=2.5, CT=570):
    """
    Calculate the Curie Point Depth.

    Read in the output from the geosoft raidal powerspectrum.
    Use the output in wavenumbers not wavelengths.

    Loops through each file in directory and then incorporates the centre
    coords of each window to the results.

    Incorporates the fractal correction "Beta" from Bansal et al 2011,
    Geophysics vol 76
    
    Geosoft power spectrum units , cycles/1000*m i.e cycles per km
    Power spectrum has been normalised by subtract average spectrum density

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
                               names=['Wavenumber', 'Num_samples', 'Ln_P',
                                      'Depth_3', 'Depth_5'],
                               comment='/')

            ave_power = pd.read_csv(datadir + dsep + fn,
                                    nrows=1, header=None, skiprows=[0, 1],
                                    delim_whitespace=True)[6]

            data['Ln_raw'] = data.Ln_p - ave_power.values
            
            data['Power'] = np.exp(data.Ln_P)

            print(data.columns)
            # Calculate centroid depth, Zo
            # centroid fit points (k*1000/2pi), may need to adjust values
            # slightly if index errors, they need to match actual values
            Zo_start = Zo_start
            Zo_end = Zo_end

            assert (Zo_start < Zo_end), 'Zo_start value must be smaller than Zo_end value'

            # calculate centroid
            # this is the y axis value Eqn7 Bansal
            data['y_centroid'] = np.log(data.Wavenumber**(Beta) *
                                        (data.Power / data.Wavenumber**2)) / (2*np.pi)

            # this is the x axis value (wavenumber/2pi) in cycles per km
            data['wavenumber_2pi'] = np.abs(data.Wavenumber) / (2*np.pi)

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

            # calculate top depth Bansal Eqn 8
            data['y_top'] = np.log(data.Wavenumber**(Beta) * data.Power) / (2*np.pi)

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
            plt.scatter(data.Wavenumber, data.Power, marker='o', s=2,
                        label="Data")
            plt.xlim(0, data.Wavenumber.max())
            plt.xlabel(r'$|k| [km^{-1}$]')
            plt.ylabel('Power')
            plt.title(fn.split('.')[0])

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