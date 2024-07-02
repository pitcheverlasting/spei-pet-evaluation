#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script provides plotting functions for comparison with Sun et al. (2023)
"""

from pylab import *
from netCDF4 import Dataset
import pandas as pd
import matplotlib.cbook as cbook

import DataHandeling
import Plotting

wet_home = '/server/user/Data/PROJECT/WET_US'
dir19 = '%s/run_2019' % (wet_home)
dir21 = '%s/run_2021' % (wet_home)
dir24 = '%s/run_2024' % (wet_home)
extdir = '/server/user/Data/PET_SW_Sun'
workspace = '%s/workspace_2024' %(wet_home)
figdir = '%s/Figure/202405' %wet_home
para_dir = '/server/user/Data/forcing/NLDAS2/parameter/CONUS'

# dims_025 = DataHandeling.Set_dims(0.25, 112, 232, 25.125, -124.875)
mk025file = '/server/user/Masks/0.25deg/us_0.25deg.nc'
mk25 = Dataset(mk025file).variables['data'][0,: :]
# set the mask missing values to nan. otherwise it will set to zero when hstacking
mk25[mk25 == 0] = np.nan

methods = ['ow.025', 'pt.025', 'faoshort.025', 'faotall.025', 'BL.ga_LC.gs_K1995.alb_CLM.025', 'BL.ga_CH.gs_K1995.alb_CLM.025', 'TS.SW.alb_CLM.025', 'Sun']
years = [(1981, 2017)] * 7 + [(1982, 2015)]  # (1980, 2022)
stmons = [12] * 7 + [0] # 24
edmons = [420] * 7 + [408] # 432
petdirs = [dir19+'/monthly', dir19+'/monthly', dir19+'/monthly', dir19+'/monthly', dir21+'/monthly', dir24+'/monthly', dir24+'/monthly', extdir]
labels = ['PET-OW', 'PET-PT', 'PET-RC-short', 'PET-RC-tall', 'PET-LC-K', 'PET-CH-K', 'PET-SW',  'PET-Sun'] #, 'PET-ERA5']
cols = ['SkyBlue', 'Orange', 'DarkOrange', 'OrangeRed', 'lightseagreen', 'Green', 'Sienna', 'RoyalBlue']
flags = [False] * 7 + [True]
variables = ['pe'] * 7 + ['PET'] # 'pet']

cmap = plt.get_cmap('nipy_spectral')
avg_cmap = Plotting.truncate_colormap(cmap, 0.2, 0.9)  # 0.15, 0.85)

def read_ETdata_usmask(ncdata):
	if len(ncdata.shape) > 2:
		print("only input 2D array for read_ncdata_mask_for_mapping()!")
		exit()
	return ncdata[12:,:][::-1] * mk25

def read_corrdata_usmask(ncdata):
	if len(ncdata.shape) > 2:
		print("only input 2D array for read_ncdata_mask_for_mapping()!")
		exit()
	return ncdata[:100,:] * mk25


def Map_climatology_8panels():

	captions = [chr(i) for i in range(ord('a'),ord('z')+1)]

	fig = plt.figure(figsize=(15, 8))
	for i in xrange(len(methods)):
		file = '%s/pet.%s.monthly.%s-%s.nc' % (petdirs[i], methods[i], years[i][0], years[i][1])
		data = Dataset(file).variables[variables[i]][stmons[i]: edmons[i], :]
		gs = DataHandeling.mon2growseason(data, 1982, 2015)
		dataplot = read_ETdata_usmask(nanmean(gs, axis=0))
		if i < len(methods)-1: dataplot = dataplot*30.4 # from month to daily
		ax = fig.add_subplot(3, 3, i+1)
		if i%3 == 0: latlabel = True
		else: latlabel = False
		lonlabel = True
		ax.m = Plotting.setup_USmap(latlabel, lonlabel)

		im = ax.m.imshow(dataplot, vmin=0, vmax=250, cmap=avg_cmap, alpha=0.7)
		# im = ax.m.imshow(season_avg[i][::-1], vmin=0, vmax=250, cmap=avg_cmap, alpha=0.7)
		# if i == 0:
		# 	ax.m.colorbar(im, location='bottom', pad='15%')
		ax.set_title('(%s) %s' % (captions[i], labels[i]), fontsize=16)

	cax = fig.add_axes([0.69, 0.2, 0.28, 0.04])
	cb = fig.colorbar(im, cax, orientation='horizontal')  # adjust the size
	# cb.ax.tick_params(labelsize=12)   # change the colorbar fontsize
	loc = np.arange(0, 250.1, 25)
	llabels = ['%d' %ll for ll in loc]
	llabels[-1] = '>250'
	cb.set_ticks(loc)
	cb.ax.set_xticklabels(llabels, fontsize=14)

	fig.tight_layout()
	fig.subplots_adjust(left=0.04, bottom=0.02, wspace=0.0)
	# plt.savefig('%s/FigC1_PET_all_Sun_spatial_8panels.pdf' %(figdir), dpi=300)
	# plt.savefig('%s/FigC1_PET_all_Sun_spatial_8panels.png' %(figdir), dpi=300)
	plt.show()
	return

# Map_climatology_8panels()
# exit()


def Time_series_PET_compare():
	"Draw time series to compare with PET Sun and ERA, the common period is 1982-2015, 34 years"

	filenames = []
	for method, dir, year in zip(methods, petdirs, years):
		filenames.append('%s/pet.%s.monthly.%s-%s.nc' % (dir, method, year[0], year[1]))

	def get_ann_ts(file, var, stmon, edmon, flag=True):
		# convert mm/d to mm/mon
		us = Dataset(file).variables[var][stmon:edmon,12:,:] * mk25[::-1]

		data = us.reshape(34, 12, -1)
		if flag == True:
			ann = sum(data, axis=1)
		else:
			ann = sum(data*30, axis=1)
		ann[ann.mask==True] = np.nan
		median = nanmean(ann, axis=1)

		return median

	fig, ax = plt.subplots(figsize=(12,4.5))

	# For PET
	for i in range(0,len(methods)):
		ts = get_ann_ts(filenames[i], variables[i], stmons[i], edmons[i], flag=flags[i])
		ax.plot(ts, color=cols[i], linewidth=1.5, label='(%s %d mm yr$^{-1}$)' % (labels[i], nanmean(ts)), alpha=0.9)

	ax.set_xticks(arange(-2, 35, 5))
	ax.set_xticklabels(arange(1980, 2016, 5))
	ax.set_ylim([500, 3000])
	plt.minorticks_on()
	ax.tick_params(axis='both', labelsize=16)
	ax.set_ylabel('Annual PET (mm yr$^{-1}$)', fontsize=18)

	plt.legend(loc=0, ncol=2, fontsize=14, frameon=False)

	fig.tight_layout()
	fig.subplots_adjust(left=0.1)
	plt.savefig('%s/FigC2_PET_LC_Sun_ERA_ann_time_series_US.pdf' %(figdir), dpi=300)

	# plt.show()

	return

# Time_series_PET_compare()