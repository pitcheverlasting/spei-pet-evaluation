#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script provides plotting functions for time series of PET/P/SPEI
"""
__author__ = "Liqing Peng"
__copyright__ = "Copyright (C) 2024 Liqing Peng"
__license__ = "MIT"
__version__ = "2024.05"

from netCDF4 import Dataset
from pylab import *
from matplotlib import colors
import os, DataHandeling
import scipy
from scipy.stats import spearmanr, pearsonr
import pandas as pd

rc('font', family='serif')

wet_home = '/server/user/Data/PROJECT/WET_US'
workspace = '%s/workspace_2024' %(wet_home) # 2021
# figdir = '%s/Figure/202012' %wet_home
figdir = '%s/Figure/202405' %wet_home
mk025file = '/server/user/Masks/0.25deg/us_0.25deg.nc'
mk25 = Dataset(mk025file).variables['data'][0,: :]
# set the mask missing values to nan. otherwise it will set to zero when hstacking
mk25[mk25 == 0] = np.nan

petfiles = ['control']
dir21 = '%s/run_2021' % (wet_home)
dir20 = '%s/run_2020' % (wet_home)
dir19 = '%s/run_2019' % (wet_home)
dir24 = '%s/run_2024' % (wet_home)
dirold = '%s/run' % (wet_home)
extdir = '/server/user/Data/PET_SW_Sun'
SM = 'SMsurf'
steps = [1, 3, 6, 12]

methods = ['ow', 'pt', 'faoshort', 'faotall', 'BL.ga_LC.gs_K1995.alb_CLM', 'BL.ga_CH.gs_K1995.alb_CLM', 'TS.SW.alb_CLM']
petdirs = [dir19, dir19, dir19, dir19, dir21,  dir24, dir24]
labels = ['PET-OW', 'PET-PT', 'PET-RC-short', 'PET-RC-tall', 'PET-LC', 'PET-CH', 'PET-SW']
labels2 = ['OW', 'PT', 'RC-short', 'RC-tall', 'LC', 'CH', 'SW']
# linecolors = ['Brown', 'Orange', 'DarkOrange', 'OrangeRed', 'Crimson', 'Red', 'Maroon']
linecolors = ['Dodgerblue', 'Orange', 'DarkOrange', 'OrangeRed', 'lightseagreen', 'green', 'Sienna']
adjustments = [0, -150, -150, 60, -120, 50, 50]

speidirs = petdirs


def prepare_Time_series_ETm_Pre():
	"Draw time series"
	variables = ['pr', 'pe']
	datanames = ['P', 'PET']

	filenames = ['%s/monthly/pre.025.monthly.1981-2017.nc'  % (dirold)]
	for method, dir in zip(methods, petdirs):
		filenames.append('%s/monthly/pet.%s.025.monthly.1981-2017.nc' % (dir, method))

	def get_ann_ts(file, var, flag=True):
		# convert mm/d to mm/mon
		us = Dataset(file).variables[var][:,12:,:] * mk25[::-1]
		data = us.reshape(37, 12, -1)
		if flag == True:
			ann = mean(data, axis=1)
		else:
			ann = sum(data*30, axis=1)
		ann[ann.mask==True] = np.nan
		# median = nanpercentile(ann, 50, axis=1)
		median = nanmean(ann, axis=1)

		return median

	df = {}
	i = 0
	df['P (mm d-1)'] = get_ann_ts(filenames[i], variables[i], flag=True)

	# For PET
	for i in range(1,len(methods)+1):
		df['%s (mm d-1)'%labels[i-1]] = get_ann_ts(filenames[i], variables[1], flag=False)

	df_final = pd.DataFrame(df, index=pd.date_range('1981', '2017', freq='AS'))
	df_final.to_csv('%s/annual_P_PET_US_1981-2017.csv' %workspace)

	return

# prepare_Time_series_ETm_Pre()
# exit()

def Time_series_ETm_Pre():
	"Draw time series"
	variables = ['pr', 'pe']
	cols = ['RoyalBlue', 'Crimson']

	filenames = ['%s/monthly/pre.025.monthly.1981-2017.nc'  % (dirold)]
	for method, dir in zip(methods, petdirs):
		filenames.append('%s/monthly/pet.%s.025.monthly.1981-2017.nc' % (dir, method))

	def get_ann_ts(file, var, flag=True):

		# convert mm/d to mm/mon
		us = Dataset(file).variables[var][:,12:,:] * mk25[::-1]

		data = us.reshape(37, 12, -1)
		if flag == True:
			ann = mean(data, axis=1)
		else:
			ann = sum(data*30, axis=1)
		ann[ann.mask==True] = np.nan
		# median = nanpercentile(ann, 50, axis=1)
		median = nanmean(ann, axis=1)

		return median

	fig, ax = plt.subplots(figsize=(12,4.5))

	### for precipitation
	i = 0
	pre = get_ann_ts(filenames[i], variables[i], flag=False)
	map = nanmean(pre)
	ax.plot(pre, color=cols[i], linewidth=3, label='$P$ (MA = %d mm yr$^{-1}$)' %map)
	# trend
	z = np.polyfit(range(37), pre, 1)
	p = np.poly1d(z)
	result = scipy.stats.linregress(range(37), pre)
	slp, R, P = result[0], result[2], result[3]
	ax.plot(range(37), p(range(37)), '--', c=cols[i], label='$P$ (trend = %.02f mm yr$^{-1}$ yr$^{-1}$)' %slp, linewidth=2, alpha=0.5)

	ax.set_xticks(arange(-1, 38, 5))
	ax.set_xticklabels(arange(1980, 2018, 5))
	ax.set_ylim([0, 2600])
	plt.minorticks_on()
	ax.tick_params(axis='both', labelsize=16)
	ax.set_ylabel('Annual P or PET (mm yr$^{-1}$)', fontsize=18)

	# For PET
	AED_ensemble = []
	for i in range(1,len(methods)+1):
		AED_ensemble.append(get_ann_ts(filenames[i], variables[1], flag=False))
		ax.plot(get_ann_ts(filenames[i], variables[1], flag=False), color=linecolors[i-1], linewidth=1.5, label='_skip_legend_', alpha=0.9)
		ax.text(19, get_ann_ts(filenames[i], variables[1], flag=False)[19]+adjustments[i-1], '%s' %labels[i-1], color=linecolors[i-1], fontsize=14)


	AED_ensemble = array(AED_ensemble)
	# median = nanpercentile(AED_ensemble, 50, axis=0)
	median = nanmean(AED_ensemble, axis=0)
	mape = nanmean(median)

	top = nanmax(AED_ensemble, axis=0)
	bottom = nanmin(AED_ensemble, axis=0)
	ax.plot(median, color=cols[1], linewidth=3, label=r'$\overline{PET}$ (MA = %d mm yr$^{-1}$)' %mape, alpha=0.8)
	ax.fill_between(range(0, 37), top, bottom, color=cols[1], alpha=0.1)

	# trend
	z = np.polyfit(range(37), median, 1)
	p = np.poly1d(z)
	result = scipy.stats.linregress(range(37), median)
	slp, R, P = result[0], result[2], result[3]
	ax.plot(range(37), p(range(37)), '--', c=cols[1], label=r'$\overline{PET}$ (trend = %.02f mm yr$^{-1}$ yr$^{-1}$)' %slp, linewidth=2, alpha=0.5)

	ymin, ymax = ax.get_ylim()
	ax.text(1, ymin + 0.9 * (ymax - ymin), '(a)', fontsize=18, ha='left')

	plt.legend(loc=1, ncol=2, frameon=False)
	fig.tight_layout()
	fig.subplots_adjust(left=0.1)
	# plt.savefig('%s/Fig5a_PET_Pre_ann_time_series_US.pdf' %(figdir), dpi=300)
	# plt.savefig('%s/Fig5a_PET_Pre_ann_time_series_US.png' %(figdir), dpi=300)

	plt.show()

	return

# Time_series_ETm_Pre()
# exit()


def prepare_Time_series_SPEI(step):
	"Draw time series"
	variables = [SM, 'spei']

	filenames = []
	filenames.append('%s/run/monthly/%s.monthly.MAmon%02d.1981-2017.nc' % (wet_home, SM, step))

	for method, dir in zip(methods, speidirs):
		filenames.append('%s/monthly/spei.%s.mon%02d.nc' % (dir, method, step))

	def get_monthly_ts(file, var, mask_flag=True):
		if mask_flag == True:
			us = Dataset(file).variables[var][:, 12:, :] * mk25[::-1]
		else:
			us = Dataset(file).variables[var][:, 12:, :] * mk25
		data = us.reshape(37 * 12, -1)
		data[data.mask == True] = np.nan
		mean = nanmean(data, axis=1)
		return mean

	df = {}

	df['SM %s-month MA'%step] = get_monthly_ts(filenames[0], variables[0], mask_flag=False)
	for i in range(1, len(filenames)):
		df['SPEI-%s %s'%(step, labels[i-1])] = get_monthly_ts(filenames[i], variables[1])

	df_final = pd.DataFrame(df, index=pd.date_range('1981-01-01', '2017-12-31', freq='MS'))
	# print(df_final)
	df_final.to_csv('%s/monthly_SM_SPEI-%s_US_1981-2017.csv' %(workspace, step))

	return


# [prepare_Time_series_SPEI(step) for step in steps]
# exit()


def Time_series_SPEI():

	fig, ax = plt.subplots(2,2, figsize=(20, 7))

	captions = [chr(i) for i in range(ord('a'), ord('z') + 1)]
	def plot_step(ax, step, caption):
		input = pd.read_csv('%s/monthly_SM_SPEI-%s_US_1981-2017.csv' % (workspace, step))

		# For ax2
		ax2 = ax.twinx()
		ax2.plot(input['SM %s-month MA' % step], color='Grey', label='SMsurf', linewidth=2, alpha=0.8)
		if step == 12:
			ax2.set_ylim([0.278, 0.328])
		elif step == 6:
			ax2.set_ylim([0.22, 0.38])

		ax2.set_ylabel('SMsurf (m$^3$ m$^{-3}$)', color='k', fontsize=16, labelpad=30, rotation=270)
		ax2.tick_params(axis='y', colors='k')

		for i in range(len(speidirs)):
			ax.plot(input['SPEI-%s %s'%(step, labels[i])], label='SPEI %s'%(labels2[i]), color=linecolors[i], linewidth=1.5)

		# For ax
		ax.set_ylim([-1.5, 1.5])
		ax.tick_params(axis='y') #, colors='DarkRed')
		ax.set_ylabel('SPEI-%s' %step, fontsize=18) #, color='DarkRed')
		ax.axhline(y=0, color='grey', linestyle='-', linewidth=1)

		ymin, ymax = ax.get_ylim()
		ax.text(12, ymin + 0.9 * (ymax - ymin), '(%s)'%caption, fontsize=16, ha='left')
		plt.xticks(arange(-1, 444, 12*5), arange(1980, 2017, 5))

		for axis in [ax, ax2]:
			axis.tick_params(axis='both', labelsize=16)
			# plt.minorticks_on() # This messed up the axis

		if step == 12:
			lines, labs = ax.get_legend_handles_labels()
			lines2, labs2 = ax2.get_legend_handles_labels()
			ax2.legend(lines + lines2, labs + labs2, loc=0, ncol=4, frameon=False, fontsize=11)

		return ax

	plot_step(ax[0, 0], 1, captions[1])
	plot_step(ax[0, 1], 3, captions[2])
	plot_step(ax[1, 0], 6, captions[3])
	plot_step(ax[1, 1], 12, captions[4])

	fig.tight_layout()
	fig.subplots_adjust(left=0.05, right=0.93)
	plt.savefig('%s/Fig5_SPEI_SM_monthly_time_series_US.pdf' %(figdir), dpi=300)
	plt.savefig('%s/Fig5_SPEI_SM_monthly_time_series_US.png' %(figdir), dpi=300)

	# plt.show()

	return

# Time_series_SPEI()
# exit()