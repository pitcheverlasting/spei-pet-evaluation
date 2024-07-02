#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script calculates the correlation between SM and SPEI, used after R SPEI package is run.
SPEI package could be installed under R 3.X, not R 4.X
Next: plot_climatology.py/plot_correlations.py/plot_time_series.py
"""
__author__ = "Liqing Peng"
__copyright__ = "Copyright (C) 2024 Liqing Peng"
__license__ = "MIT"
__version__ = "2024.05"


import os, DataHandeling
from pylab import *
from netCDF4 import Dataset
import pandas as pd
from scipy.stats import spearmanr, pearsonr

wet_home = '/server/user/Data/PROJECT/WET_US'
dims_025 = DataHandeling.Set_dims(0.25, 112, 232, 25.125, -124.875)
# exps = ['run_2021']
exps = ['run_2024']
SM = 'SMsurf'

petfiles = []
model_modes = ['BL', 'TS']
albedo_modes = ['CLM']
ga_modes = ['LC', 'CH']
gs_modes = ['K1995', 'Zhou2006']

for albedo_mode in albedo_modes:
	for model_mode in model_modes:
		if model_mode == 'BL':
			for gs_mode in gs_modes:
				for ga_mode in ga_modes:
					filename = '%s.ga_%s.gs_%s.alb_%s' %(model_mode, ga_mode, gs_mode, albedo_mode)
					petfiles.append(filename)
		else:
			filename = '%s.SW.alb_%s' %(model_mode, albedo_mode)
			petfiles.append(filename)


def get_start_end_date(styr, edyr):
	"To select growing season"
	stdy = datetime.datetime(styr, 1, 1); eddy = datetime.datetime(edyr, 12, 31)
	st_idx = (stdy.year - 1981) * 12 + stdy.month - 1
	ed_idx = (eddy.year - 1981) * 12 + eddy.month - 1
	#relativedelta(datetime.datetime(1984, 1, 1)).months  # find out the start date from the gridded datasets
	dates = pd.date_range((str(styr)+'-01-01'), (str(edyr)+'-12-31'), freq='M')
	ddays = pd.to_datetime(dates); dyears = ddays.year; dmonths = ddays.month
	return [stdy, eddy, st_idx, ed_idx, ddays, dyears, dmonths]

stdy, eddy, st_idx, ed_idx, ddays, dyears, dmonths = get_start_end_date(1981, 2017)


def calculate_monthly_growsingseason_correlation_spei_SM(expdir, file, step):
	"Get growing season monthly data"
	mondir = '%s/%s/monthly' % (wet_home, expdir)
	corrdir = '%s/%s/evaluation/spei/growingseason' % (wet_home, expdir)
	if not os.path.exists(corrdir):
		os.makedirs(corrdir)

	"SPEI different time scales"
	outfile = '%s/spei.%s.mon%02d.nc' % (mondir, file, step)
	corrfile = '%s/corr.spei.%s.%s.mon%02d.1981-2017.nc' % (corrdir, SM, file, step)

	index = (dmonths>=4) & (dmonths<=10)

	data = Dataset(outfile).variables['spei'][:]
	data = data[index, :, :]

	smfile = '%s/run/monthly/%s.monthly.MAmon%02d.1981-2017.nc' % (wet_home, SM, step)
	surf = Dataset(smfile).variables[SM][:]
	surf = surf[index ,:, :]

	mask = (~isnan(surf)).sum(0)
	corr = np.empty((2, 112, 232)); corr.fill(np.nan)

	y_valid, x_valid = np.where(mask>10)

	for y, x in zip(y_valid, x_valid):
		A = data[:, 111-y, x] # get the flip of the lat grid
		B = surf[:, y, x]

		index = (~isnan(A) & ~isnan(B) & (A!=0) & (B!=0))

		if index.sum() > 30:
			# corr[:2, y, x] = spearmanr(A[index], B[index])
			corr[:2, y, x] = pearsonr(A[index], B[index])
		else:
			corr[:2, y, x] = np.nan
		print y, x

	DataHandeling.Create_Write_NETCDF_File(dims_025, corrfile, 'corr', 'corr', corr, datetime.datetime(2000, 1, 1), 2, 'days')

	return

for exp in exps[:1]:
	for file in petfiles:
		for step in [1,3,6,12]:
			calculate_monthly_growsingseason_correlation_spei_SM(exp, file, step)


