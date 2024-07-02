#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script prepares resampled daily PET to monthly for SPEI analysis
Next: R SPEI (see run_SPEI_R.txt)
Then: SM_SPEI_correlation.py
"""
__author__ = "Liqing Peng"
__copyright__ = "Copyright (C) 2024 Liqing Peng"
__license__ = "MIT"
__version__ = "2024.05"


import os
from pylab import *

cdo = '/server/user/miniconda/bin/cdo'
wet_home = '/server/user/Data/PROJECT/WET_US'
sample = '%s/run/monthly/SMsurf.monthly.1981-2017.nc' % (wet_home)

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
# print(petfiles)
# exit()

exps = ['run_2024']

def daily2multiscale_for_spei(expdir, file):

	mondir = '%s/%s/monthly' % (wet_home, expdir)
	anndir = '%s/%s/yearly' % (wet_home, expdir)
	tempdir = '%s/%s/monthly_0125' % (wet_home, expdir)

	if not os.path.exists(mondir):
		os.makedirs(mondir)

	if not os.path.exists(tempdir):
		os.makedirs(tempdir)

	if not os.path.exists(anndir):
		os.makedirs(anndir)


	# Daily to Monthly
	for year in xrange(1981, 2018):
		petfile = '%s/%s/pet.%s.%04d.nc' % (wet_home, expdir, file, year)
		monfile = '%s/pet.%s.monthly.%04d.nc' % (tempdir, file, year)
		os.system('%s -P 20 monmean %s %s' %(cdo, petfile, monfile))

	# Merge monthly to one file
	petfile = '%s/pet.%s.monthly.*.nc' % (tempdir, file)
	monfile = '%s/pet.%s.monthly.1981-2017.nc' % (tempdir, file)
	os.system('%s -P 10 mergetime %s %s' % (cdo, petfile, monfile))

	# Remap 0.125 to 0.25
	petoutfile = '%s/pet.%s.025.monthly.1981-2017.nc' % (mondir, file)
	os.system('%s -P 20 remapcon,%s %s %s' % (cdo, sample, monfile, petoutfile))
	# os.system('rm %s' % (monfile))

	# Yearly file
	annfile = '%s/pet.%s.yearly.1981-2017.nc' % (anndir, file)
	os.system('%s -P 10 yearmonmean %s %s' % (cdo, petoutfile, annfile))


	return


for exp in exps[:1]:
	for file in petfiles:
		daily2multiscale_for_spei(exp, file)

# # merge the datasets for comparison
# # Merge monthly to one file
# extdir = '/server/user/Data/PET_SW_Sun'
# petfile = '%s/SW_PET_PT_PE_*_remap.nc' % (extdir)
# monfile = '%s/pet.Sun.monthly.1982-2015.nc' % (extdir)
# # https://code.mpimet.mpg.de/boards/2/topics/12741?r=12746, mergetime and select only one variable
# # os.system('%s -P 10 mergetime -apply,-selname,PET [ %s ] %s' % (cdo, petfile, monfile)) # this may be for newer CDO
# os.system('%s -P 10 mergetime %s %s' % (cdo, petfile, monfile))
