#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script runs the PET models
Next: prepare_monthly_ETmax.py
"""
__author__ = "Liqing Peng"
__copyright__ = "Copyright (C) 2024 Liqing Peng"
__license__ = "MIT"
__version__ = "2024.05"

from pylab import *
import pandas as pd
from netCDF4 import Dataset
import time, DataHandeling
from ET_library import PET


##########################################  IO: path, variables, dimensions

# forcing directories
dailyforcing_dir = '/server/user/Data/forcing/NLDAS2/daily_forcing/CONUS'
wet_home = '/server/user/Data/PROJECT/WET_US'
run_dir = 'run_2024'
dims_nldas = DataHandeling.Set_dims(0.125, 224, 464, 25.0625, -124.9375)

# meteorological forcing:
rsds = {'dirname': 'rsds', 'filename': '_daily_', 'varname': 'rsds', 'longname':'short-wave_downwards_surface_radiation_flux', 'unit':'W/m2'}
rlds = {'dirname': 'rlds', 'filename': '_daily_', 'varname': 'rlds', 'longname':'long-wave_downwards_surface_radiation_flux', 'unit':'W/m2'}
pres = {'dirname': 'ps', 'filename': '_daily_', 'varname': 'ps', 'longname':'surface_pressure', 'unit':'Pa'}
tavg = {'dirname': 'tavg_scale', 'filename': '_daily_', 'varname': 'tavg', 'longname':'average_air_temperature', 'unit':'K'}
tmax = {'dirname': 'tmax_scale', 'filename': '_daily_', 'varname': 'tmax', 'longname':'maximum_air_temperature', 'unit':'K'}
tmin = {'dirname': 'tmin_scale', 'filename': '_daily_', 'varname': 'tmin', 'longname':'minimum_air_temperature', 'unit':'K'}
rhum = {'dirname': 'hurs_scale', 'filename': '_daily_', 'varname': 'hurs', 'longname':'relative_humidity', 'unit':'%'}
wind = {'dirname': 'wind', 'filename': '_daily_', 'varname': 'wind', 'longname':'wind', 'unit':'m/s'}
prec = {'dirname': 'pr', 'filename': '_daily_', 'varname': 'pr', 'longname':'precipitation', 'unit': 'kg m-2 s-1'}

# Albedo
albd = {'dirname': 'albedo', 'filename': '.glass.daily.', 'varname': 'albedo', 'longname':'GLASS albedo 8day to daily', 'unit':''}

# vegetation forcing:
para_dir = '/server/user/Data/forcing/NLDAS2/parameter/CONUS'
#1. IGBP: /server/user/Data/forcing/NLDAS2/parameter/CONUS/land_cover_climatology_nldas.nc
lc = {'filename': 'nldas_land_cover_climatology.nc', 'varname': 'landcover', 'longname':'"MODIS global land cover climatology (0.5km to 0.125deg)', 'unit':''}
chmin = {'filename': 'nldas_canopy_height_Lang_minimum.nc', 'varname': 'chmin', 'longname':'Minimum canopy height combining Lang et al. (2023) and look up table', 'unit':'m'}
chmax = {'filename': 'nldas_canopy_height_Lang_maximum.nc', 'varname': 'chmax', 'longname':'Maximum canopy height combining Lang et al. (2023) and look up table', 'unit':'m'}
z0m = {'filename': 'nldas_roughness_length_momentum_CH_Lang.nc', 'varname': 'Z0m', 'longname':'Roughness length for momentum', 'unit':'m'}
d0 = {'filename': 'nldas_zeroplane_displacement_height_CH_Lang.nc', 'varname': 'd0', 'longname':'Zeroplane displacement height', 'unit':'m'}
z0h = {'filename': 'nldas_roughness_length_heat_CH_Lang.nc', 'varname': 'Z0h', 'longname': 'Roughness length for heat', 'unit':'m'}
z0g = {'filename': 'nldas_roughness_length_ground.nc', 'varname': 'z0g', 'longname': 'Roughness length for ground', 'unit':'m'}
albclm = {'filename': 'nldas_albedo_climatology_366day_1982-2012.nc', 'varname': 'albedo', 'longname': '366 daily climatology albedo', 'unit':''}
laiclm = {'filename': 'nldas_LAI_climatology_366day_1982-2016.nc', 'varname': 'LAI', 'longname': '366 daily climatology LAI', 'unit':'m2 m-2'}
gstmax = {'filename': 'nldas_maximum_stomatal_conductance.nc', 'varname': 'gst_max', 'longname': 'Maximum stomatal conductance', 'unit':'mm/s'}
laimax = {'filename': 'nldas_maximum_LAI.nc', 'varname': 'LAImax', 'longname': 'Maximum LAI', 'unit':'m2 m-2'}
rstmin = {'filename': 'nldas_stomatal_resistance_minimum.nc', 'varname': 'rstmin', 'longname': 'Minimum stomatal resistance from land cover look up table', 'unit':'s/m'}
kb = {'filename': 'nldas_kB_land_cover.nc', 'varname': 'kB', 'longname': 'kB -1 term', 'unit':''}


def Read_forcing_daily(year, day):
	vars = ['SWin', 'LWin', 'Albedo', 'Tavg', 'Tmax', 'Tmin', 'RHum', 'Pres', 'Wind']
	dicts = [rsds, rlds, albd, tavg, tmax, tmin, rhum, pres, wind]
	datasets = {}

	for var, dict in zip(vars, dicts):
		file = '%s/%s/%s%s%04d.nc' %  (dailyforcing_dir, dict['dirname'], dict['varname'], dict['filename'], year)
		datasets[var] = Dataset(file).variables[dict['varname']][day, :, :]

	return datasets

def Read_Prec_daily(year, day):
	dict = prec
	file = '%s/%s/%s%s%04d.nc' %  (dailyforcing_dir, dict['dirname'], dict['varname'], dict['filename'], year)
	data = Dataset(file).variables[dict['varname']][day, :, :] * 60 * 60 * 24

	return data

def Read_parameter():
	# dictname in the ET_library
	vars = ['Landcover', 'z0m', 'd0p', 'z0h', 'z0g', 'gstmax', 'rstmin', 'chmin', 'chmax', 'LAImax', 'kB']
	dicts = [lc, z0m, d0, z0h, z0g, gstmax, rstmin, chmin, chmax, laimax, kb]
	datasets = {}

	for var, dict in zip(vars, dicts):
		file = '%s/%s' %  (para_dir, dict['filename'])
		datasets[var] = Dataset(file).variables[dict['varname']][0, :, :]

	return datasets


def Read_parameter_climatology(day, input):
	vars2 = ['Albedo_clm', 'LAI_clm']  # climatology
	dicts2 = [albclm, laiclm]

	for var, dict in zip(vars2, dicts2):
		file = '%s/%s' %  (para_dir, dict['filename'])
		input[var] = Dataset(file).variables[dict['varname']][day, :, :]

	return input


def Calculate_PE_daily_method_run(year):
	"""Calculate PET for each daily time step
	Input:
		year
	Output:
		monthly PE writing to netcdf
	"""

	veg_param = Read_parameter()

	# Read in daily data
	date = datetime.datetime(year, 1, 1, 0, 0)
	day = (date - datetime.datetime(year, 1, 1)).days

	while date <= datetime.datetime(year, 12, 31, 0, 0):

		T0 = time.time()
		input = Read_forcing_daily(year, day)
		veg_param = Read_parameter_climatology(day, veg_param)

		T1 = time.time()
		print("Read in forcing: %04d-%02d ......%s seconds " % (year, day, T1-T0))

		# There are 4 dimensions
		# 1. Partitioning: big leaf vs two source
		# 2. GA: open water, FAO GA, land cover GA
		# 3. GC/GS: open water, FAO GC, land cover LAI varying GC (only LAI climatology)
		# 4. Albedo: open water albedo, FAO albedo, varying albedo

		# Set up albedo scenarios
		cls_alb_clim = PET(input, veg_param, albedo_mode='CLM', lai_mode='CLM')


		#### 5 Experiments
		ga_modes = ['LC', 'CH']
		gs_modes = ['K1995', 'Zhou2006']
		for gs_mode in gs_modes:
			for ga_mode in ga_modes:
				# Big leaf
				filename = wet_home + '/%s' %(run_dir) + '/pet.BL.ga_%s.gs_%s.alb_CLM.%04d' %(ga_mode, gs_mode, year) + '.nc'
				print(year, filename)
				DataHandeling.Write_NETCDF_File(dims_nldas, filename, 'pe', 'PET Big Leaf', cls_alb_clim.Penman_Monteith_Big_Leaf(ga_mode=ga_mode, gs_mode=gs_mode), date, "hours", "mm day-1")

		# Two source
		filename = wet_home + '/%s' %(run_dir) + '/pet.TS.SW.alb_CLM.%04d' %(year) + '.nc'
		DataHandeling.Write_NETCDF_File(dims_nldas, filename, 'pe', 'PET Shuttleworth Wallace', cls_alb_clim.Shuttleworth_Wallace_Canopy_Soil(), date, "hours", "mm day-1")
		# print(year, filename)

		print("Daily PET saved!......%s seconds " % (time.time()-T1))

		date = date + relativedelta(days=1)
		day = day + 1
		print(date)

	return

if len(sys.argv)==1: styr=1981; edyr=2017
else: styr = int(sys.argv[1]); edyr = int(sys.argv[2])
[Calculate_PE_daily_method_run(year) for year in xrange(styr, edyr+1)]
