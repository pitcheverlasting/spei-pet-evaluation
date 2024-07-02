#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This library provides analysis tools for:
1. Lat-lon
2. NetCDF
3. Time series
"""
__author__ = "Liqing Peng"
__copyright__ = "Copyright (C) 2024 Liqing Peng"
__license__ = "MIT"
__version__ = "2024.05"

from pylab import *
import pandas as pd
import netCDF4 as netcdf
import os

###########################################################################
# Lat-Lon toolkits
###########################################################################
def Set_dims(res, nlat, nlon, minlat, minlon):

	dims = {}
	dims['minlat'] = minlat; dims['minlon'] = minlon
	dims['nlat'] = nlat; dims['nlon'] = nlon
	dims['res'] = res
	dims['maxlat'] = dims['minlat'] + dims['res'] * (dims['nlat'] - 1)
	dims['maxlon'] = dims['minlon'] + dims['res'] * (dims['nlon'] - 1)
	dims['undef'] = -9.99e+08

	return dims


def Create_Write_NETCDF_File_datetime(dims, file, vars, varname, data, tinitial, nt, datetimes):
	"The difference is this function will create and write the data at one time"
	nlat = dims['nlat']
	nlon = dims['nlon']
	res = dims['res']
	minlon = dims['minlon']
	minlat = dims['minlat']
	undef = dims['undef']

	# Prepare the netcdf file
	# Create file
	f = netcdf.Dataset(file, 'w', format='NETCDF4')

	# Define dimensions
	f.createDimension('lon', nlon)
	f.createDimension('lat', nlat)
	f.createDimension('time', nt)

	# Longitude
	f.createVariable('lon', 'd', ('lon',))
	f.variables['lon'][:] = np.linspace(minlon, minlon+res*(nlon-1), nlon)
	f.variables['lon'].units = 'degrees_east'
	f.variables['lon'].long_name = 'Longitude'
	f.variables['lon'].res = res

	# Latitude
	f.createVariable('lat', 'd', ('lat',))
	f.variables['lat'][:] = np.linspace(minlat, minlat+res*(nlat-1), nlat)
	f.variables['lat'].units = 'degrees_north'
	f.variables['lat'].long_name = 'Latitude'
	f.variables['lat'].res = res

	# Time
	f.createVariable('time', 'f4', ('time', ))
	f.variables['time'].units = 'hours since %04d-%02d-%02d %02d:00:00.0' % (tinitial.year,tinitial.month,tinitial.day,tinitial.hour)
	f.variables['time'].calendar = 'gregorian'
	f.variables['time'][:] = array([netcdf.date2num(dt, f.variables['time'].units, f.variables['time'].calendar) for dt in datetimes])
	f.variables['time'].long_name = 'Time'

	# Data
	if type(vars) is str:
		datafield = f.createVariable(vars, 'f', ('time', 'lat', 'lon',), fill_value=undef, zlib=True)
		f.variables[vars].long_name = varname
		datafield[:] = data

	else:
		for v, var in enumerate(vars):
			datafield = f.createVariable(var, 'f', ('time', 'lat', 'lon',), fill_value=undef, zlib=True)
			f.variables[var].long_name = varname[v]
			datafield[:] = data[v, :, :]

	f.sync()
	f.close()

	return f



def Create_NETCDF_File(dims, file, vars, varname, tinitial, freq, unit):
	nlat = dims['nlat']
	nlon = dims['nlon']
	res = dims['res']
	minlon = dims['minlon']
	minlat = dims['minlat']
	undef = dims['undef']

	# Prepare the netcdf file
	# Create file
	f = netcdf.Dataset(file, 'w', format='NETCDF4')

	# Define dimensions
	f.createDimension('lon', nlon)
	f.createDimension('lat', nlat)
	f.createDimension('t', None)

	# Longitude
	f.createVariable('lon', 'd', ('lon',))
	f.variables['lon'][:] = np.linspace(minlon, minlon + res * (nlon - 1), nlon)
	f.variables['lon'].units = 'degrees_east'
	f.variables['lon'].long_name = 'Longitude'
	f.variables['lon'].res = res

	# Latitude
	f.createVariable('lat', 'd', ('lat',))
	f.variables['lat'][:] = np.linspace(minlat, minlat + res * (nlat - 1), nlat)
	f.variables['lat'].units = 'degrees_north'
	f.variables['lat'].long_name = 'Latitude'
	f.variables['lat'].res = res

	# Time
	f.createVariable('t', 'f4', ('t', ))
	f.variables['t'].units = '%s since %04d-%02d-%02d %02d:00:00.0' % (
		freq, tinitial.year, tinitial.month, tinitial.day, tinitial.hour)
	f.variables['t'].calendar = 'gregorian'
	f.variables['t'].long_name = 'Time'

	# Data
	if type(vars) is str:
		f.createVariable(vars, 'f', ('t', 'lat', 'lon',), fill_value=undef, zlib=True)
		f.variables[vars].long_name = varname
		f.variables[vars].units = unit

	else:
		for v, var in enumerate(vars):
			f.createVariable(var, 'f', ('t', 'lat', 'lon',), fill_value=undef, zlib=True)
			f.variables[var].long_name = varname[v]

	f.sync()
	f.close()

	return f


def Write_NETCDF_File_Step(file, vars, newdata, time):
	# Open the existing file
	f = netcdf.Dataset(file, 'a')
	date_time = f.variables['t']
	posCnt = len(date_time)
	date_time[posCnt] = netcdf.date2num(time, date_time.units, date_time.calendar)

	if type(vars) is str:
		f.variables[vars][posCnt, :, :] = newdata
	else:
		for v, var in enumerate(vars):
			f.variables[var][posCnt, :, :] = newdata[v, :, :]

	f.sync()
	f.close()

	return f


def Write_NETCDF_File(dims, file, vars, varname, data, time, freq, unit):
	if os.path.isfile(file) is False:
		Create_NETCDF_File(dims, file, vars, varname, time, freq, unit)
	Write_NETCDF_File_Step(file, vars, data, time)

	return



def Create_Write_NETCDF_File(dims, file, vars, varname, data, tinitial, nt, freq):
	"The difference is this function will create and write the data at one time"
	nlat = dims['nlat']
	nlon = dims['nlon']
	res = dims['res']
	minlon = dims['minlon']
	minlat = dims['minlat']
	undef = dims['undef']
	t = np.arange(0, nt)

	# Prepare the netcdf file
	# Create file
	f = netcdf.Dataset(file, 'w', format='NETCDF4')

	# Define dimensions
	f.createDimension('lon', nlon)
	f.createDimension('lat', nlat)
	f.createDimension('t', nt)

	# Longitude
	f.createVariable('lon', 'd', ('lon',))
	f.variables['lon'][:] = np.linspace(minlon, minlon+res*(nlon-1), nlon)
	f.variables['lon'].units = 'degrees_east'
	f.variables['lon'].long_name = 'Longitude'
	f.variables['lon'].res = res

	# Latitude
	f.createVariable('lat', 'd', ('lat',))
	f.variables['lat'][:] = np.linspace(minlat, minlat+res*(nlat-1), nlat)
	f.variables['lat'].units = 'degrees_north'
	f.variables['lat'].long_name = 'Latitude'
	f.variables['lat'].res = res

	# Time
	f.createVariable('t', 'f4', ('t', ))
	f.variables['t'][:] = t
	f.variables['t'].units = '%s since %04d-%02d-%02d %02d:00:00.0' % (freq, tinitial.year,tinitial.month,tinitial.day,tinitial.hour)
	f.variables['t'].long_name = 'Time'

	# Data
	if type(vars) is str:
		datafield = f.createVariable(vars, 'f', ('t', 'lat', 'lon',), fill_value=undef, zlib=True)
		f.variables[vars].long_name = varname
		datafield[:] = data

	else:
		for v, var in enumerate(vars):
			datafield = f.createVariable(var, 'f', ('t', 'lat', 'lon',), fill_value=undef, zlib=True)
			f.variables[var].long_name = varname[v]
			datafield[:] = data[v, :, :]

	f.sync()
	f.close()

	return f



def Create_pairdata_NETCDF_File(file, vars, varname, data, tinitial, tstep):

	nt = data.shape[1]
	undef = -9.99e+08

	t = np.arange(0, nt)

	# Prepare the netcdf file
	# Create file
	f = netcdf.Dataset(file, 'w', format='NETCDF4')

	# Define dimensions
	f.createDimension('time', len(t))

	# Time
	f.createVariable('t', 'd', ('time', ))
	f.variables['t'][:] = t
	f.variables['t'].units = '%s since %04d-%02d-%02d %02d:00:00.0' % (tstep,tinitial.year,tinitial.month,tinitial.day,tinitial.hour)
	f.variables['t'].long_name = 'Time'


	# Data
	if type(vars) is str:
		datafield = f.createVariable(vars, 'f', ('time',), fill_value=undef, zlib=True)
		f.variables[vars].long_name = varname
		datafield[:] = data

	else:
		for v, var in enumerate(vars):
			datafield = f.createVariable(var, 'f', ('time', ), fill_value=undef, zlib=True)
			f.variables[var].long_name = varname[v]
			datafield[:] = data[v, :]


	f.sync()
	f.close()

	return f


from scipy import stats
def Trend_Slope_Intercept_KendallTau(DATA, mask):

	""" Use Theil-Sen slope to calculate the median (nonparametric) trend """
	tstep, glat, glon = DATA.shape
	x = np.arange(tstep, dtype=float)
	stat = np.empty((3, glat, glon))
	stat.fill(np.nan)

	for i in xrange(0, glat):
		print 'lat ', i
		for j in xrange(0, glon):
			if np.isnan(mask[i, j])==False:
				stat[0, i, j] = stats.mstats.theilslopes(DATA[:, i, j], alpha=0.95)[0]  # Recall that mstats use alpha=0.95 but stats use alpha=0.05
				stat[1, i, j] = mean(DATA[:, i, j]) - stat[0, i, j] * mean(x)
				stat[2, i, j] = stats.mstats.kendalltau(np.arange(1, tstep+1), DATA[:, i, j])[1]

	return stat

###########################################################################
# Time series toolkits (migrated from ETmax_new_analysis.py)
###########################################################################
def mon24season(mon, months):
	"return the seasonal data"

	MAM = mon[((months>=3) & (months<=5)), :]
	JJA = mon[((months>=6) & (months<=8)), :]
	SON = mon[((months>=9) & (months<=11)), :]
	DJF = mon[((months>=12) | (months<=2)), :]
	ALL = mon

	return [ALL, MAM, JJA, SON, DJF]

def mon2growseason(mon, styr, edyr):
	"return the seasonal data"
	dates = pd.date_range((str(styr)+'-01-01'), (str(edyr)+'-12-31'), freq='M')
	ddays = pd.to_datetime(dates); dyears = ddays.year; dmonths = ddays.month

	gs = mon[((dmonths>=4) & (dmonths<=10)), :]

	return gs