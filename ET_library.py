#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This library provides formulas for different PET/Ga/Gs methods
"""
__author__ = "Liqing Peng"
__copyright__ = "Copyright (C) 2024 Liqing Peng"
__license__ = "MIT"
__version__ = "2024.05"

# import numpy as np
from pylab import *

############################# Unit conversion ###############################

def Convert_Unit_Temp(input):
	"# kelvin to degree"
	data = input - 273.16
	return data

def Convert_Unit_Rad(input):
	"# W/m2 to MJ/m2/day"
	watt2jule = 10e5/86400.0
	data = input / float(watt2jule)
	return data

def Convert_2m_Wind(input):
	"""Wind: m/s"""
	z_forcing = 10  # 10 m wind field
	return (4.87 / (np.log(67.8 * z_forcing - 5.42))) * input


class PET:

	def __init__(self, INPUT, vegetation_parameter, albedo_mode='CLM', lai_mode='CLM'):

		self.Tmax = Convert_Unit_Temp(INPUT['Tmax'])  # Temperature (Celcius)
		self.Tmin = Convert_Unit_Temp(INPUT['Tmin'])  # Temperature (Celcius)
		self.Tavg = Convert_Unit_Temp(INPUT['Tavg'])  # Temperature (Celcius)
		self.RH = INPUT['RHum']    # Relative humidity (%)
		self.Wind = INPUT['Wind']    # Wind speed (m/s)
		self.Pres = INPUT['Pres']  # Pressure (Pa)
		self.LWin = INPUT['LWin']  # Radiation (W/m2)
		self.SWin = INPUT['SWin']  # Radiation (W/m2)

		# Vegetation dynamics or climatology
		if albedo_mode == 'DYN':
			self.Albedo = INPUT['Albedo']
		elif albedo_mode == 'CLM':
			self.Albedo = vegetation_parameter['Albedo_clm']
		else:
			self.Albedo = albedo_mode

		# Vegetation dynamics or climatology
		if lai_mode == 'DYN':
			self.LAI = INPUT['LAI']
			self.LAI[self.LAI<0] = 0.0
		else:
			self.LAI = vegetation_parameter['LAI_clm']

		self.landcover = vegetation_parameter['Landcover']
		self.canopy_height_min = vegetation_parameter['chmin']
		self.canopy_height_max = vegetation_parameter['chmax']
		self.d0p = vegetation_parameter['d0p']
		self.z0m = vegetation_parameter['z0m']
		self.z0h = vegetation_parameter['z0h']
		self.z0g = vegetation_parameter['z0g']
		self.gstmax = vegetation_parameter['gstmax'] / 1000.0  # m/s
		self.rstmin = vegetation_parameter['rstmin']  # s/m
		self.LAImax = vegetation_parameter['LAImax']
		self.kB = vegetation_parameter['kB']
		self.z_forcing = 10.0

		# Run the basic calculations
		self.Calculate_Rnet()
		self.Calculate_Saturated_Vapor_Pressure()
		self.Calculate_Slope_Saturation_Vapor_Pressure_Curve()
		self.Calculate_VPD()
		self.Calculate_Density_Moist_Air()
		self.Calculate_Psychrometric_Constant()

	############################ Meteorological parameters ############################
	def Calculate_Rnet(self):
		"Rnet"
		# stefan_b = 4.903e-9  # [MJ K-4 m-2 day-1]
		stefan_b = 5.670367e-8  # [W K-4 m-2]
		self.SWout = self.Albedo * self.SWin
		self.LWout = stefan_b * ((self.Tmax + 273.16) ** 4 + (self.Tmin + 273.16) ** 4) / 2.0
		self.Rn = self.SWin - self.SWout + self.LWin - self.LWout  # (W/m2)

		return

	def Calculate_Saturated_Vapor_Pressure(self):
		"DELTA"
		# Because the relationship between saturated vapor pressure and temperature is not linear, using Tmax and Tmin is preferable to estimating estar
		estar_max = 0.6108 * np.exp((17.27 * self.Tmax) / (237.3 + self.Tmax))   # (kPa)
		estar_min = 0.6108 * np.exp((17.27 * self.Tmin) / (237.3 + self.Tmin))   # (kPa)
		self.estar = (estar_max + estar_min) / 2 * 1000.0   # (Pa)

		return

	def Calculate_Slope_Saturation_Vapor_Pressure_Curve(self):

		DELTA_max = 4098 * 0.6108 * np.exp((17.27 * self.Tmax) / (237.3 + self.Tmax)) / (237.3 + self.Tmax) ** 2  # (kPa/degC)
		DELTA_min = 4098 * 0.6108 * np.exp((17.27 * self.Tmin) / (237.3 + self.Tmin)) / (237.3 + self.Tmin) ** 2  # (kPa/degC)
		self.DELTA = (DELTA_max + DELTA_min) / 2    # (kPa/degC)

		return

	def Calculate_VPD(self):
		"Humidity related output"
		self.e = self.estar * self.RH / 100.0  # (Pa)
		self.VPD = self.estar - self.e    	   # (Pa)

		# q = r/(1+r)
		# r = 0.622e/(Pres-e)
		self.eps = 0.622
		self.SHum = self.eps * self.e / (self.Pres - (1-self.eps) * self.e)   # (kg/kg)

		return

	def Calculate_Density_Moist_Air(self):
		"Temperature related"
		Rs = 287.058  		#(J kg-1 K-1)
		# Calculate Virtual Temperature
		self.Tv = (self.Tavg + 273.16) / (1 - (1 - self.eps) * self.e / self.Pres)
		self.rhoa = (self.Pres + self.e) / self.Tv / Rs    # (kg/m3)
		# Pa = J/m3

		return

	def Calculate_Psychrometric_Constant(self):
		"Gamma"
		self.lv = 2.501 - 0.002361 * self.Tavg         # (MJ/kg)
		self.cs = 1.005 + 1.82 * self.SHum             # (kJ/kg/K)
		# kJ/kg/K * Pa *kg /MJ = Pa/K/k
		self.gamma = self.cs * self.Pres / self.eps / self.lv * 10 ** (-6)   # (kPa/K)

		return

	def Calculate_Vegetation_Fraction(self):

		self.fveg = 1 - np.exp(-0.5 * self.LAI)

		return

	###################### Bigleaf: Aerodynamic Conductance ######################

	def ga_roughness(self):
		"""
		# Required unit, converted from the calculation
		Wind: m/s
		Z0m
		Z0h
		d0
		:return: ga (m/s)
		"""
		Zm, Zh = self.z_forcing, self.z_forcing
		kappa = 0.4
		log_term = np.log((Zm - self.d0p)/self.z0m) * np.log((Zh - self.d0p)/self.z0h)
		wind_term = kappa ** 2 * self.Wind
		ga = wind_term / log_term

		return ga

	def ga_canopy_height(self):
		self.Calculate_canopy_height()
		Zm = self.canopy_height + 2
		Zh = self.canopy_height + 2
		kappa = 0.4
		d0p = self.canopy_height * 2.0 / 3.0
		z0m = self.canopy_height / 8.0
		z0h = z0m / np.exp(self.kB)
		log_term = np.log((Zm - d0p) / z0m) * np.log((Zh - d0p) / z0h)

		# Zhou 2006, Brutsaert, 1982 p59,
		# Weather station roughness length Federer, 1996
		z0w = 0.005
		# Fetch at weather station
		Fw = 5000.0
		# Height of internal boundary layer
		z_boundary = 0.334 * (Fw ** 0.875) * (z0w ** 0.125)
		# Wind speed at the reference height
		u_ref = self.Wind * np.log(z_boundary / z0w) / np.log(z_boundary / z0m) * np.log(
			(Zm - d0p) / z0m) / np.log(self.z_forcing / z0w)
		wind_term = kappa ** 2 * u_ref
		ga = wind_term / log_term

		return ga

	def ga_open_water(self):
		"""
		Open water aerodynamic resistance equation
		# Required unit, converted from the calculation
		Pres: Pa
		Wind: m/s
		rhoa: kg/m3
		cs: J/kg/K
		VPD: kPa
		ga: m/s
		lv: J/kg
		:return: ga (m/s)
		"""
		# Convert 10m NLDAS-2 wind to 2m wind
		wind2m = Convert_2m_Wind(self.Wind)
		# self.ra = self.rhoa * (1000.0) * 86400 / self.Pres * self.eps * self.lv / (6.43 * (1 + 0.536 * wind2m))
		ga = 6.43 * (1 + 0.536 * wind2m) * self.Pres / (86400/1000.0) / self.rhoa / self.eps / self.lv / 1e6  # convert MJ/kg to J/kg
		# ga filter
		# ga[(ga < 0)|(ga > 3.0/20.0)] = np.nan

		return ga

	def ga_reference_crop(self):
		"""
		FAO reference crop aerodynamic resistance equation
		# Required unit, converted from the calculation
		Wind: m/s
		:return: ga (m/s)
		"""
		# Convert 10m NLDAS-2 wind to 2m wind
		wind2m = Convert_2m_Wind(self.Wind)
		ga = wind2m / 208.0
		# ga filter
		# ga[(ga < 0)|(ga > 3.0/20.0)] = np.nan

		return ga

	########################### Bigleaf: Surface Conductance ######################

	def gsmax_K1995(self):

		LAI = np.copy(self.LAI)
		# give a buffer range for 3 in Kelliher 1995
		# The Shuttle Worth use LAI=4 for closed canopy threshold
		LAI[LAI > 4] = 4

		gsmax = self.gstmax * LAI

		# set the nonvegetated grid (at a specific time of the year or location) to nan
		self.nonveg_index = (LAI == 0) | (self.landcover == 0) | (self.landcover == 13) | (self.landcover == 15) | (self.landcover == 16)
		gsmax[self.nonveg_index] = np.nan

		return gsmax

	def gsmax_Yan2012(self):

		LAI = np.copy(self.LAI)
		# give a buffer range for 3 in Kelliher 1995
		# The Shuttle Worth use LAI=4 for closed canopy threshold
		LAI[LAI > 4] = 4
		# rst_min to Canopy conductance

		gsmax = 1.0/self.rstmin * LAI

		# set the nonvegetated grid (at a specific time of the year or location) to nan
		self.nonveg_index = (LAI == 0) | (self.landcover == 0) | (self.landcover == 13) | (self.landcover == 15) | (self.landcover == 16)
		gsmax[self.nonveg_index] = np.nan

		return gsmax

	def gsmax_Zhou2006(self):
		# Zhou 2006
		LAI_effective = np.copy(self.LAI)
		LAI_effective[(LAI_effective > 2) & (LAI_effective <= 4)] = LAI_effective[(LAI_effective > 2) & (LAI_effective <= 4)] / 2
		LAI_effective[LAI_effective > 4] = 2
		# rst_min to Canopy conductance
		gsmax = 1.0 / self.rstmin * LAI_effective

		# set the nonvegetated grid (at a specific time of the year or location) to nan
		self.nonveg_index = (LAI_effective == 0) | (self.landcover == 0) | (self.landcover == 13) | (self.landcover == 15) | (self.landcover == 16)
		gsmax[self.nonveg_index] = np.nan

		return gsmax

	def gs_reference_crop_short(self):
		# Make sure it is float
		gs = 1.0/70.0
		return gs

	########################### Two Source SW Parameters ######################

	def Calculate_canopy_height(self):
		self.canopy_height = self.canopy_height_min + (self.canopy_height_max - self.canopy_height_min) * self.LAI/self.LAImax
		self.canopy_height[(self.LAImax == 0) | (self.canopy_height <= 0)] = 0

	def Calculate_height_parameters(self):
		self.Calculate_canopy_height()
		self.z_ref = self.canopy_height + 2
		z0c = np.empty(self.canopy_height.shape); z0c.fill(np.nan)
		Cd = np.empty(self.canopy_height.shape); Cd.fill(np.nan)
		# n: eddy diffusivity decay constant of vegetation
		self.n = np.empty(self.canopy_height.shape); self.n.fill(np.nan)
		self.d0 = np.empty(self.canopy_height.shape); self.d0.fill(np.nan)
		self.z0 = np.empty(self.canopy_height.shape); self.z0.fill(np.nan)

		z0c[self.canopy_height<=1] = 0.13 * self.canopy_height[self.canopy_height<=1]
		self.n[self.canopy_height<=1] = 2.5

		z0c[(self.canopy_height > 1) & (self.canopy_height < 10)] = 0.139 * self.canopy_height[(self.canopy_height > 1) & (self.canopy_height < 10)] - 0.009 * self.canopy_height[(self.canopy_height > 1) & (self.canopy_height < 10)] ** 2
		self.n[(self.canopy_height > 1) & (self.canopy_height < 10)] = 2.306 + 0.194 * self.canopy_height[(self.canopy_height > 1) & (self.canopy_height < 10)]

		z0c[self.canopy_height >= 10] = 0.05 * self.canopy_height[self.canopy_height >= 10]
		self.n[self.canopy_height >= 10] = 4.25

		Cd[self.canopy_height == 0] = 1.4 * 10 ** (-3)
		Cd[self.canopy_height > 0] = (-1 + np.exp(0.909 - 3.03 * z0c[self.canopy_height > 0] / self.canopy_height[self.canopy_height > 0])) ** 4 / 4.0

		self.d0[self.LAI >= 4] = self.canopy_height[self.LAI >= 4] - z0c[self.LAI >= 4] / 0.3
		self.d0[self.LAI < 4] = 1.1 * self.canopy_height[self.LAI < 4] * np.log(1 + (Cd[self.LAI < 4] * self.LAI[self.LAI < 4]) ** 0.25)

		z01 = 0.3 * (self.canopy_height - self.d0)
		z02 = self.z0g + 0.3 * self.canopy_height * (Cd * self.LAI) ** 0.5
		self.z0 = np.minimum(z01, z02)
		print(self.z0.shape)


	def Convert_Wind_reference_height(self):
		"The wind speed observed at the weather stations is converted to the reference height using a logarithmic profile in that the internal boundary layer heights over the weather ground and canopy surface are matched and a step change in surface roughness from z0 to z0w is assumed"
		# Zhou 2006, Brutsaert, 1982 p59,
		# Weather station roughness length Federer, 1996
		z0w = 0.005
		# Fetch at weather station
		Fw = 5000.0
		# Height of internal boundary layer
		z_boundary = 0.334 * (Fw ** 0.875) * (z0w ** 0.125)
		# Wind speed at the reference height
		u_ref = self.Wind * np.log(z_boundary/z0w) / np.log(z_boundary/self.z0) * np.log((self.z_ref - self.d0)/self.z0) / np.log(self.z_forcing/z0w)
		return u_ref

	def Calculate_aerodynamic_resistances_between_surface_reference(self):
		kappa = 0.4
		self.Calculate_height_parameters()
		self.u_ref = self.Convert_Wind_reference_height()
		ustar = kappa * self.u_ref / np.log((self.z_ref - self.d0)/self.z0)
		Kh = kappa * ustar * (self.canopy_height - self.d0)
		# 1. aerodynamic resistance leaving the surface before being incorporated in the canopy (s/m)
		self.r_s_a = self.canopy_height * np.exp(self.n) / self.n / Kh * (np.exp(- self.n * self.z0g / self.canopy_height) - np.exp(- self.n * (0.13*self.canopy_height + 0.63 * self.canopy_height)/self.canopy_height))

		# 2. aerodynamic resistance between the canopy height to the reference height (s/m)
		self.r_a_a = 1 / kappa / ustar * np.log((self.z_ref - self.d0)/(self.canopy_height - self.d0)) + self.canopy_height / self.n / Kh * (np.exp(self.n * (1 - (0.13 * self.canopy_height + 0.63 * self.canopy_height)/self.canopy_height)) - 1)


	def Calculate_canopy_bulk_boundary_layer_resistance(self):
		# S & W 1985, Brisson 1998, 50 s/m
		R_boundary = 50
		shielding_factor_boundary = 0.5
		# Bulk boundary layer resistance of the vegetative elements in the canopy
		self.r_c_a = R_boundary * shielding_factor_boundary / self.LAI


	def Calculate_canopy_bulk_stomatal_resistance(self):
		# Use rst_min from the look up table
		# shielding_factor_canopy = 0.5
		LAI_effective = np.copy(self.LAI)
		LAI_effective[(LAI_effective > 2) & (LAI_effective <=4)] = LAI_effective[(LAI_effective > 2) & (LAI_effective <=4)]/2
		LAI_effective[LAI_effective > 4] = 2

		# self.r_c_s = self.rstmin * shielding_factor_canopy / LAI_effective
		self.r_c_s = self.rstmin / LAI_effective
		# Tourula and Heikinheimo, 1998 limit, from Zhou 2006
		self.r_c_s[self.r_c_s > 50000] = 50000
		# set the nonvegetated grid (at a specific time of the year or location) to nan
		self.nonveg_index = (LAI_effective == 0) | (self.landcover == 0) | (self.landcover == 13) | (self.landcover == 15) | (self.landcover == 16)
		self.r_c_s[self.nonveg_index] = np.nan


	############################# PET algorithms ############################

	def Penman_Monteith_Big_Leaf(self, ga_mode='OW', gs_mode='OW'):
		"""
		Classical Penman-Monteith equation for big leaf model
		# Required unit, converted from the calculation
		DELTA: kPa/K
		gamma: kPa/K
		Rn: W/m2
		rhoa: kg/m3
		cs: J/kg/K
		VPD: kPa
		ga: m/s
		lv: J/kg
		:return: Penman-PET
		"""
		if ga_mode == 'OW':
			ga = self.ga_open_water()
		elif ga_mode == 'FAO':
			ga = self.ga_reference_crop()
		elif ga_mode == 'LC':
			ga = self.ga_roughness()
		elif ga_mode == 'CH':
			ga = self.ga_canopy_height()
		else:
			print('ga_mode missing!')
			exit()

		epsilon = self.DELTA/self.gamma
		nominator = epsilon * self.Rn + self.rhoa * (self.cs * 1000.0) / self.gamma * (self.VPD / 1000.0) * ga

		# Big leaf model gs = gc
		if gs_mode == 'OW':
			denominator =  epsilon + 1
		elif gs_mode == 'FAO':
			gs = self.gs_reference_crop_short()
			denominator = epsilon + 1 + ga/gs
		elif gs_mode == 'K1995':
			gs = self.gsmax_K1995()
			denominator = epsilon + 1 + ga / gs
		elif gs_mode == 'Yan2012':
			gs = self.gsmax_Yan2012()
			denominator = epsilon + 1 + ga/gs
		elif gs_mode == 'Zhou2006':
			gs = self.gsmax_Zhou2006()
			denominator = epsilon + 1 + ga/gs
		else:
			print('gs_mode missing!')
			exit()

		# PET: kg/m2/s, convert to mm/d, * 86400
		ETp_PM = nominator / denominator / (self.lv * 1e6) * 86400
		ETp_PM[ETp_PM<0.0] = 0.0

		return ETp_PM


	def Shuttleworth_Wallace_Canopy_Soil(self):
		"""
		Shuttleworth - Wallace equation
		# Required unit, converted from the calculation
		DELTA: kPa/K
		gamma: kPa/K
		Rn: W/m2
		rhoa: kg/m3
		cs: J/kg/K
		VPD: kPa
		ga: m/s
		lv: J/kg
		:return: Penman-PET
		"""
		# Get FVEG
		self.Calculate_Vegetation_Fraction()

		# set up the resistances
		# r^s_a - aerodynamic resistance leaving the surface before being incorporated in the canopy (s/m)
		# r^a_a - aerodynamic resistance between the canopy height to the reference height (s/m)
		self.Calculate_aerodynamic_resistances_between_surface_reference()
		# r^c_a - bulk boundary layer resistance of the vegetative elements in the canopy (s/m)
		self.Calculate_canopy_bulk_boundary_layer_resistance()
		# r^c_s - bulk stomatal resistance for the vegetation, well known (s/m)
		self.Calculate_canopy_bulk_stomatal_resistance()
		# r^s_s - surface resistance of the subtrate (s/m). If it is saturated surface, leave it to zero
		r_s_s = 0.0

		# For ETp canopy
		ETp_canopy_nominator = self.DELTA * self.Rn + (self.rhoa * (self.cs * 1000.0) * (self.VPD / 1000.0) - self.DELTA * self.r_c_a * self.Rn * (1 - self.fveg)) / (self.r_a_a + self.r_c_a)
		ETp_canopy_denominator = self.DELTA + self.gamma * (1 + self.r_c_s / (self.r_a_a + self.r_c_a))

		# PET: kg/m2/s, convert to mm/d, * 86400
		ETp_canopy = ETp_canopy_nominator / ETp_canopy_denominator / (self.lv * 1e6) * 86400

		# For ETp soil
		ETp_soil_nominator = self.DELTA * self.Rn + (self.rhoa * (self.cs * 1000.0) * (self.VPD / 1000.0) - self.DELTA * self.r_c_a * self.Rn * self.fveg) / (self.r_a_a + self.r_s_a)
		ETp_soil_denominator = self.DELTA + self.gamma * (1 + r_s_s / (self.r_a_a + self.r_s_a))
		ETp_soil = ETp_soil_nominator / ETp_soil_denominator / (self.lv * 1e6) * 86400

		# For Cc and Cs
		R_a = (self.DELTA + self.gamma) * self.r_a_a
		R_s = (self.DELTA + self.gamma) * self.r_s_a + self.gamma * r_s_s
		R_c = (self.DELTA + self.gamma) * self.r_c_a + self.gamma * self.r_c_s

		C_canopy = 1.0 / (1 + R_c * R_a / R_s / (R_c + R_a))
		C_soil = 1.0 / (1 + R_s * R_a / R_c / (R_s + R_a))

		# W/m2 -> kg/s/m2
		ETp_SW = C_canopy * ETp_canopy + C_soil * ETp_soil

		# ETp_PM[ETp_PM<0.0] = 0.0

		return ETp_SW



	##### Penman and open-water
	def Penman(self):
		"""
		Penman equation
		# Required unit, converted from the calculation
		DELTA: kPa/K
		gamma: kPa/K
		Rn: W/m2
		rhoa: kg/m3
		cs: J/kg/K
		VPD: kPa
		ga: m/s
		lv: J/kg
		:return: Penman-PET
		"""

		ga = self.ga_roughness()
		energy_term = self.DELTA * self.Rn
		aero_term = self.rhoa * (self.cs * 1000.0) * (self.VPD / 1000.0) * ga
		denominator = (self.lv * 1e6) * (self.DELTA + self.gamma)
		# PET: kg/m2/s, convert to mm/d, * 86400
		ETp_Penman = (energy_term + aero_term) / denominator * 86400
		ETp_Penman[ETp_Penman<0.0] = 0.0

		return ETp_Penman


	def Penman_openwater_shuttleworth(self):
		"""
		This is the original Penman openwater equation from Shuttleworth, 2013
		# Required unit, converted from the calculation
		DELTA: kPa/K
		gamma: kPa/K
		Rn: MJ/m2/day
		rhoa: kg/m3
		cs: J/kg/K
		VPD: kPa
		ga: m/s
		lv: MJ/kg
		:return: Penman-openwater-PET
		"""
		A = Convert_Unit_Rad(self.Rn)
		wind2m = Convert_2m_Wind(self.Wind)
		PET_R = (self.DELTA / (self.DELTA + self.gamma)) * A / self.lv
		PET_A = (self.gamma / (self.DELTA + self.gamma)) * ((6.43 * (1 + 0.536 * wind2m) * (self.VPD / 1000.0)) / self.lv)
		ETp_ow = PET_R + PET_A
		ETp_ow[ETp_ow<0] = 0.0

		return ETp_ow

	def Penman_openwater_ga(self):
		"""
		This is based on the Penman equation, but the ga module is inverted from the Penman open water equation in Shuttleworth, see details in GCB, 2019, SI Eq.S2, S6.
		# Required unit, converted from the calculation
		DELTA: kPa/K
		gamma: kPa/K
		Rn: W/m2
		rhoa: kg/m3
		cs: J/kg/K
		VPD: kPa
		ga: m/s
		lv: J/kg
		:return: Penman-openwater-PET
		"""
		ga = self.ga_open_water()
		energy_term = self.DELTA * self.Rn
		aero_term = self.rhoa * (self.cs * 1000.0) * (self.VPD / 1000.0) * ga
		denominator = (self.lv * 1e6) * (self.DELTA + self.gamma)
		# PET: kg/m2/s, convert to mm/d, * 86400
		ETp_ow = (energy_term + aero_term) / denominator * 86400
		ETp_ow[ETp_ow<0.0] = 0.0

		return ETp_ow

	##### FAO refererence crops
	def FAO_reference_shortcrop(self):
		"""
		This is the original FAO equation
		# Required unit, converted from the calculation
		DELTA: kPa/K
		gamma: kPa/K
		Rn/A: MJ/m2/day
		Tavg: Celcius
		VPD: kPa
		wind2m: m/s
		:return: Penman-openwater-PET
		"""
		A = Convert_Unit_Rad(self.Rn)
		wind2m = Convert_2m_Wind(self.Wind)
		# 1/lambda is approximate to 1/2.45 = 0.408
		# PET_N = self.DELTA * A / self.lv + 900 * self.gamma / (self.Tavg + 273) * wind2m * (self.VPD / 1000.0)
		PET_N = 0.408 * self.DELTA * A + 900 * self.gamma / (self.Tavg + 273) * wind2m * (self.VPD / 1000.0)
		# Set rs = 0 s/m or gs = inifity m/s
		# PET_D = self.DELTA + self.gamma * 1.0 # Then this will be the same as the PET FAO ga
		PET_D = self.DELTA + self.gamma * (1 + 0.34 * wind2m)

		ETp_FAO = PET_N / PET_D
		ETp_FAO[ETp_FAO<0] = 0.0

		return ETp_FAO


	def FAO_reference_tallcrop(self):
		"""
		This is the original FAO equation
		# Required unit, converted from the calculation
		DELTA: kPa/K
		gamma: kPa/K
		Rn/A: MJ/m2/day
		Tavg: Celcius
		VPD: kPa
		wind2m: m/s
		:return: Penman-openwater-PET
		"""
		A = Convert_Unit_Rad(self.Rn)
		wind2m = Convert_2m_Wind(self.Wind)
		# 1/lambda is approximate to 1/2.45 = 0.408
		# PET_N = self.DELTA * A / self.lv + 900 * self.gamma / (self.Tavg + 273) * wind2m * (self.VPD / 1000.0)
		PET_N = 0.408 * self.DELTA * A + 1600 * self.gamma / (self.Tavg + 273) * wind2m * (self.VPD / 1000.0)
		# Set rs = 0 s/m or gs = inifity m/s
		# PET_D = self.DELTA + self.gamma * 1.0 # Then this will be the same as the PET FAO ga
		PET_D = self.DELTA + self.gamma * (1 + 0.38 * wind2m)

		ETp_FAO = PET_N / PET_D
		ETp_FAO[ETp_FAO<0] = 0.0

		return ETp_FAO


	def Penman_FAO_ga(self):
		"""
		This is based on the Penman equation, but the ga module is inverted from the FAO reference crop formulation. Note that there is no Rs/Gs in the Penman equation, which means this is fundamental different from the FAO, which has Rs=70 s/m.
		# Required unit, converted from the calculation
		DELTA: kPa/K
		gamma: kPa/K
		Rn: W/m2
		rhoa: kg/m3
		cs: J/kg/K
		VPD: kPa
		ga: m/s
		lv: J/kg
		:return: Penman-reference crop-PET
		"""
		ga = self.ga_reference_crop()
		energy_term = self.DELTA * self.Rn
		aero_term = self.rhoa * (self.cs * 1000.0) * (self.VPD / 1000.0) * ga
		denominator = (self.lv * 1e6) * (self.DELTA + self.gamma)
		# PET: kg/m2/s, convert to mm/d, * 86400
		ETp_FAO = (energy_term + aero_term) / denominator * 86400
		ETp_FAO[ETp_FAO<0.0] = 0.0

		return ETp_FAO


	def Priestley_Taylor(self):
		"""
		This is the original P-T equation
		# Required unit, converted from the calculation
		DELTA: kPa/K
		gamma: kPa/K
		Rn/A: MJ/m2/day
		lv: MJ/kg
		:return: Penman-openwater-PET
		"""
		alpha = 1.26
		A = Convert_Unit_Rad(self.Rn)
		ETp_PT = alpha * self.DELTA / (self.DELTA + self.gamma) * A / self.lv
		ETp_PT[ETp_PT<0.0] = 0.0

		return ETp_PT


	def Equilibrium_Evaporation(self):
		"""
		This is the equilibrium evaoporation
		# Required unit, converted from the calculation
		DELTA: kPa/K
		gamma: kPa/K
		Rn/A: MJ/m2/day
		lv: MJ/kg
		:return: Penman-openwater-PET
		"""
		A = Convert_Unit_Rad(self.Rn)
		ETp_eq = self.DELTA / (self.DELTA + self.gamma) * A / self.lv
		ETp_eq[ETp_eq<0.0] = 0.0

		return ETp_eq

	def Control_Potential_Evaporation(self):
		"""
		This is the control scenario, 0
		:return: Penman-control
		"""
		ETp_control = np.zeros((self.Rn).shape)

		return ETp_control
