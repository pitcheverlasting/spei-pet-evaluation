#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script provides plotting functions for correlation bar plots
"""
__author__ = "Liqing Peng"
__copyright__ = "Copyright (C) 2024 Liqing Peng"
__license__ = "MIT"
__version__ = "2024.05"

from pylab import *
from netCDF4 import Dataset
import pandas as pd
import matplotlib.cbook as cbook
import Plotting
import itertools, pickle

wet_home = '/server/user/Data/PROJECT/WET_US'
workspace = '%s/workspace_2024' %(wet_home)  # 2021
figdir = '%s/Figure/202405' %wet_home  #202012
para_dir = '/server/user/Data/forcing/NLDAS2/parameter/CONUS'

# dims_025 = DataHandeling.Set_dims(0.25, 112, 232, 25.125, -124.875)
mk025file = '/server/user/Masks/0.25deg/us_0.25deg.nc'
mk25 = Dataset(mk025file).variables['data'][0,: :]
# set the mask missing values to nan. otherwise it will set to zero when hstacking
mk25[mk25 == 0] = np.nan

SM = 'SMsurf'
steps = [1, 3, 6, 12]
dir24 = '%s/run_2024/evaluation/spei/growingseason' % (wet_home)
dir21 = '%s/run_2021/evaluation/spei/growingseason' % (wet_home)
dir20 = '%s/run_2020/evaluation/spei/growingseason' % (wet_home)
dir19 = '%s/run_2019/evaluation/spei/growingseason' % (wet_home)

# methods = ['control', 'ow', 'pt', 'faoshort', 'BL.ga_LC.gs_K1995.alb_CLM', 'TS.SW.alb_CLM']
# dirs = [dir20, dir19, dir19, dir19, dir24,  dir24]

# methodnames = ['PET = 0', 'PET-OW', 'PET-PT', 'PET-RC', 'PET-LC', 'PET-SW']
# methodcolors = ['DodgerBlue', 'SkyBlue', 'gold', 'crimson', 'cornflowerblue', 'limegreen']
# methodcolors = ['DodgerBlue', 'cornflowerblue', 'goldenrod', 'crimson', 'green', 'LightSeaGreen']

# 2024 June
methods = ['control', 'ow', 'pt', 'faoshort', 'faotall',
		   'BL.ga_LC.gs_K1995.alb_CLM', 'BL.ga_LC.gs_Zhou2006.alb_CLM',
		   'BL.ga_CH.gs_K1995.alb_CLM', 'BL.ga_CH.gs_Zhou2006.alb_CLM',
		   'TS.SW.alb_CLM']
dirs = [dir20, dir19, dir19, dir19, dir19,
		dir21, dir21,
		dir24, dir24, dir24]
methodnames = ['PET = 0', 'OW', 'PT', 'RC-short', 'RC-tall',
			   'LC-K', 'LC-Z', 'CH-K', 'CH-Z', 'SW']
methodcolors = ['Dodgerblue', 'SkyBlue', 'Orange', 'DarkOrange', 'OrangeRed',
				'lightseagreen', 'lightseagreen', 'green', 'green', 'Sienna']

final_methods_list = [1, 2, 3, 4, 5, 7, 9]
methods_fig7_list = [4, 5, 7, 9]

lc = {'filename': 'nldas_land_cover_climatology.nc', 'varname': 'landcover', 'longname':'"MODIS global land cover climatology (0.5km to 0.125deg)', 'unit':''}


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


def energy_water_limited_threshold():
	pfile = '%s/run/yearly/pre.025.yearly.1981-2017.nc'  % (wet_home)
	petfile = '%s/run_2019/yearly/pet.ow.yearly.1981-2017.nc' % (wet_home)
	p = Dataset(pfile).variables['pr'][0, :]
	pe = Dataset(petfile).variables['pe'][0, :]
	aridity = (pe/p)

	return aridity


def prepare_correlation_array(step):
	dataplot = []
	for method, dir, methodname in zip(methods, dirs, methodnames):
		corrfile = '%s/corr.spei.%s.%s.mon%02d.1981-2017.nc' % (dir, SM, method, step)
		corr = Dataset(corrfile).variables['corr'][0, :]
		dataplot.append(read_corrdata_usmask(corr))
	dataplot = array(dataplot)

	dataplot.dump('%s/array_corr_SM_SPEI%02d' % (workspace, step))
	return

# [prepare_correlation_array(step) for step in steps]


def corr_relative_barplot_parameterization(step):
	"This is for Fig.4: the bar plot performance"
	# step = 1 # old fig is 1
	n = 8 # remove PT
	positions = arange(0, n)
	lctypes = ['EBF', 'DBF', 'MF', 'ENF', 'WSA', 'SAV', 'CSH', 'OSH', 'GRA', 'CRO', 'MOS']
	foresttypes = ['EBF', 'DBF', 'MF', 'ENF', 'WSA']
	nonforesttypes = ['SAV', 'CSH', 'OSH', 'GRA', 'CRO', 'MOS']

	# methods = ['control', 'ow', 'pt', 'faoshort', 'faotall',
	# 		   'BL.ga_LC.gs_K1995.alb_CLM', 'BL.ga_LC.gs_Zhou2006.alb_CLM',
	# 		   'BL.ga_CH.gs_K1995.alb_CLM', 'BL.ga_CH.gs_Zhou2006.alb_CLM',
	# 		   'TS.SW.alb_CLM']

	methodnames = ['OW', 'RC-short', 'RC-tall',
				   'LC-K', 'LC-Z', 'CH-K', 'CH-Z', 'SW']
	colors_left = ['skyblue', 'DarkOrange', 'DarkOrange',
				   'lightseagreen', 'lightseagreen',
				   'g', 'g',
				   'Sienna']
	colors_right = ['lightskyblue', 'Gold', 'Gold',
				   'aquamarine', 'aquamarine',
				   'lightgreen', 'lightgreen',
				   'Sandybrown']

	data_array = load('%s/array_corr_SM_SPEI%02d' %(workspace, step))
	lcdata = Dataset('%s/nldas_land_cover_climatology_0.25deg.nc' % (para_dir)).variables[lc['varname']][0, :100, :]

	def find_lc_group(corrdata, lctypes):
		index = (lcdata == igbp_colors[lctypes[0]][2])
		if len(lctypes) > 1:
			for lctype in lctypes[1:]:
				index = index | (lcdata == igbp_colors[lctype][2])

		data_all = corrdata[index]
		data_all = data_all[~isnan(data_all)]
		return data_all

	def plot_forest_type(ax, caption, xlimits, lcgroup, ylabel):
		corr_control_list = []
		corr_ow_list = []
		data_list = []

		for im in range(1,10,1):
			ow_control = data_array[1, :] - data_array[0, :]
			method_ow = data_array[im, :] - data_array[1, :]
			if im != 2: # skip PT method
				if np.median(method_ow[~isnan(method_ow)]) > 0:
					corr_control_list.append(np.mean(find_lc_group(ow_control, lcgroup)))
					corr_ow_list.append(np.mean(find_lc_group(method_ow, lcgroup)))
					data_list.append(np.mean(find_lc_group(ow_control, lcgroup)) + np.median(find_lc_group(method_ow, lcgroup)))
				else:
					method_control = data_array[im, :] - data_array[0, :]
					corr_control_list.append(np.mean(find_lc_group(method_control, lcgroup)))

					corr_ow_list.append(0)
					data_list.append(np.nan)

		ax.barh(positions, np.array(corr_control_list), align='center', color=np.array(colors_left), edgecolor=np.array(colors_left))
		ax.barh(positions, np.array(corr_ow_list), left=np.array(corr_control_list), align='center', color=np.array(colors_right), edgecolor=np.array(colors_right))

		# print out
		df = pd.DataFrame(np.array(corr_control_list)+np.array(corr_ow_list), index=methodnames)
		print(df)

		ax.scatter(data_list, positions, marker='o', facecolor='k', edgecolor='k', s=30, zorder=2)
		ax.set_ylim([-1, n])
		ax.set_xlim([xlimits[0], xlimits[1]])
		ax.set_xticks(arange(xlimits[0], xlimits[1]+0.001, 0.02))
		ax.set_xlabel('$\Delta$R', fontsize=20)
		ax.invert_yaxis()
		ax.set_yticks(positions)
		if ylabel == True:
			ax.set_yticklabels(methodnames, fontsize=16)
		else:
			ax.set_yticklabels([])
		ax.tick_params(axis='x', labelsize=16)
		ax.set_title('%s'%caption, fontsize=18)


		for side in ax.spines.keys():
			ax.spines[side].set_linewidth(1.5)
		ax.set_axisbelow(True)
		ax.xaxis.grid(which='major', color='silver', linestyle='-', linewidth=0.5)
		ax.xaxis.set_ticks_position('bottom')

		xmin, xmax = ax.get_xlim()
		ymin, ymax = ax.get_ylim()
		ax.text(0.005, ymin + 0.88 * (ymax - ymin), 'Reference', ha='left', fontsize=16)

		return ax

	fig, ax = plt.subplots(1, 3, figsize=(12, 6))

	plot_forest_type(ax[0], '(a) CONUS', (0, 0.08), lctypes, True)
	plot_forest_type(ax[1],'(b) Forested', (0, 0.08), foresttypes, False)
	plot_forest_type(ax[2],'(c) Nonforested', (0, 0.08), nonforesttypes, False)

	plt.tight_layout(rect=[0, 0, 1, 1])
	plt.subplots_adjust(bottom=0.1, left=0.1)
	plt.show()
	# plt.savefig('%s/corr_diff_8methods_SPEI%s.pdf' %(figdir, step), dpi=300)
	# plt.savefig('%s/corr_diff_8methods_SPEI%s.png' %(figdir, step), dpi=300)

# corr_relative_barplot_parameterization(1)
# [corr_relative_barplot_parameterization(step) for step in steps]
# exit()


def Map_correlation_difference_final_methods():
	"""Updated Fig. 6
	"""
	fig = plt.figure(figsize=(10, 8))
	for i, mi in enumerate(final_methods_list):
		for j, step in enumerate(steps):
			dataplot = load('%s/array_corr_SM_SPEI%02d' % (workspace, step))
			ax = fig.add_subplot(7, 4, i*4 + j + 1)
			ax.m = Plotting.setup_USmap(False, False)
			if i == 0:
				im1 = ax.m.imshow(dataplot[mi,:], vmin=0, vmax=1, cmap=Plotting.cmap_spectral_sub, alpha=0.7)
				ax.set_title('SPEI-%s' % step, fontsize=14)
			else: # 1: OW, 5: dir21
				im2 = ax.m.imshow(dataplot[mi,:]-dataplot[1,:], vmin=-0.2, vmax=0.2, cmap=Plotting.cmap_seismic_sub, alpha=0.7)

			if j == 0:
				if i == 0:
					ax.set_ylabel('%s' % methodnames[mi], labelpad=10, fontsize=12)
				elif (i == 1) | (i == 6):
					ax.set_ylabel('%s - OW' % (methodnames[mi]), labelpad=10, fontsize=12)
				else:
					ax.set_ylabel('%s\n - OW' % (methodnames[mi]), labelpad=10, fontsize=12)

	cax = fig.add_axes([0.03, 0.03, 0.38, 0.02])
	cb1 = fig.colorbar(im1, cax, orientation='horizontal')  # adjust the size
	cb1.ax.set_title('R', fontsize=12)
	cb1.ax.tick_params(labelsize=12)   # change the colorbar fontsize

	cax = fig.add_axes([0.46, 0.03, 0.5, 0.02])
	cb2 = fig.colorbar(im2, cax, orientation='horizontal')  # adjust the size
	cb2.ax.tick_params(labelsize=12)   # change the colorbar fontsize
	loc = np.arange(-0.2, 0.201, 0.04)
	labels = ['%.02f' %ll for ll in loc]
	labels[-1] = '  >0.20'
	labels[0] = '<-0.20'
	cb2.set_ticks(loc)
	cb2.ax.set_xticklabels(labels, fontsize=12)
	cb2.ax.set_title('$\Delta$R', fontsize=12)

	fig.tight_layout()
	fig.subplots_adjust(left=0.07, bottom=0.08, right=0.98, wspace=0.0, hspace=0)
	plt.savefig('%s/Fig6_corr_diff_Rmethods7_Cscales4.pdf' %(figdir), dpi=300)
	plt.savefig('%s/Fig6_corr_diff_Rmethods7_Cscales4.png' %(figdir), dpi=300)
	# plt.show()
	return

# Map_correlation_difference_final_methods()
# exit()

def prepare_methods_4scales_corr_for_boxplot():
	"Make a boxplots for 5 methods"

	r1, r3, r6, r12 = [], [], [], []
	list1 = [r1, r3, r6, r12]
	dr_1, dr_3, dr_6, dr_12 = [], [], [], []
	list2 = [dr_1, dr_3, dr_6, dr_12]
	rdr_1, rdr_3, rdr_6, rdr_12 = [], [], [], []
	list3 = [rdr_1, rdr_3, rdr_6, rdr_12]

	for step, r, dr, rdr in zip([1, 3, 6, 12], list1, list2, list3):
		dataplot = load('%s/array_corr_SM_SPEI%02d' % (workspace, step))
		corr0 = dataplot[0,:]
		corr1 = dataplot[1,:]
		# difference of ow - control experiment
		ow_corr_increment = corr1 - corr0
		for i, method, dir, methodname in zip(range(6), methods, dirs, methodnames):
			corr_value = dataplot[i, :]
			corr_diff = corr_value - corr1
			# corr_relative = corr_diff / ow_corr_increment
			r.append(corr_value[~isnan(corr_value)])
			dr.append(corr_diff[~isnan(corr_diff)])
			rdr.append(ow_corr_increment[~isnan(ow_corr_increment)])

		with open('%s/corr_value_SM_SPEI_6methods_SPEI%02d.pkl' %(workspace, step), "wb") as fp:
			pickle.dump(r, fp)
		with open('%s/corr_diff_SM_SPEI_6methods_SPEI%02d.pkl' %(workspace, step), "wb") as fp:
			pickle.dump(dr, fp)
		with open('%s/corr_ow-control_SM_SPEI_6methods_SPEI%02d.pkl' %(workspace, step), "wb") as fp:
			pickle.dump(rdr, fp)

	return

# prepare_methods_4scales_corr_for_boxplot()
# exit()


igbp_colors = {'WB': ('Water', 'Lightblue', 0),
		'ENF': ('Evergreen Needleleaf Forests', 'DarkGreen', 1),
		'EBF': ('Evergreen Broadleaf Forests', 'ForestGreen', 2),
		'DNF': ('Deciduous Needleleaf Forests', 'GreenYellow', 3),
		'DBF': ('Deciduous Broadleaf Forests', 'DarkOliveGreen', 4),
		'MF': ('Mixed Forests', 'Olive', 5),
		'CSH': ('Closed Shrublands', 'DarkKhaki', 6),
		'OSH': ('Open Shrublands', 'Khaki', 7),
		'WSA': ('Woody Savannas', 'SaddleBrown', 8),
		'SAV': ('Savannas', 'BurlyWood', 9),
		'GRA': ('Grasslands', 'YellowGreen', 10),
		'WET': ('Permanent Wetlands', 'PaleGreen', 11),
		'CRO': ('Croplands', 'Gold', 12),
		'URB': ('Urban and Built-Up Lands', 'Grey', 13),
		'MOS': ('Cropland/vegetation Mosaic', 'DarkGoldenrod', 14),
		'SNO': ('Snow and Ice', 'White', 15),
		'BSV': ('Barren or Sparsely Vegetated', 'Wheat', 16)}


def compare_methods_corr_boxplot_aridity_landcover(step):
	"""This is for Figure 7: evaluate the methods across land cover and aridity
	"""
	# Global
	positions_lc = arange(0, 11) # land cover
	positions_cl = arange(0, 3) # aridity
	width = 0.07  # Adjust
	# four methods -> 3 methods
	wspace = (1 - width * 4) / (4 + 1)
	widths = arange(wspace + width / 2, 1.0, wspace + width)
	if step == 1:
		caps = ['a', 'b']; ymin, ymax = -0.1, 0.15
	else:
		caps = ['c', 'd']; ymin, ymax = -0.3, 0.3

	cltypes = ['CONUS', 'Arid', 'Humid']
	lctypes = ['EBF', 'DBF', 'MF', 'ENF', 'WSA', 'SAV', 'CSH', 'OSH', 'GRA', 'CRO', 'MOS']

	lcdata = Dataset('%s/nldas_land_cover_climatology_0.25deg.nc' % (para_dir)).variables[lc['varname']][0, :100, :]
	aridity = read_ETdata_usmask(energy_water_limited_threshold())
	corrfile = '%s/corr.spei.%s.ow.mon%02d.1981-2017.nc' % (dir19, SM, step)
	corr1 = Dataset(corrfile).variables['corr'][0, :]

	def find_lc_group(corrdata, lctypes):
		index = (lcdata == igbp_colors[lctypes[0]][2])
		if len(lctypes) > 1:
			for lctype in lctypes[1:]:
				index = index | (lcdata == igbp_colors[lctype][2])

		data_all = corrdata[index]
		data_all = data_all[~isnan(data_all)]
		return data_all

	def find_aridity_group(corrdata):
		data_arid = corrdata[aridity>2]
		data_humid = corrdata[aridity<=2]

		data_all = corrdata[~isnan(corrdata)]
		data_arid = data_arid[~isnan(data_arid)]
		data_humid = data_humid[~isnan(data_humid)]

		return data_all, data_arid, data_humid

	# Plot
	fig = plt.figure(figsize=(16, 5))
	gs = matplotlib.gridspec.GridSpec(1, 2, width_ratios=[3, 1])
	ax1 = plt.subplot(gs[0])
	ax2 = plt.subplot(gs[1])

	im = 0
	bxpstats_all, bxpstats_arid, bxpstats_humid = [], [], []
	for i, (method, dir, methodname) in enumerate(zip(methods, dirs, methodnames)):
		if i in methods_fig7_list:
			corrfile = '%s/corr.spei.%s.%s.mon%02d.1981-2017.nc' % (dir, SM, method, step)
			corr = Dataset(corrfile).variables['corr'][0, :]

			corr_all, corr_arid, corr_humid = find_aridity_group(read_corrdata_usmask(corr - corr1))
			bxpstats_all.append(corr_all)
			bxpstats_arid.append(corr_arid)
			bxpstats_humid.append(corr_humid)

			bxpstats_lc = []
			for lctype in lctypes:
				corr_lc = find_lc_group(read_corrdata_usmask(corr - corr1), [lctype])
				bxpstats_lc.append(corr_lc)

			### ax1
			plot_loc = positions_lc + widths[im]
			bp1 = ax1.violinplot(bxpstats_lc, positions=plot_loc, showmedians=False, showmeans=True, showextrema=False,
									 widths=width * 2.2)
			for ib, box in enumerate(bp1['bodies']):
				box.set(edgecolor=methodcolors[i], linewidth=0.1, facecolor=methodcolors[i], alpha=0.8)

			v = bp1['cmeans']
			v.set_edgecolor('k')
			v.set_linewidth(2)

			# Black dot means mean
			medians = [np.percentile(bxp, 50, axis=0) for bxp in bxpstats_lc]
			ax1.scatter(plot_loc, medians, marker='o', color='k', s=18, zorder=3, lw=2)

			ax1.set_xlim([-0.2, 11 + 0.2])
			ax1.set_xticks(arange(0.5, 11 + 0.01, 1))
			ax1.set_xticklabels(lctypes, fontsize=16)

			im = im + 1

	### ax2
	for group, bxpstats in enumerate([bxpstats_all, bxpstats_arid, bxpstats_humid]):
		plot_loc = positions_cl[group] + widths
		bp2 = ax2.violinplot(bxpstats, positions=plot_loc, showmedians=False, showmeans=True, showextrema=False, widths=width*2)
		for ib, box in enumerate(bp2['bodies']):
			box.set(edgecolor=methodcolors[methods_fig7_list[ib]], linewidth=0.1, facecolor=methodcolors[methods_fig7_list[ib]], alpha=0.8)

		v = bp2['cmeans']
		v.set_edgecolor('k')
		v.set_linewidth(2)

		medians = [np.percentile(bxp, 50, axis=0) for bxp in bxpstats]
		ax2.scatter(plot_loc, medians, marker='o', color='k', s=18, zorder=3, lw=2)

		ax2.set_xlim([-0.2, 3 + 0.2])
		ax2.set_xticks(arange(0.5, 3 + 0.01, 1))
		ax2.set_xticklabels(cltypes, fontsize=16)


	for i, ax in enumerate([ax1, ax2]):
		ax.set_ylim([ymin, ymax])
		if i == 0:
			ax.set_ylabel('$\Delta$R', fontsize=22)

		for side in ax.spines.keys():
			ax.spines[side].set_linewidth(1.5)
		ax.set_axisbelow(True)
		ax.yaxis.grid(which='major', color='silver', linestyle='-', linewidth=1)
		ax.xaxis.set_ticks_position('bottom')
		ax.tick_params(axis='both', labelsize=16)

	# Add lines
	ax1.vlines([5, 9], -1, 1, colors='k')
	ax1.hlines(0, -1, 12, colors='k')
	ax2.hlines(0, -1, 4, colors='k')
	xmin, xmax = ax1.get_xlim()
	ymin, ymax = ax1.get_ylim()
	ax1.text(xmin + 0.03 * (xmax - xmin), ymin + 0.9 * (ymax - ymin), '(%s)' % caps[0], ha='left', fontsize=18)
	ax1.text(xmin + 0.2 * (xmax - xmin), ymin + 0.9 * (ymax - ymin), 'Forests', ha='left', fontsize=16)
	ax1.text(xmin + 0.5 * (xmax - xmin), ymin + 0.9 * (ymax - ymin), 'Shrublands & Grasslands', ha='left', fontsize=16)
	ax1.text(xmin + 0.85 * (xmax - xmin), ymin + 0.9 * (ymax - ymin), 'Croplands', ha='left', fontsize=16)

	xmin, xmax = ax2.get_xlim()
	ymin, ymax = ax2.get_ylim()
	ax2.text(xmin + 0.1 * (xmax - xmin), ymin + 0.9 * (ymax - ymin), '(%s)' % caps[1], ha='left', fontsize=18)

	# draw temporary lines and use them to create a legend
	lines = [plt.plot([0, 0], c=methodcolors[i], lw=10, alpha=0.7) for i in methods_fig7_list]
	lines = list(itertools.chain(*lines))
	# Set up the legend
	plt.figlegend(lines, ['PET-RC-tall', 'PET-LC', 'PET-CH', 'PET-SW'], loc='lower center', bbox_to_anchor=(0.35, 0.1), numpoints=1, frameon=False,
				  ncol=2,
				  markerscale=2,
				  fontsize=16)
	[line.set_visible(False) for line in lines]
	plt.tight_layout(rect=[0, 0, 1, 1])
	plt.suptitle('SPEI-%s'%step, fontsize=20)
	plt.subplots_adjust(bottom=0.09, top=0.9)
	plt.savefig('%s/Fig7_corr_diff_3methods_lc_aridity_SPEI%s.pdf' %(figdir, step), dpi=300)
	plt.savefig('%s/Fig7_corr_diff_3methods_lc_aridity_SPEI%s.png' %(figdir, step), dpi=300)

	# plt.show()

# compare_methods_corr_boxplot_aridity_landcover(1)
# compare_methods_corr_boxplot_aridity_landcover(12)
# exit()

