#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script provides plotting functions for:
1. initial assessments of different parameterizations
2. map of correlation differences
"""
__author__ = "Liqing Peng"
__copyright__ = "Copyright (C) 2024 Liqing Peng"
__license__ = "MIT"
__version__ = "2024.05"

from pylab import *
from netCDF4 import Dataset
import pandas as pd
import matplotlib.cbook as cbook
import Plotting, DataHandeling
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

maskfile = '/server/user/Masks/0.125deg/us_nldas_0.125deg.nc'
mask = Dataset(maskfile).variables['data'][0,: :]
# set the mask missing values to nan. otherwise it will set to zero when hstacking
mask[mask == 0] = np.nan

SM = 'SMsurf'
steps = [1, 3, 6, 12]
dir24 = '%s/run_2024/monthly_0125' % (wet_home)
dir21 = '%s/run_2021/monthly_0125' % (wet_home)
dir20 = '%s/run_2020/monthly_0125' % (wet_home)
dir19 = '%s/run_2019/monthly_0125' % (wet_home)

methodcolors = ['DodgerBlue', 'cornflowerblue', 'goldenrod', 'crimson', 'green', 'LightSeaGreen']

lc = {'filename': 'nldas_land_cover_climatology.nc', 'varname': 'landcover', 'longname':'"MODIS global land cover climatology (0.5km to 0.125deg)', 'unit':''}

cmap = plt.get_cmap('nipy_spectral')
avg_cmap = Plotting.truncate_colormap(cmap, 0.2, 0.9)  # 0.15, 0.85)

captions = [chr(i) for i in range(ord('a'), ord('z') + 1)]

def Map_climatology_initial():
	"Copy from ETmax_new_analysis"
	# conv_lists = ['ow', 'faotall', 'faoshort', 'pt']
	# conv_titles = ['Open Water', 'FAO Tall reference crop', 'FAO Short reference crop', 'Priestley-Taylor']
	# conv_dirs = ['run_2019', 'run_2019', 'run_2019', 'run_2019']

	# Initial experiments
	gaow_lists = ['ow', 'BL.ga_OW.gs_OW.alb_CLM', 'BL.ga_OW.gs_LAI.alb_OW', 'BL.ga_OW.gs_LAI.alb_CLM'] # BL.ga_OW.gs_OW.alb_OW == OW
	gaow_titles = ['OW', 'Ga OW|Gs OW|'+ r'$\alpha$ CLM', 'Ga OW|Gs LAI|'+ r'$\alpha$ OW', 'Ga OW|Gs LAI|'+r'$\alpha$ CLM']
	gaow_dirs = ['run_2019', 'run_2020', 'run_2020', 'run_2020']

	galc_lists = ['pt', 'BL.ga_LC.gs_OW.alb_CLM', 'BL.ga_LC.gs_LAI.alb_OW', 'BL.ga_LC.gs_LAI.alb_CLM']
	galc_titles = ['PT', 'Ga LC|Gs OW|' + r'$\alpha$ CLM', 'Ga LC|Gs LAI|' + r'$\alpha$ OW', 'Ga LC|Gs LAI|' + r'$\alpha$ CLM']
	galc_dirs = ['run_2019', 'run_2020', 'run_2020', 'run_2020']

	gafao_lists = ['faoshort', 'faotall', 'BL.ga_FAO.gs_OW.alb_FAO', 'BL.ga_FAO.gs_LAI.alb_FAO']
	gafao_titles = ['RC-short', 'RC-tall', 'Ga RC|Gs OW|'+ r'$\alpha$ RC', 'Ga RC|Gs LAI|'+ r'$\alpha$ RC']
	gafao_dirs = ['run_2019', 'run_2019', 'run_2020', 'run_2020']

	names = gaow_titles + galc_titles + gafao_titles
	files = gaow_lists + galc_lists + gafao_lists
	dirs = gaow_dirs + galc_dirs + gafao_dirs


	def save_mean():
		dataplot = []
		for i in xrange(12):
			# run_2019
			pe = Dataset('%s/%s/monthly_0125/pet.%s.monthly.1981-2017.nc'% (wet_home, dirs[i], files[i])).variables['pe'][:, :200,:] #* mask[:200, :]
			gs = DataHandeling.mon2growseason(pe, 1981, 2017)
			gs_clm = np.mean(gs, axis=0)
			dataplot.append(gs_clm)
		dataplot = array(dataplot)*mask[:200, :]*30.0
		dataplot.dump('%s/gs_mean_ETm_12methods_initial' %workspace)
		return

	# save_mean()
	# exit()

	dataplot = load('%s/gs_mean_ETm_12methods_initial' %workspace)

	fig = plt.figure(figsize=(15, 7.5))
	for i in xrange(12):
		ax = fig.add_subplot(3, 4, i+1)
		if i == 0 or i == 4 or i == 8: latlabel = True
		else: latlabel = False
		if i >=8: lonlabel = True
		else: lonlabel = False
		ax.m = Plotting.setup_USmap(latlabel, lonlabel)
		im = ax.m.imshow(dataplot[i], vmin=0, vmax=250, cmap=avg_cmap, alpha=0.7)
		# im = ax.m.imshow(season_avg[i][::-1], vmin=0, vmax=250, cmap=avg_cmap, alpha=0.7)
		# if i == 0:
		# 	ax.m.colorbar(im, location='bottom', pad='15%')
		ax.set_title('(%s) %s' % (captions[i], names[i]), fontsize=16)

	cax = fig.add_axes([0.27, 0.04, 0.48, 0.03])
	cb = fig.colorbar(im, cax, orientation='horizontal')  # adjust the size
	# cb.ax.tick_params(labelsize=12)   # change the colorbar fontsize
	loc = np.arange(0, 250.1, 25)
	labels = ['%d' %ll for ll in loc]
	labels[-1] = '>250'
	cb.set_ticks(loc)
	cb.ax.set_xticklabels(labels, fontsize=14)

	fig.tight_layout()
	fig.subplots_adjust(left=0.04, bottom=0.07, wspace=0.0)
	plt.savefig('%s/Fig2_PET_gs_climatology_12panels_initial.pdf' %(figdir), dpi=300) # png for showing
	# plt.show()
	return

# Map_climatology_initial()
# exit()


def correlation_diff_PET_method_summary_plot_initial():
	"Copy from ETmax_new_analysis"
	gaow_lists = ['ow', 'BL.ga_OW.gs_OW.alb_CLM', 'BL.ga_OW.gs_LAI.alb_OW', 'BL.ga_OW.gs_LAI.alb_CLM']  # BL.ga_OW.gs_OW.alb_OW == OW

	gaow_dirs = ['run_2019', 'run_2020', 'run_2020', 'run_2020']

	galc_lists = ['pt', 'BL.ga_LC.gs_OW.alb_CLM', 'BL.ga_LC.gs_LAI.alb_OW', 'BL.ga_LC.gs_LAI.alb_CLM']
	galc_dirs = ['run_2019', 'run_2020', 'run_2020', 'run_2020']

	gafao_lists = ['faoshort', 'faotall', 'BL.ga_FAO.gs_OW.alb_FAO', 'BL.ga_FAO.gs_LAI.alb_FAO']
	gafao_dirs = ['run_2019', 'run_2019', 'run_2020', 'run_2020']

	files = gaow_lists + galc_lists + gafao_lists
	dirs = gaow_dirs + galc_dirs + gafao_dirs

	def save_data():
		corr_avg = np.empty((4, 12, 4)); corr_avg.fill(np.nan)

		for j, scale in enumerate([1, 3, 6, 12]):
			base = Dataset('%s/run_2019/evaluation/spei/growingseason/corr.spei.SMsurf.%s.mon%02d.1981-2017.nc' % (
			wet_home, 'ow', scale)).variables['corr'][0, :100] * mk25

			for i in xrange(len(files)):
				"SPEI different time scales"
				corr = Dataset('%s/%s/evaluation/spei/growingseason/corr.spei.SMsurf.%s.mon%02d.1981-2017.nc' % (
				wet_home, dirs[i], files[i], scale)).variables['corr'][0, :100] * mk25
				diff = corr - base
				corr_avg[0, i, j] = nanmean(corr)
				corr_avg[1, i, j] = nanmedian(corr)
				corr_avg[2, i, j] = np.float(sum(diff > 0.003)) / sum(~diff.mask) * 100
				corr_avg[3, i, j] = np.float(sum(diff < -0.003)) / sum(~diff.mask) * 100

		corr_avg.dump('%s/spei_SMsurf_PET_12methods_gs_corr_4scale' % workspace)

		return

	# save_data()

	def get_diff(corr_avg):
		# One, Seasonal Gs - no Gs
		dif_gsclm_gsoff = []
		dif_gsclm_gsoff.append(corr_avg[2] - corr_avg[0]) # BL.ga_OW.gs_LAI.alb_OW - ow
		dif_gsclm_gsoff.append(corr_avg[3] - corr_avg[1]) # BL.ga_OW.gs_LAI.alb_CLM - BL.ga_OW.gs_OW.alb_CLM
		dif_gsclm_gsoff.append(corr_avg[7] - corr_avg[5]) # BL.ga_LC.gs_LAI.alb_CLM - BL.ga_LC.gs_OW.alb_CLM
		dif_gsclm_gsoff.append(corr_avg[11] - corr_avg[10]) # BL.ga_FAO.gs_LAI.alb_FAO - BL.ga_FAO.gs_OW.alb_FAO

		# Two, Rough Ga - open water Ga
		dif_gargh_gaow = []
		dif_gargh_gaow.append(corr_avg[5] - corr_avg[1]) # BL.ga_LC.gs_OW.alb_CLM - BL.ga_OW.gs_OW.alb_CLM
		dif_gargh_gaow.append(corr_avg[6] - corr_avg[2]) # BL.ga_LC.gs_LAI.alb_OW - BL.ga_OW.gs_LAI.alb_OW
		dif_gargh_gaow.append(corr_avg[7] - corr_avg[3]) # BL.ga_LC.gs_LAI.alb_CLM - BL.ga_OW.gs_LAI.alb_CLM

		# Three, Seasonal albedo - constant albedo
		dif_albclm_albc = []
		dif_albclm_albc.append(corr_avg[1] - corr_avg[0]) # BL.ga_OW.gs_OW.alb_CLM - ow
		dif_albclm_albc.append(corr_avg[3] - corr_avg[2]) # BL.ga_OW.gs_LAI.alb_CLM - BL.ga_OW.gs_LAI.alb_OW
		dif_albclm_albc.append(corr_avg[7] - corr_avg[6]) # BL.ga_LC.gs_LAI.alb_CLM - BL.ga_LC.gs_LAI.alb_OW

		# Last, Consistent - inconsistent
		dif_cons_incons = []
		# ow is a consistent method, but should not be considered here
		# dif_cons_incons.append(corr_avg[0] - corr_avg[6])
		# dif_cons_incons.append(corr_avg[0] - corr_avg[7])
		# dif_cons_incons.append(corr_avg[0] - corr_avg[8])
		# dif_cons_incons.append(corr_avg[0] - corr_avg[11])
		# dif_cons_incons.append(corr_avg[0] - corr_avg[17])
		# galc
		dif_cons_incons.append(corr_avg[7] - corr_avg[5])
		dif_cons_incons.append(corr_avg[7] - corr_avg[6])
		# garc
		dif_cons_incons.append(corr_avg[11] - corr_avg[10])
		dif_cons_incons.append(corr_avg[8] - corr_avg[10])
		dif_cons_incons.append(corr_avg[9] - corr_avg[10])

		return dif_gargh_gaow, dif_gsclm_gsoff, dif_albclm_albc, dif_cons_incons


	### Draw figure: bar plot (remove the scatters)
	toplabels = ['Surface roughness', 'Surface conductance', 'Albedo', 'Consistency']
	botlabels = ['Rough -\nOpen water', 'Seasonal -\nInfinite', 'Seasonal -\nConstant',
				 'Consistent surface -\nInconsistent surface']

	corr_avg = load('%s/spei_SMsurf_PET_12methods_gs_corr_4scale' % workspace)[1, :,
			   :]  # get the mean correlation 0, or median correlation 1

	dataplot_1 = get_diff(corr_avg[:, 0])
	dataplot_3 = get_diff(corr_avg[:, 1])
	dataplot_6 = get_diff(corr_avg[:, 2])
	dataplot_12 = get_diff(corr_avg[:, 3])


	fig, ax = plt.subplots(figsize=(10, 4), sharey=True)

	def barplot_diff(ax, dataplot, position, boxcolor, label):
		datamn = [mean(dp) for dp in dataplot]
		datasd = [std(dp) for dp in dataplot]

		bp = ax.bar(position, datamn, yerr=datasd, width=0.15, color='none', linewidth=2.5, edgecolor=boxcolor,
					ecolor='k', capsize=5, label=label)
		ax.set_xticklabels(botlabels)
		plt.hold(True)

		return bp

	boxWidth = 0.15
	r1 = np.arange(0.16, 4, 1)
	r2 = [x + boxWidth + 0.02 for x in r1]
	r3 = [x + boxWidth + 0.02 for x in r2]
	r4 = [x + boxWidth + 0.02 for x in r3]

	bp1 = barplot_diff(ax, dataplot_1, r1, 'salmon', 'SPEI-1')
	bp2 = barplot_diff(ax, dataplot_3, r2, 'gold', 'SPEI-3')
	bp3 = barplot_diff(ax, dataplot_6, r3, 'cornflowerblue', 'SPEI-6')
	bp4 = barplot_diff(ax, dataplot_12, r4, 'limegreen', 'SPEI-12')

	ax.tick_params(axis='x', labelsize=14)
	ax.tick_params(axis='y', labelsize=15)
	ax.legend(loc=0, ncol=2, frameon=False)

	ax.set_xlim([0, 4])
	ax.set_ylim([-0.03, 0.03])

	ax.set_xticks(np.arange(0.5, 4.1, 1))

	plot_loc = [0.4, 1.4, 2.6, 3.6]
	for k in xrange(4):
		ax.text(plot_loc[k], 0.032, toplabels[k], ha='center', fontsize=14)
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
	ax.set_axisbelow(True)
	ax.axhline(y=0, color='k', linestyle='--', linewidth=1)
	ax.yaxis.grid(which='major', color='silver', linestyle='--', linewidth=0.5)
	ax.set_ylabel('$\Delta R$', labelpad=2, fontsize=20)

	fig.tight_layout()
	fig.subplots_adjust(left=0.12, top=0.9, right=0.98)

	# plt.show()
	plt.savefig('%s/Fig3_corr_median_summary_4features_barplot.pdf'%figdir, dpi=300)
	plt.savefig('%s/Fig3_corr_median_summary_4features_barplot.png'%figdir, dpi=300)

	return

# correlation_diff_PET_method_summary_plot_initial()
# exit()


