from matplotlib.patches import Polygon as MplPolygon
from pylab import *
from matplotlib import colors
from mpl_toolkits.basemap import Basemap, cm
import os, DataHandeling
# from netCDF4 import Dataset

# wet_home = '/server/user/Data/PROJECT/WET_US'
# # figdir = wet_home + '/Figure/201906'
# figdir = wet_home + '/Figure/202004'
# mondir = '%s/run/monthly/0p125' % (wet_home)
# workspace = '%s/workspace' %(wet_home)
# para_dir = '/server/user/Data/forcing/NLDAS2/parameter/CONUS'
# lc = {'filename': 'nldas_land_cover_climatology.nc', 'varname': 'landcover', 'longname':'"MODIS global land cover climatology (0.5km to 0.125deg)', 'unit':''}
# maskfile = '/server/user/Masks/0.125deg/us_nldas_0.125deg.nc'
# mask = Dataset(maskfile).variables['data'][0,: :]

# mk025file = '/server/user/Masks/0.25deg/us_0.25deg.nc'
# mk25 = Dataset(mk025file).variables['data'][0,: :]

igbp_keys = ['WB', 'ENF', 'EBF', 'DNF', 'DBF', 'MF', 'CSH', 'OSH', 'WSA', 'SAV', 'GRA', 'WET', 'CRO', 'URB', 'MOS', 'SNO', 'BSV']

igbp_colors = {'WB': ('Water', 'Lightblue'),
		'ENF': ('Evergreen Needleleaf Forests', 'DarkGreen'),
		'EBF': ('Evergreen Broadleaf Forests', 'ForestGreen'),
		'DNF': ('Deciduous Needleleaf Forests', 'GreenYellow'),
		'DBF': ('Deciduous Broadleaf Forests', 'DarkOliveGreen'),
		'MF': ('Mixed Forests', 'Olive'),
		'CSH': ('Closed Shrublands', 'DarkKhaki'),
		'OSH': ('Open Shrublands', 'Khaki'),
		'WSA': ('Woody Savannas', 'SaddleBrown'),
		'SAV': ('Savannas', 'BurlyWood'),
		'GRA': ('Grasslands', 'YellowGreen'),
		'WET': ('Permanent Wetlands', 'PaleGreen'),
		'CRO': ('Croplands', 'Gold'),
		'URB': ('Urban and Built-Up Lands', 'Grey'),
		'MOS': ('Cropland/vegetation Mosaic', 'DarkGoldenrod'),
		'SNO': ('Snow and Ice', 'White'),
		'BSV': ('Barren or Sparsely Vegetated', 'Wheat')}

dims_025 = DataHandeling.Set_dims(0.25, 112, 232, 25.125, -124.875)
dims_nldas = DataHandeling.Set_dims(0.125, 224, 464, 25.0625, -124.9375)

rc('font', family='serif')

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
	new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
	return new_cmap

cmap_spectral = plt.get_cmap('nipy_spectral')
cmap_spectral_sub = truncate_colormap(cmap_spectral, 0.05, 0.95) #0.15, 0.85)

cmap_gist = plt.get_cmap('gist_rainbow')
cmap_gist_sub = truncate_colormap(cmap_spectral, 0.15, 0.85)

cmap_rdbu = plt.get_cmap('RdBu')
cmap_rdbu_sub = truncate_colormap(cmap_rdbu,0.05, 0.95)

cmap_seismic = plt.get_cmap('seismic')
cmap_seismic_sub = truncate_colormap(cmap_seismic,0.05, 0.95)

cmap_bwr = plt.get_cmap('bwr')
cmap_bwr_sub = truncate_colormap(cmap_bwr,0, 1)

def setup_USmap(latlabel=True, lonlabel=True):
	m = Basemap(llcrnrlon=-125, llcrnrlat=25, urcrnrlon=-67, urcrnrlat=50, projection='cyl', resolution='c')
	shp_info = m.readshapefile('/server/user/Masks/Shapefile/US_states/states','states',drawbounds=True, linewidth=0.45,color='gray')

	# draw boundaries
	if latlabel == True:
		m.drawparallels(np.arange(30, 55, 10), labels=[1, 0, 0, 0], dashes=[1,1], linewidth=0.25, color='0.5')  # only left ytick
	else:
		m.drawparallels(np.arange(30, 55, 10), labels=[0, 0, 0, 0], dashes=[1,1], linewidth=0.25, color='0.5')  # only left ytick

	if lonlabel == True:
		m.drawmeridians(np.arange(-120, -40, 20), labels=[0, 0, 0, 1], dashes=[1,1], linewidth=0.25, color='0.5')  # only bottom xtick
	else:
		m.drawmeridians(np.arange(-120, -40, 20), labels=[0, 0, 0, 0], dashes=[1,1], linewidth=0.25, color='0.5')  # only bottom xtick


	# m.fillcontinents(color='0.9', lake_color='lightblue', zorder=0)
	# m.drawcoastlines(color='0.6', linewidth=0.5)
	m.drawcountries(color='0.6', linewidth=0.5)
	# m.drawstates(color='0.5', linewidth=0.2)

	ax = plt.gca()

	for shapedict,state in zip(m.states_info, m.states):
		if shapedict['STATE_NAME'] not in ['Alaska', 'Hawaii']: continue
		poly = MplPolygon(state,facecolor='gray',edgecolor='gray')
		ax.add_patch(poly)

	return m

def setup_CONUSmap(latlabel=True):
	m = Basemap(llcrnrlon=-125, llcrnrlat=25, urcrnrlon=-67, urcrnrlat=50, projection='cyl', resolution='c')
	# draw boundaries
	if latlabel == True:
		m.drawparallels(np.arange(30, 55, 10), labels=[1, 0, 0, 0], dashes=[1,1], linewidth=0.25, color='0.5') #, fontsize=15)  # only left ytick
	else:
		m.drawparallels(np.arange(30, 55, 10), labels=[0, 0, 0, 0], dashes=[1,1], linewidth=0.25, color='0.5') #, fontsize=15)  # only left ytick

	m.drawmeridians(np.arange(-120, -40, 20), labels=[0, 0, 0, 1], dashes=[1,1], linewidth=0.25, color='0.5') #, fontsize=15)  # only bottom xtick

	# m.drawmapboundary(fill_color='white', zorder=0)
	m.fillcontinents(color='0.8', lake_color='white', zorder=0)

	m.drawcoastlines(color='0.6', linewidth=0.5)
	m.drawcountries(color='0.3', linewidth=0.7)
	m.drawstates(color='0.5', linewidth=0.2)

	return m

def ContourMap(ax, data, level, unit, res, latlabelflag):

	# ax.m = setup_USmap(latlabelflag)

	# use contourf for the hatch area
	# res = 0.125
	lons = np.arange(-125 + res/2.0, -67., res)
	lats = np.arange(25 + res/2.0, 50., res)
	x, y = np.meshgrid(lons, lats)
	X, Y = ax.m(x, y)

	colors = [plt.cm.seismic(i) for i in [0.2, 0.4, 0.6, 0.8]]

	im = ax.m.contourf(X, Y, data, alpha=0.8, levels=level, colors=colors, hold='on', extend='both')  # edgecolor='0.6', linewidth=0, shading='gouraud',vmin=cblevs[0], vmax=cblevs[-1],
	ax.m.colorbar(im, location='bottom', extend='both', pad='24%', boundaries=level) #ticks=cbticks) # ,this will cut down the boundary in the colorbar, # default size='5%'

	plt.xlabel(unit, fontsize=14, labelpad=18)

	return im