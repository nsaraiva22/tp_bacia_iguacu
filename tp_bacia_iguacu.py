# -*- coding: utf-8 -*-
"""
Created on Sun Dec 18 16:51:11 2022
@author: N.P.Saraiva

Este script lê um arquivo de precipitação diária prevista pelo modelo GEFS e
gera a precipitação média prevista (15 dias) para a sub-bacia do Iguaçu, a
partir da geometria da bacia dipónivel pelo arquivo shapefile

INSTRUÇÕES:
-> linha 28: inserir o caminho de onde estão os arquivos a serem lidos
"""

import os
import warnings
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy, cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader

warnings.filterwarnings('ignore')

path = 'D:/tempo_ok/arq_questao_3/'
os.chdir(path)

# abrindo arquivo netcdf
precip = xr.open_dataset('precip-total_GEFSav_glob_20221212T00.nc')

# dimensões do dado "tp" (total preciptation)
dimen = precip.tp.dims

# salvando os dados das dimensões como variáveis
time   = precip.tp[dimen[0]].data
lat    = precip.tp[dimen[1]].data
lon    = precip.tp[dimen[2]].data
t_pre  = precip.tp.data

# encontrando a data final e inicial do arquivo netcdf
data_i = np.datetime_as_string(time[0], unit='D').split('-')
data_i = str(data_i[2] + '/' + data_i[1] + '/' + data_i[0])
data_f = np.datetime_as_string(time[-1], unit='D').split('-')
data_f = str(data_f[2] + '/' + data_f[1] + '/' + data_f[0])

# calculando a média no tempo com a estrutura do xarray
t_pre_mean = precip.tp.mean('time')

# plot da temperatura média no tempo (apenas para visualização)
precip.tp.mean('time').plot()

###############################################################################
# plot da precipitação total média para a bacia do iguaçu
plt.rcParams['font.family'] = 'Times New Roman'

fig = plt.figure(figsize=(7,6))
fig.suptitle('Precipitação Média Prevista (15 dias)', size=16)

# projeção equidistante
ax = plt.axes(projection=ccrs.PlateCarree())

# Adiciona linha de costa, terra, borda
ax.coastlines(resolution='10m', color='gray')
ax.add_feature(cartopy.feature.LAND, facecolor='silver')
ax.add_feature(cartopy.feature.BORDERS, edgecolor='black')

# Define a área do plote
limites = [-51, -54, -24, -27]
# [min. lon, max. lon, min. lat, max. lat]
ax.set_extent(limites, crs=ccrs.PlateCarree())

# adicionando o dado
data_p = ax.imshow(t_pre_mean, origin='lower', extent=limites, vmin=0,\
    vmax=10, cmap='BuPu')

fig.colorbar(data_p, extend='both').set_label('$\mathregular{kg.m^{-2}}$')

# add um shapefile
shapefile_provinces = list(shpreader.Reader(path+'baixo_iguacu.shp')\
    .geometries())
ax.add_geometries(shapefile_provinces, ccrs.PlateCarree(), edgecolor='gray',\
    facecolor='none', linewidth=0.8)
    
plt.title('%(i)s a %(f)s' %{'i':data_i, 'f':data_f}, fontweight='bold',\
    fontsize=12, loc='left')
plt.title('GEFS', fontweight='bold', fontsize=9, loc='right')

gl = ax.gridlines(linestyle='--', xlocs=np.arange(-54, -50.5, 0.5),\
    ylocs=np.arange(-27, -23.5, 0.5), draw_labels=True)
gl.top_labels = False
gl.right_labels = False

plt.savefig(path+'precip_media.png', bbox_inches='tight', dpi=100)

fig.show() 
