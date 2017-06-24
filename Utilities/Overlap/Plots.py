#!/usr/bin/python
import matplotlib
matplotlib.use('GTKAgg') 

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist as axisartist
from matplotlib.dates import DayLocator, HourLocator, DateFormatter
import matplotlib.colors as mc
from netCDF4 import num2date, date2num
import warnings

################# Parameters ###############################################
NMAX     = 1E6
ratio    = 1E-3                                                              # Aspect ratio km/min
Debug    = True
############################################################################

### Color map creation
colors1    = plt.cm.jet(np.linspace(0, 1, 256))
startcolor = np.array([1,1,1,1])
endcolor   = colors1[0]
cmap       = mc.LinearSegmentedColormap.from_list('own',[startcolor,endcolor])
colors2    = cmap(np.linspace(0, 1, 72))
colors     = np.vstack((colors2, colors1))
my_cmap    = mc.LinearSegmentedColormap.from_list('colormap', colors)

def resampling(x, data, z, zmax):
  NX = len(x)
  NZ = len(z)
  
  day_0 = x[0]
  units = day_0.strftime("minutes since %Y-%m-%d 00:00:00")
  xm    = date2num(x,units=units)
  
  DX = xm[1]-xm[0]
  DZ = z[1]-z[0]
  DZinv = 1.0/DZ
  if zmax>0:
    NZ0 = int(round(zmax * DZinv))
    if NZ0<NZ: NZ=NZ0
  
  N  = NX*NZ
  if Debug: print "Original data: {} points".format(N)
  wz = 1
  wx = 1
  while N>NMAX:
    if DX*DZinv*ratio>1: 
      wz = 2*wz
      DZinv = 0.5 * DZinv
    else: 
      wx = 2*wx
      DX = 2*DX
    N = 0.5*N

  NZ = NZ - NZ%wz
  NX = NX - NX%wx
  data = data[:NZ,:NX]
  z    = z[:NZ]
  xm   = xm[:NX]
  
  if wz>1:
    if Debug: print "Using vertical rebin with wz={}".format(wz)
    NZ = NZ/wz
    data_wz = np.full((NZ,NX),np.nan)
    z_wz    = np.mean(z.reshape(-1, wz), axis=1)
    for it in range(NX):
      if not np.all(np.isnan(data[:,it])):
        data_wz[:,it] = np.nanmean(data[:,it].reshape(-1, wz), axis=1)
  else:
    data_wz = data
    z_wz    = z

  if wx>1:
    if Debug: print "Using horizontal rebin with wx={}".format(wx)
    NX = NX/wx
    data_wzwx = np.full((NZ,NX),np.nan)
    xm_wx     = np.mean(xm.reshape(-1, wx), axis=1)
    x_wx      = num2date(xm_wx,units=units)
    with warnings.catch_warnings():
      # I expect to see RuntimeWarnings in this block
      warnings.simplefilter("ignore", category=RuntimeWarning)
      for iz in range(NZ):
        data_wzwx[iz,:] = np.nanmean(data_wz[iz,:].reshape(-1, wx), axis=1)
  else:
    data_wzwx = data_wz
    x_wx = x
  return x_wx, data_wzwx, z_wz

def get_figure(automatic):
  fig = plt.figure(figsize=(15,6))
  ax  = fig.add_subplot(axisartist.Subplot(fig, "111"))
  
  if not automatic:
    ax.xaxis.set_major_locator( DayLocator() )
    ax.xaxis.set_minor_locator( HourLocator([12]) )
  ax.xaxis.set_major_formatter( DateFormatter('%a %d - %H:%M') )
  ax.autoscale_view()
  ax.axis["bottom"].minor_ticklabels.set_rotation(40)
  ax.axis["bottom"].major_ticklabels.set_rotation(40)
  ax.axis["bottom"].minor_ticklabels.set_ha('right')
  ax.axis["bottom"].major_ticklabels.set_ha('right')
  ax.axis["bottom"].label.set_pad(80)
  ax.axis["bottom"].toggle(ticklabels=True)
  
  ax.axis[:].major_ticks.set_tick_out(True)
  
  ax.fmt_xdata = DateFormatter('%Y-%m-%d %H:%M')
  fig.autofmt_xdate()
  return fig, ax
  
def show_raw(x,data,z,zmax=0.0,vmax=20):
  x_low, data_low, z_low = resampling(x,data,z,zmax)
  fig, ax = get_figure(True)

  masked = np.ma.array (data_low, mask=np.isnan(data_low))
  
  CS    = ax.pcolormesh(x_low,z_low,masked, cmap=my_cmap, vmin=0, vmax=vmax)
  cbar  = plt.colorbar(CS)
  cbar.set_label(r"Signal strength [$mv \times km^2$]", labelpad=-63)

  plt.xlabel('')
  plt.ylabel('Height AGL [km]')
  #plt.show()
