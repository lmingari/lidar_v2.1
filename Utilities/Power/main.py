#!/usr/bin/python

import numpy as np
from netCDF4 import Dataset, num2date
from datetime import datetime, timedelta
from ConfigParser import SafeConfigParser
import Plots

################# Parameters ###############################################
cfgfile       = "../../parameters.cfg"
block         = "Cordoba"
Debug         = True
############################################################################

config = SafeConfigParser()
config.read(cfgfile)

#Read parameters from an external file (string case)
station    = config.get(block, "prefix")
ncpath_raw = config.get("Paths", "ncpath_raw")
ncfile_raw = config.get("Paths", "ncfile_raw")

### Read Time, Data (mV), Heigth (km)
ds       = Dataset(ncpath_raw+station+ncfile_raw)
ch1      = ds.variables["ch1"][:]
ch2      = ds.variables["ch2"][:]
ch3      = ds.variables["ch3"][:]
z        = ds.variables["alt"][:]
times    = ds.variables["time"]
day_0    = datetime(year=ds.YEAR, month=ds.MONTH, day=ds.DAY)
height   = ds.HEIGHT
x        = num2date(times[:],units=times.units)
ds.close()

ch1      = ch1.T
ch2      = ch2.T
ch3      = ch3.T

### Number of profiles and vertical levels
NX = len(x)
NZ = len(z)

### Background removed and range-corrected signal intensity
for it in range(NX):
  coeff = np.polyfit(z[-1000:], ch1[-1000:,it], 1)
  pol   = np.poly1d(coeff)
  ch1[:,it] = (ch1[:,it] - pol(z))*z**2
  ###
  coeff = np.polyfit(z[-1000:], ch2[-1000:,it], 1)
  pol   = np.poly1d(coeff)
  ch2[:,it] = (ch2[:,it] - pol(z))*z**2
  ###
  coeff = np.polyfit(z[-1000:], ch3[-1000:,it], 1)
  pol   = np.poly1d(coeff)
  ch3[:,it] = (ch3[:,it] - pol(z))*z**2

if Debug:
  print "Max Value for channel 1: {} mV km2".format(np.nanmax(ch1))
  print "Max Value for channel 2: {} mV km2".format(np.nanmax(ch2))
  print "Max Value for channel 3: {} mV km2".format(np.nanmax(ch3))

def power(x, filename):
  print "**** Function power ****"
  c1     = np.ones_like(x, dtype=np.float)
  c2     = np.ones_like(x, dtype=np.float)
  c3     = np.ones_like(x, dtype=np.float)
  dtype1 = np.dtype([('time', 'object'), ('cvis_par', 'f8'), ('cvis_per', 'f8'), ('cir', 'f8')])
  fmt    = "%Y-%m-%d_%H:%M"
  try:
    data = np.loadtxt(filename, dtype=dtype1, converters={0: lambda x: datetime.strptime(x, fmt)})
    for item in data:
      if item['time']<=x[0]: 
        continue
      elif item['time']>=x[-1]: 
        break
      else:
        if Debug: 
          print "Applying power correction for ", item["time"]
          print "cvis_par={} - cvis_per={} - cir={}".format(item['cvis_par'],item['cvis_per'],item['cir'])
        for ix in range(NX):
          if item['time']<x[ix]:
            c1[ix] = c1[ix] * item['cvis_par']
            c2[ix] = c2[ix] * item['cvis_per']
            c3[ix] = c3[ix] * item['cir']
  except IOError:
    if Debug: print "Power file not found: {}".format(filename)
  return c1, c2, c3
  
#Plots.show_raw(x,ch3,z,zmax=18.,vmax=1.)

c1, c2, c3 = power(x,"power.cfg")

for ix in range(NX):
  ch1[:,ix] = ch1[:,ix] * c1[ix]
  ch2[:,ix] = ch2[:,ix] * c2[ix]
  ch3[:,ix] = ch3[:,ix] * c3[ix]
  
DZ = z[1]-z[0]
DZinv = 1.0/DZ
iz = int(round(0.6 * DZinv))-1
Plots.xy(x,ch2[iz,:])
#Plots.show_raw(x,ch1,z,zmax=18.,vmax=10.)
