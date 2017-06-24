#!/usr/bin/python

import numpy as np
from netCDF4 import Dataset, num2date
from datetime import datetime, timedelta
from ConfigParser import SafeConfigParser
import Plots

################# Parameters ###############################################
cfgfile       = "../../parameters.cfg"
block         = "Comodoro"
Debug         = True
############################################################################

config = SafeConfigParser()
config.read(cfgfile)

#Read parameters from an external file (string case)
station    = config.get(block, "prefix")
ncpath_out = config.get("Paths", "ncpath_out")
ncfile_out = config.get("Paths", "ncfile_out")
path       = ncpath_out + block + '/'
ncfile_out = "06_2017.nc"

if Debug: print "Opening file: ", path+ncfile_out
### Read Time, Data (mV), Heigth (km)
ds       = Dataset(path+ncfile_out)
bsc532   = ds.variables["bsc532"][:]
bsc1064  = ds.variables["bsc1064"][:]
dep      = ds.variables["dep"][:]
ext_d    = ds.variables["ext_d"][:]
ext_s    = ds.variables["ext_s"][:]
zb       = ds.variables["zb"][:]
zt       = ds.variables["zt"][:]
zinv     = ds.variables["zinv"][:]
z1       = ds.variables["alt1"][:]
z2       = ds.variables["alt2"][:]
times    = ds.variables["time"]
x        = num2date(times[:],units=times.units)
ds.close()

### Number of profiles and vertical levels
NX  = len(x)
NZ1 = len(z1)
NZ2 = len(z2)
  
#Plots.show_beta(x,bsc532.T,z1)
#Plots.show_beta(x,bsc1064.T,z1)
#Plots.show_dep(x,dep.T,z1)
#Plots.show_alpha(x,ext_d.T,z2,zb,zt,zinv,zmax=9.)
Plots.show_alpha(x,ext_s.T,z2,zb,zt,zinv,zmax=9.)


