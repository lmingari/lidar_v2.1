#!/usr/bin/python

from netCDF4   import Dataset, num2date, date2num
from datetime  import datetime, timedelta, time
import numpy as np
import os

################# Parameters ###############################################
Debug    = True
############################################################################

def save_ncd(filename, station, x, lz, absc532, absc1064, ldep, dust, sphere, zb, zt, invtop, NZ1, NZ2):
  if Debug:
    print "**********************"
    print "Creating output file: {}".format(filename)
  ### Dimensions
  ncfile = Dataset(filename,'w',format="NETCDF3_CLASSIC") 
  ncfile.createDimension("time", None)
  ncfile.createDimension("alt1", NZ1)
  ncfile.createDimension("alt2", NZ2)

  ### Coordinate variables
  time  = ncfile.createVariable("time","f4",("time",))
  alt1  = ncfile.createVariable("alt1","f4",("alt1",))
  alt2  = ncfile.createVariable("alt2","f4",("alt2",))

  ### Coordinate variable attributes
  alt1.units = "km"
  alt1.description = "altitude"
  alt1[:] = lz[:NZ1]

  alt2.units = "km"
  alt2.description = "altitude"
  alt2[:] = lz[:NZ2]
  
  day_0 = x[0]
  time.units = day_0.strftime("minutes since %Y-%m-%d 00:00:00")
  time.description = "time after 0000UTC"
  time[:] = date2num(x,units=time.units)

  ### Variables
  bsc1 = ncfile.createVariable("bsc532","f4",("time","alt1",))
  bsc1.units = "/sr /km"
  bsc1.description = "Attenuated Backscatter coefficient (532 nm)"
  bsc1[:,:] = 1000.0*(absc532[:NZ1,:]).T

  bsc2 = ncfile.createVariable("bsc1064","f4",("time","alt1",))
  bsc2.units = "/sr /km"
  bsc2.description = "Attenuated Backscatter coefficient (1064 nm)"
  bsc2[:,:] = 1000.0*(absc1064[:NZ1,:]).T

  dep1 = ncfile.createVariable("dep","f4",("time","alt1",))
  dep1.units = ""
  dep1.description = "Volume Depolarization ratio"
  dep1[:,:] = (ldep[:NZ1,:]).T

  ext1 = ncfile.createVariable("ext_d","f4",("time","alt2",))
  ext1.units = "/km"
  ext1.description = "Extinction coefficient - Non spherical particles"
  ext1[:,:] = 1000.0*(dust[:NZ2,:]).T

  ext2 = ncfile.createVariable("ext_s","f4",("time","alt2",))
  ext2.units = "/km"
  ext2.description = "Extinction coefficient - Spherical particles"
  ext2[:,:] = 1000.0*(sphere[:NZ2,:]).T

  zb1 = ncfile.createVariable("zb","f4",("time",))
  zb1.units = "km"
  zb1.description = "Cloud Bottom"
  zb1[:] = zb

  zt1 = ncfile.createVariable("zt","f4",("time",))
  zt1.units = "km"
  zt1.description = "Cloud Top"
  zt1[:] = zt

  zinv1 = ncfile.createVariable("zinv","f4",("time",))
  zinv1.units = "km"
  zinv1.description = "Inversion height"
  zinv1[:] = invtop

### Global Attributes
  ncfile.TITLE   = "LIDAR products"
  ncfile.YEAR    = day_0.year
  ncfile.MONTH   = day_0.month
  ncfile.DAY     = day_0.day
  ncfile.STATION = station

  ncfile.close()
  if Debug: print "Done!"
    
def update_ncd(filename, x, absc532, absc1064, ldep, dust, sphere, zb, zt, invtop, NZ1, NZ2):
  if Debug:
    print "**********************"
    print "Updating file: {}".format(filename)
  NT = len(x)
  ncfile = Dataset(filename,'a',format="NETCDF3_CLASSIC") 
  #Dimensions
  NT_file  = len(ncfile.dimensions["time"])
  NZ1_file = len(ncfile.dimensions["alt1"])
  NZ2_file = len(ncfile.dimensions["alt2"])
  #Variables
  t     = ncfile.variables["time"]
  t_end = num2date(t[-1], units = t.units)
  v1    = ncfile.variables["bsc532"]
  v2    = ncfile.variables["bsc1064"]
  v3    = ncfile.variables["dep"]
  v4    = ncfile.variables["ext_d"]
  v5    = ncfile.variables["ext_s"]
  v6    = ncfile.variables["zb"]
  v7    = ncfile.variables["zt"]
  v8    = ncfile.variables["zinv"]
  
  print NZ1_file, NZ2_file

  #Check file consistence  
  if x[0]==t_end and NZ1==NZ1_file and NZ2==NZ2_file:
    for it in range(NT):
      t[NT_file+it-1]    = date2num(x[it], units = t.units)
      v1[NT_file+it-1,:] = absc532[:NZ1,it]
      v2[NT_file+it-1,:] = absc1064[:NZ1,it]
      v3[NT_file+it-1,:] = ldep[:NZ1,it]
      v4[NT_file+it-1,:] = dust[:NZ2,it]
      v5[NT_file+it-1,:] = sphere[:NZ2,it]
      v6[NT_file+it-1]   = zb[it]
      v7[NT_file+it-1]   = zt[it]
      v8[NT_file+it-1]   = invtop[it]
  else:
    print "Error: File is not consistent with input data"

  ncfile.close()
  if Debug: print "Done!"
    
def monthly_ncd(path, station, x, lz, absc532, absc1064, ldep, dust, sphere, zb, zt, invtop, NZ1, NZ2):
  if Debug: print "**** Function monthly_ncd ****"
  NX = len(x)
  i1 = NX
  fname = ""
  for ix in range(NX):
    fname_new = x[ix].strftime("%m_%Y") + ".nc"
    if fname != fname_new:
      if fname and ix>i1+1:
        #Save data to fname
        if os.path.isfile(path+fname):
          #Append data
          update_ncd(path+fname, 
                     x[i1:ix], 
                     absc532[:,i1:ix], 
                     absc1064[:,i1:ix], 
                     ldep[:,i1:ix], 
                     dust[:,i1:ix], 
                     sphere[:,i1:ix], 
                     zb[i1:ix], 
                     zt[i1:ix], 
                     invtop[i1:ix], 
                     NZ1, NZ2)
        else:
          #Create a new file
          save_ncd(path+fname, station, 
                   x[i1:ix], lz, 
                   absc532[:,i1:ix], 
                   absc1064[:,i1:ix], 
                   ldep[:,i1:ix], 
                   dust[:,i1:ix], 
                   sphere[:,i1:ix], 
                   zb[i1:ix], 
                   zt[i1:ix], 
                   invtop[i1:ix], 
                   NZ1, NZ2)
            
      #Inquiry information on fname_new
      if os.path.isfile(path+fname_new):
        if Debug: print "Found file: ", fname_new
        ncfile = Dataset(path+fname_new,'r',format="NETCDF3_CLASSIC") 
        t      = ncfile.variables["time"]
        x_in   = num2date(t[-1], units = t.units)
        ncfile.close()
      else:
        if Debug: print "Not found file: ", fname_new
        x_in = x[ix]
      i1 = ix
      fname = fname_new
    else:
      if x[i1]<x_in: i1=ix

  #Special case, for the last time:
  if os.path.isfile(path+fname):
    #Append data
    update_ncd(path+fname, 
               x[i1:], 
               absc532[:,i1:], 
               absc1064[:,i1:], 
               ldep[:,i1:], 
               dust[:,i1:], 
               sphere[:,i1:], 
               zb[i1:], 
               zt[i1:], 
               invtop[i1:],
               NZ1, NZ2)
  else:
    #Create a new file
    save_ncd(path+fname, station, 
             x[i1:], lz, 
             absc532[:,i1:], 
             absc1064[:,i1:], 
             ldep[:,i1:], 
             dust[:,i1:], 
             sphere[:,i1:], 
             zb[i1:], 
             zt[i1:], 
             invtop[i1:], 
             NZ1, NZ2)
