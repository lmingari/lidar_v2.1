#!/usr/bin/python

from LoadLicel import LoadLicel
from netCDF4   import Dataset, num2date, date2num
from datetime  import datetime, timedelta, time
import numpy as np
import os

################# Parameters ###############################################
dt       = timedelta(minutes=1)                                            # Delta time (1 min). NO CHANGE
#FileSize = 82492                                                           # Filesize in bytes
Debug    = True
############################################################################

def findtimes(t1,t2,sampling,path,prefix,ncfile,FileSize):
  print "Using File size: {}".format(FileSize)
  if os.path.isfile(ncfile):
    if Debug: print "Found file: ", ncfile
    ncfile = Dataset(ncfile,'r',format="NETCDF3_CLASSIC") 
    t      = ncfile.variables["time"]
    t1_out = num2date(t[-1], units = t.units)
    ncfile.close()
  else:
    if Debug: print "Not found file: ", ncfile
    t = t1.date()
    SearchFile = True
    while t<=t2.date() and SearchFile:
      folder = t.strftime("%Y%m%d/")
      if os.path.exists(path + folder): 
        filename = prefix + t.strftime("%y") + format(t.month,'X') + t.strftime("%d")
        filelist_folder = os.listdir(path + folder)
        filelist_folder = filter(lambda x: filename in x , filelist_folder)
        filelist_folder = filter(lambda x: os.stat(path + folder + x).st_size == FileSize, filelist_folder)
        if filelist_folder:
          if t1.date() == t:
            tt       = t1
          else:
            tt       = datetime.combine(t, time())
            dt_shift = timedelta( minutes = int((tt-t1).total_seconds()//60.0) % sampling )
            tt       = tt - dt_shift
          ChangeFolder = False
          while SearchFile and not ChangeFolder:
            for it in range(sampling):
              t0 = tt + it * dt
              if t0.date()>t: 
                ChangeFolder = True
                break
              else:
                filename = prefix + t0.strftime("%y") + format(t0.month,'X') + t0.strftime("%d%H.%M")
              if any([filename in item for item in filelist_folder]):
                SearchFile = False
                t1_out = tt
                break
            tt = tt + dt*sampling
      t = t + timedelta(days=1)
    if SearchFile: t1_out = t2

  t = t2.date()
  SearchFile = True
  while t>=t1_out.date() and SearchFile:
    folder = t.strftime("%Y%m%d/")
    if os.path.exists(path + folder): 
      filename = prefix + t.strftime("%y") + format(t.month,'X') + t.strftime("%d")
      filelist_folder = os.listdir(path + folder)
      filelist_folder = filter(lambda x: filename in x , filelist_folder)
      filelist_folder = filter(lambda x: os.stat(path + folder + x).st_size == FileSize, filelist_folder)
      if filelist_folder:
        if t2.date() == t: tt = t2
        else: tt = datetime.combine(t+timedelta(days=1), time())
        dt_shift = timedelta( minutes = int((tt-t1).total_seconds()//60.0) % sampling )
        tt       = tt - dt_shift
        ChangeFolder = False
        while SearchFile and not ChangeFolder:
          for it in reversed(range(sampling)):
            t0 = tt + it * dt
            if t0.date()<t: 
              ChangeFolder = True
              break
            else:
              filename = prefix + t0.strftime("%y") + format(t0.month,'X') + t0.strftime("%d%H.%M")
            if any([filename in item for item in filelist_folder]):
              SearchFile = False
              t2_out = tt
              break
          tt = tt - dt*sampling
    t = t - timedelta(days=1)
  if SearchFile: t2_out = t1_out
            
  return t1_out, t2_out
  
def get_data(t1,t2,sampling,path,prefix,FileSize):
  x_time=[]
  t=t1
  while t<=t2:
    x_time.append(t)
    t=t+sampling*dt
  NT=len(x_time)
  x_time    = np.array(x_time)
  n_data    = np.zeros(NT)
  FileFound = False
  t = t1.date()
  while t<=t2.date():
    folder = t.strftime("%Y%m%d/")
    if os.path.exists(path + folder):
      filename = prefix + t.strftime("%y") + format(t.month,'X') + t.strftime("%d")
      filelist_folder = os.listdir(path + folder)
      filelist_folder = filter(lambda x: filename in x and len(x)==15, filelist_folder)
      filelist_folder = filter(lambda x: os.stat(path + folder + x).st_size == FileSize, filelist_folder)
      if filelist_folder:
        print "Reading folder: ", folder
        for file_item in filelist_folder:
          file_date = datetime.strptime(folder+file_item[:-4], "%Y%m%d/"+filename+"%H.%M")
          i_time = int((file_date-t1).total_seconds()//(60.0*sampling))
          if 0<=i_time<NT:
            n_data[i_time] = n_data[i_time] + 1
            # Get Licel data
            x = LoadLicel(path + folder + file_item)
            if not FileFound:
              height = x.globalparameters.HeightASL
              z      = x.channel[0].Range * 0.001      # Height in kilometers
              NZ     = len(z)
              y1     = np.zeros((NT,NZ))
              y2     = np.zeros((NT,NZ))
              y3     = np.zeros((NT,NZ))
            for ichannel in range(x.globalparameters.Channels):
              if not x.channel[ichannel].isPhotonCounting:
                if x.channel[ichannel].Wavelength == 532 and not x.channel[ichannel].isPolarized:
                  y1[i_time,:] = y1[i_time,:] + x.channel[ichannel].Signal
                  if Debug and not FileFound:
                    print "************ ch1 ************"
                    print "Using channel: {} nm".format(x.channel[ichannel].Wavelength)    # Show data information
                    for key, item in x.globalparameters.__dict__.items():
                        print key + ": ", item
                elif x.channel[ichannel].Wavelength == 532 and x.channel[ichannel].isPolarized:
                  y2[i_time,:] = y2[i_time,:] + x.channel[ichannel].Signal
                  if Debug and not FileFound:
                    print "************ ch2 ************"
                    print "Using channel: {} nm".format(x.channel[ichannel].Wavelength)    # Show data information
                    for key, item in x.globalparameters.__dict__.items():
                        print key + ": ", item
                elif x.channel[ichannel].Wavelength == 1064 and not x.channel[ichannel].isPolarized:
                  y3[i_time,:] = y3[i_time,:] + x.channel[ichannel].Signal
                  if Debug and not FileFound:
                    print "************ ch3 ************"
                    print "Using channel: {} nm".format(x.channel[ichannel].Wavelength)    # Show data information
                    for key, item in x.globalparameters.__dict__.items():
                        print key + ": ", item
            if not FileFound: FileFound = True
    t = t + timedelta(days=1)
  for it in range(NT):
    if n_data[it]==0:
      y1[it,:] = np.nan
      y2[it,:] = np.nan
      y3[it,:] = np.nan        
    else:
      #print x_time[it], n_data[it], np.max(y1[it,:])
      y1[it,:] = y1[it,:] / n_data[it]
      y2[it,:] = y2[it,:] / n_data[it]
      y3[it,:] = y3[it,:] / n_data[it]
  return x_time, y1, y2, y3, z, height
  
def createncd(x, data1, data2, data3, z, height, ncfilename):
  NT = len(x)
  NZ = len(z)
  print "**********************"
  print "Creating file: {}".format(ncfilename)
  #Dimensions
  ncfile = Dataset(ncfilename,'w',format="NETCDF3_CLASSIC") 
  ncfile.createDimension("time", None)
  ncfile.createDimension("alt", NZ)

  #Coordinate variables
  time = ncfile.createVariable("time","f4",("time",))
  altitude = ncfile.createVariable("alt","f4",("alt",))

  #Coordinate variable attributes
  altitude.units = "km"
  altitude.description = "altitude"
  altitude[:] = z

  day_0 = datetime(x[0].year,x[0].month,x[0].day)
  time.units = day_0.strftime("minutes since %Y-%m-%d 00:00:00")
  time.description = "time after 0000UTC"
  time[:] = [item.total_seconds()/60.0 for item in x-day_0]

  #Variables
  ch1 = ncfile.createVariable("ch1","f4",("time","alt",))
  ch1.units = "mV"
  ch1.description = "532-nm channel - component: o"
  ch1[:,:] = data1

  ch2 = ncfile.createVariable("ch2","f4",("time","alt",))
  ch2.units = "mV"
  ch2.description = "532-nm channel - component: p"
  ch2[:,:] = data2

  ch3 = ncfile.createVariable("ch3","f4",("time","alt",))
  ch3.units = "mV"
  ch3.description = "1064-nm channel"
  ch3[:,:] = data3

#Global Attributes
  ncfile.TITLE  = "LIDAR data from Licel"
  ncfile.YEAR   = day_0.year
  ncfile.MONTH  = day_0.month
  ncfile.DAY    = day_0.day
  ncfile.HEIGHT = height

  ncfile.close()
  print "Done!"
  print "**********************"
  
def updatencd(x, data1, data2, data3, z, ncfilename):
  NT = len(x)
  NZ = len(z)
  print "**********************"
  print "Updating file: {}".format(ncfilename)
  ncfile = Dataset(ncfilename,'a',format="NETCDF3_CLASSIC") 
  #Dimensions
  NT_file = len(ncfile.dimensions["time"])
  NZ_file = len(ncfile.dimensions["alt"])
  #Variables
  t       = ncfile.variables["time"]
  t_end   = num2date(t[-1], units = t.units)
  ch1     = ncfile.variables["ch1"]
  ch2     = ncfile.variables["ch2"]
  ch3     = ncfile.variables["ch3"]
  
  #Check file consistence  
  if NZ==NZ_file and x[0]==t_end:
    for it in range(NT):
      t[NT_file+it-1]     = date2num(x[it], units = t.units)
      ch1[NT_file+it-1,:] = data1[it,:]
      ch2[NT_file+it-1,:] = data2[it,:]
      ch3[NT_file+it-1,:] = data3[it,:]
  else:
    print "Error: File is not consistent with input data"

  ncfile.close()
  print "Done!"
  print "**********************"
