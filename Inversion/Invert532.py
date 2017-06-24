#!/usr/bin/python

import Calibration
import CloudDetect
import Physics
import Plots
import InOut
import numpy as np
from netCDF4 import Dataset, num2date
from datetime import datetime, timedelta
from ConfigParser import SafeConfigParser
import warnings
import os.path

################# Parameters ###############################################
Debug         = True
############################################################################

def invert(station_block, cfgfile):
  config = SafeConfigParser()
  config.read(cfgfile)

  #Read parameters from an external file (integer case)
  block      = station_block
  wz         = config.getint(block, "wz")
  wsmooth    = config.getint(block, "wsmooth")

  #Read parameters from an external file (float case)
  s1         = config.getfloat(block, "s1")
  depol      = config.getfloat(block, "depol")
  leakrate   = config.getfloat(block, "leakrate")
  mdr        = config.getfloat(block, "mdr")
  dep_s      = config.getfloat(block, "dep_s")
  dep_d      = config.getfloat(block, "dep_d")
  maxint     = config.getfloat(block, "maxint")
  clgrad     = config.getfloat(block, "clgrad")
  clth       = config.getfloat(block, "clth")
  rth1       = config.getfloat(block, "rth1")
  rth4       = config.getfloat(block, "rth4")
  rth6       = config.getfloat(block, "rth6")
  pblth      = config.getfloat(block, "pblth")

  #Read parameters from an external file (string case)
  station    = config.get(block, "prefix")
  ncpath_raw = config.get("Paths", "ncpath_raw")
  ncpath_out = config.get("Paths", "ncpath_out")
  ncfile_raw = config.get("Paths", "ncfile_raw")
  ncfile_out = config.get("Paths", "ncfile_out")
  yrfile_vis = config.get("Paths", "yrfile_vis")
  yrfile_ir  = config.get("Paths", "yrfile_ir")
  power_file = config.get("Paths", "power_file")

  ncpath_out = ncpath_out + block + '/'

  #Read parameters from an external file (boolean case)
  PlotRaw    = config.getboolean("Output", "PlotRaw")
  PlotCal    = config.getboolean("Output", "PlotCal")
  PlotBeta   = config.getboolean("Output", "PlotBeta")
  PlotDep    = config.getboolean("Output", "PlotDep")
  PlotAlpha  = config.getboolean("Output", "PlotAlpha")
  NCDout     = config.getboolean("Output", "NCDout")
  NCDmonth   = config.getboolean("Output", "NCDmonth")

  if os.path.isfile(ncpath_raw+station+ncfile_raw):
    print "Opening file: ", ncpath_raw+station+ncfile_raw
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
  else:
    print "Unable to open file: ", ncpath_raw+station+ncfile_raw
    exit()
    
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
    coeff = np.polyfit(z[-2000:], ch3[-2000:,it], 1)
    pol   = np.poly1d(coeff)
    ch3[:,it] = (ch3[:,it] - pol(z))*z**2
    
  ### Smoothing for IR channel
  if wsmooth>1:
    for it in range(NX):
      ch3[:,it] = np.convolve(ch3[:,it],np.ones(wsmooth)/wsmooth,mode='same')

  if Debug:
    print "Max Value for channel 1: {} mV km2".format(np.nanmax(ch1))
    print "Max Value for channel 2: {} mV km2".format(np.nanmax(ch2))
    print "Max Value for channel 3: {} mV km2".format(np.nanmax(ch3))

  ### Power Correction
  if Debug: print "Performing power corrections..."
  c1, c2, c3 = Calibration.power(x,ncpath_out+power_file)
  for ix in range(NX):
    ch1[:,ix] = ch1[:,ix] * c1[ix]
    ch2[:,ix] = ch2[:,ix] * c2[ix]
    ch3[:,ix] = ch3[:,ix] * c3[ix]
  
  ### Polarization correction (Experimental!)
  depol_estimated = Calibration.depol(ch1, ch2, z)
  if depol_estimated>1.0:
    if Debug: print "Using depol={} instead of depol={}".format(depol_estimated, depol)
    depol = depol_estimated
    
  para   = ch1                                       # Parallel component - Visible channel
  perp   = (ch2 - ch1 * leakrate) / depol            # Perpendicular component - Visible channel
  intvis = para + perp                               # Visible channel (total)
  intir  = ch3                                       # Infrared channel
  with np.errstate(divide='ignore', invalid='ignore'):
    dep  = np.where(para == 0, np.nan, perp / para)  # Volume linear depolarization ratio

  if PlotRaw:
    print "Plotting raw data..."
    Plots.show_raw(x,intvis,z, ncpath_out+"raw_vis.png", zmax=16)
    Plots.show_raw(x,intir, z, ncpath_out+"raw_ir.png", zmax=16)

  ### Overlap Correction
  if Debug: print "Performing overlap corrections..."
  yrvis = Calibration.overlap(ncpath_out + yrfile_vis,z)
  yrir  = Calibration.overlap(ncpath_out + yrfile_ir,z)
  for ix in range(NX):
    intvis[:,ix] = intvis[:,ix] * yrvis
    intir[:,ix]  = intir[:,ix]  * yrir

  ### Calibration
  factor  = Calibration.vis_factor(intvis,z)
  intvis  = intvis * factor       # Attenuated Backscatter Coefficient at 532 nm
  color_r = Calibration.color_ratio(intvis, intir, z, maxint)
  intir   = intir / color_r       # Attenuated Backscatter Coefficient at 1064 nm
  if Debug: 
    print "Calibration factor - Visible channel: {}".format(factor)
    print "Calibration factor - IR channel: {}".format(1.0/color_r)
    print "Vis. channel: min={}, max={}".format(np.nanmin(intvis), np.nanmax(intvis))
    print "IR channel: min={}, max={}".format(np.nanmin(intir), np.nanmax(intir))
    
  ### REBIN function: resizes a vector
  if NZ%wz==0:
    NZ   = NZ/wz
    lvis = np.full((NZ,NX),np.nan)
    lir  = np.full((NZ,NX),np.nan)
    ldep = np.full((NZ,NX),np.nan)
    lz   = np.mean(z.reshape(-1, wz), axis=1)
    DZ = lz[1]-lz[0]
    DZinv = 1.0/DZ
    with warnings.catch_warnings():
      # I expect to see RuntimeWarnings in this block
      warnings.simplefilter("ignore", category=RuntimeWarning)
      for it in range(NX):
        if not np.all(np.isnan(intvis[:,it])) and not np.all(np.isnan(intir[:,it])):
          lvis[:,it] = np.nanmean(intvis[:,it].reshape(-1, wz), axis=1)
          lir[:,it]  = np.nanmean(intir[:,it].reshape(-1, wz), axis=1)
          ldep[:,it] = np.nanmean(dep[:,it].reshape(-1, wz), axis=1)
    if Debug: print "REBIN successful. Vertical resolution: {} m".format(1000*DZ)
  else:
    raise SystemExit("STOP: Rebin not performed: verify your configuration")
    
  ### Molecular profiles
  alpha_m, beta_m = Physics.rayleigh(1000.0*lz, 532.0, height)

  ### Cloud Detection
  zb, zt          = CloudDetect.cloud_height(lvis,lir,lz,clgrad,clth)           # Detect Cloud Base above 240 m
  rf, pbl, invtop = CloudDetect.phenomena(lvis,lir,ldep,lz,zb,rth1,rth4,pblth)

  if PlotCal:
    print "Plotting calibrated signal..."
    Plots.show_cal(x,intvis,z,zb,zt,ncpath_out+"cal_vis.png")
    Plots.show_cal(x,intir,z,zb,zt,ncpath_out+"cal_ir.png")
    #Plots.show_cal(x,intvis,z,pbl,invtop,ncpath_out+"cloud_vis.png")
    #Plots.show_cloud(x,intvis,z,pbl,5*rf,ncpath_out+"cloud_vis.png")
    
  ### Inversion
  ### Vertical indexes
  iz_100 = int(round(0.10 * DZinv))-1
  iz_120 = int(round(0.12 * DZinv))-1
  iz_150 = int(round(0.15 * DZinv))-1
  iz_450 = int(round(0.45 * DZinv))-1
  iz_600 = int(round(0.60 * DZinv))-1
  iz_9km = int(round(9.00 * DZinv))-1
  iz_18km = int(round(18. * DZinv))-1

  ### Optical properties for aerosols
  ext_vis = np.full((NZ,NX), np.nan)       # Extinction coefficient - Visible
  adr     = np.full((NZ,NX), np.nan)       # Particulate depolarization ratio

  ### Use the Fernald's algorithm - High clouds caseext
  for ix in range(NX):
    if rf[ix]>0 or invtop[ix] < 0.21: continue
#    profile_ir  = lir [:,ix]
    profile_vis = 1.0*lvis[:,ix]
    if np.all(np.isnan(profile_vis)): continue
    ### In order to increase the signal-to-noise ratio at inversion height.
    ### This is the boundary condition for the inversion method and 
    ### the whole profiles depend on this boundary condition
    iz_inv = int(round(invtop[ix] * DZinv))-1
    profile_vis[iz_inv-1] = np.nanmean(profile_vis[iz_inv-2:iz_inv+1])
    ### We assume that the molecular contribution 
    ### is dominant at the inversion height
    if zb[ix]<3.0: continue
    else: bsc_ini = beta_m[iz_inv-1] * 1E-4 #Here I made a little modification to the Shimitzu code
    alpha, beta = Physics.fernald(profile_vis, lz, alpha_m, beta_m, s1, invtop[ix], bsc_ini)
    extmin = min(alpha[iz_150:1+max(iz_inv-iz_150, iz_150+1)])
    ###
    for itc in range(1,16):
      if extmin > -1E-5: break
      alpha, beta = Physics.fernald(profile_vis, lz, alpha_m, beta_m, s1, invtop[ix], bsc_ini*2**itc)
      extmin = min(alpha[iz_150:1+max(iz_inv-iz_150, iz_150+1)])
    ###
    sr1  = (beta + beta_m) / beta_m    # Backscatter ratio
    adr1 = (ldep[:,ix] * (sr1 + sr1*mdr - mdr) - mdr) / (sr1 - 1 + sr1*mdr - ldep[:,ix])   # Aerosols depolarization ratio
    ext_vis[iz_100:iz_inv,ix] = alpha[iz_100:iz_inv]
    adr[iz_100:iz_inv,ix]     = adr1[iz_100:iz_inv]
    ### Is this necessary?
    if zb[ix] < 9: 
      iz_zb = int(round(zb[ix] * DZinv))-1
      zmax  = min(zt[ix],9)
      iz_zt = int(round(zmax * DZinv))-1
      ext_vis[iz_zb:iz_zt,ix] = np.nan
      
  ### Verify variability?
  for ix in range(NX):
    profile = ext_vis[iz_120:iz_450,ix]
    if np.all(np.isnan(profile)): 
      continue
    else:
      extstd = np.nanstd(profile)
      if extstd>rth6:
        rf[ix]=6
        ext_vis[:,ix] = np.nan

  time_series = lvis[iz_600,:] / (ext_vis[iz_600,:]/s1 + beta_m[iz_600])
  with np.errstate(invalid='ignore'): mask = ext_vis[iz_600,:]>=0
  if mask.sum()>0:
    ext_int_r = np.median( time_series[mask] )
  else: 
    ext_int_r = 1e11
    
  if Debug: print "Calibration constant: {}".format(ext_int_r)

  ### Use the Fernald's algorithm - Low clouds case
  for ix in range(NX):
    if np.isnan(zb[ix]) or zb[ix]>=3 or invtop[ix] < 0.21: continue
    profile_vis = 1.0*lvis[:,ix]
    if np.all(np.isnan(profile_vis)): continue
    iz_inv      = int(round(invtop[ix] * DZinv))-1
    profile_vis[iz_inv-1] = np.nanmean(profile_vis[iz_inv-2:iz_inv+1])
    ### This boundary condition assumes a transmitance = 1 below of 1km
    ### Above 1km, we follow the same criterion as Shimitzu
    ### This is a very poor estimate
    if invtop[ix]<1.0: 
      bsc_ini = profile_vis[iz_inv] / ext_int_r - beta_m[iz_inv]
    else: 
      bsc_ini = profile_vis[iz_inv] / ext_int_r * 0.8 * np.exp(invtop[ix]/4.5) - beta_m[iz_inv]
    if bsc_ini<0:
      if Debug: print "Changing boundary condition bsc_ini = {}".format(bsc_ini)
      bsc_ini = beta_m[iz_inv-1] * 1E-4
    alpha, beta = Physics.fernald(profile_vis, lz, alpha_m, beta_m, s1, invtop[ix], bsc_ini)
    ###
    sr1  = (beta + beta_m) / beta_m
    adr1 = (ldep[:,ix] * (sr1 + sr1*mdr - mdr) - mdr) / (sr1 - 1 + sr1*mdr - ldep[:,ix])
    ###
    ext_vis[iz_100:iz_inv,ix] = alpha[iz_100:iz_inv]
    adr[iz_100:iz_inv,ix]     = adr1[iz_100:iz_inv]
    ###
    iz_zb = int(round(zb[ix] * DZinv))-1
    zmax  = min(zt[ix],9)
    iz_zt = int(round(zmax * DZinv))-1
    ext_vis[iz_zb:iz_zt,ix] = np.nan

  absc532  = lvis / ext_int_r
  absc1064 = lir  / ext_int_r
  ###
  dr = (adr - dep_s) * (1 + dep_d) / (1 + adr) / (dep_d - dep_s)
  with np.errstate(invalid='ignore'): 
    dr[dr<0]=0.0
    dr[dr>1]=1.0
  dust    = ext_vis * dr
  sphere  = ext_vis * (1.0-dr)
  bsc     = ext_vis / s1

  if Debug: print "Max. Att. Backscatter: {}".format(np.nanmax(absc532))

  if PlotBeta:
    print "Plotting attenuated backscatter coefficients..."
    Plots.show_beta(x,1000.0*absc532,lz,ncpath_out+'absc_vis.png',zmax=18)
    Plots.show_beta(x,1000.0*absc1064,lz,ncpath_out+'absc_ir.png',zmax=18)

  if PlotAlpha:
    print "Plotting extinction coefficients..."
    Plots.show_alpha(x,1000.0*dust,lz,zb,zt,invtop,ncpath_out+'dust.png',zmax=9)
    Plots.show_alpha(x,1000.0*sphere,lz,zb,zt,invtop,ncpath_out+'sphere.png',zmax=9)

  if PlotDep:
    print "Plotting depolarization ratio..."
    with np.errstate(invalid='ignore'):
      ldep[absc532<1E-7]=np.nan
    Plots.show_dep(x,ldep,lz,ncpath_out+'dep.png',zmax=18)
    
  if NCDout:
    print "Creating NetCDF file: {}".format(ncpath_out+ncfile_out)
    NZ1 = iz_18km+1
    NZ2 = iz_9km+1
    InOut.save_ncd(ncpath_out+ncfile_out, station, x, lz, absc532, absc1064, ldep, dust, sphere, zb, zt, invtop, NZ1, NZ2)
    
  if NCDmonth:
    print "Creating monthly NetCDF files"
    NZ1 = iz_18km+1
    NZ2 = iz_9km+1
    InOut.monthly_ncd(ncpath_out, station, x, lz, absc532, absc1064, ldep, dust, sphere, zb, zt, invtop, NZ1, NZ2)
