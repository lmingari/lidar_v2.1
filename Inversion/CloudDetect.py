#!/usr/bin/python
import numpy as np

################# Parameters ####################
Debug  = True
#################################################

#if not Debug: np.seterr(invalid='ignore')

def cloud_height(vis,ir,z,clgrad,clth):
  ### z in km
  NX = ir.shape[-1]
  NZ = ir.shape[0]
  DZ = z[1]-z[0]
  zb = np.full(NX,np.nan)
  zt = np.full(NX,np.nan)
  ###
  for ix in range(NX):
    profile_ir  = ir [:,ix]
    profile_vis = vis[:,ix]
    diff        = (profile_ir[1:]-profile_ir[:-1])/DZ
    if np.all(np.isnan(profile_ir)): continue
    for iz in range(NZ-2):
#      if z[iz]<0.24: continue
      if z[iz]<0.35: continue
      clgrad1 = clgrad * np.exp(-1 * (z[iz] / 3.0))
      if diff[iz]>=clgrad1/3E-2 and diff[iz+1]>=clgrad1/3E-2:
        cb1 = iz
        ct1 = iz
        while ct1<NZ-1:
          if profile_vis[ct1]<profile_vis[cb1]: break
          ct1 = ct1 + 1
        clth1 = clth * np.exp(-1. * z[iz]**2 / 90.0)
        circ  = (profile_ir [cb1:ct1]>=clth1).sum()
        cvisc = (profile_vis[cb1:ct1]>=clth1).sum()
        ###
        if min(circ, cvisc) >= (ct1 - cb1) * 0.6 and (z[ct1] - z[cb1]) >= (z[cb1] - 9.0) / 5:
          zb[ix] = z[cb1]
          zt[ix] = z[ct1]
          break
  return zb, zt

def phenomena(vis,ir,dep,z,zb,rth1,rth4,pblth):
  NX = ir.shape[-1]
  NZ = ir.shape[0]
  DZ = z[1]-z[0]
  DZinv = 1.0/DZ

  iz_120  = int(round(0.12 * DZinv))-1
  iz_240  = int(round(0.24 * DZinv))-1
  iz_300  = int(round(0.3  * DZinv))-1
  iz_450  = int(round(0.45 * DZinv))-1
  iz_600  = int(round(0.6  * DZinv))-1
  iz_1500 = int(round(1.5  * DZinv))-1
  iz_3000 = int(round(3.0  * DZinv))-1
  iz_9km  = int(round(9.0  * DZinv))-1
  iz_15km = int(round(15.0  * DZinv))-1
  iz_16km = int(round(16.0  * DZinv))-1
  iz_18km = int(round(18.0  * DZinv))-1

  rf     = np.zeros(NX)
  pbl    = np.full(NX,np.nan)
  invtop = np.full(NX,np.nan)
  nlg    = np.full(NX,9.0)
  for ix in range(NX):
    profile_ir  = ir [:,ix]
    profile_vis = vis[:,ix]
    profile_dep = dep[:,ix]
    if np.all(np.isnan(profile_ir)): continue
    diff        = (profile_ir[1:]-profile_ir[:-1])*DZinv
    with np.errstate(divide='ignore', invalid='ignore'):
      cr = np.where(profile_vis == 0, np.nan, profile_ir/profile_vis)
    ###
    std_tail = np.nanstd(profile_vis[iz_15km:iz_16km]) * 0.1
    for iz in range(iz_450,iz_9km): 
      if profile_vis[iz]<std_tail:
        nlg[ix] = z[iz]
        break
    ### Inversion must start > 600m
    zmax       = np.nanmin([zb[ix],9])
    invtop[ix] = np.nanmin([zb[ix] - DZ * 2**(zmax/1.5) - 0.15, nlg[ix], 9.0])
    ### Rain
    zmax       = np.nanmin([zb[ix],3])
    iz_max     = int(round(zmax * DZinv))-1
    with np.errstate(invalid='ignore'):
      count = np.nansum(cr[iz_240:iz_max]>=1.1)
    zmin = 0.24
    if count * DZ >= (zmax-zmin) * rth1: 
      rf[ix] = 1
      continue
    ### Fog?
    if nlg[ix]<=0.6:
      rf[ix] = 2
      continue
    ### Snow
    if np.nanmean(profile_vis[iz_120:iz_240]) >= 5E5 and np.nanmean(profile_dep[iz_120:iz_240]) >= 0.2:
      rf[ix] = 3
      continue
    zmax   = np.nanmax([zb[ix],1.5])
    iz_max = int(round(zmax * DZinv))-1
    if np.nanmin(diff[iz_240:iz_max]) < rth4/3E-2:
      rf[ix] = 4
      continue
    if np.nanmax(profile_vis[iz_240:iz_max]) >= 5E6 and np.nanmax(profile_dep[iz_240:iz_max]) >= 0.3:
      rf[ix] = 5
      continue
    ### PBL
    zmax   = np.nanmin([zb[ix]-0.15,4.5])
    iz_top = int(round(zmax * DZinv))-1
    with np.errstate(invalid='ignore'):
      for iz in range(iz_600,iz_top):
        if np.nansum(profile_vis[iz-iz_300:iz])/np.nansum(profile_vis[iz+1:iz+iz_300+1]) >= pblth and \
           max(profile_vis[iz_300:iz]) < min(profile_vis[iz_300:iz]) * 2.0:
          pbl[ix] = z[iz]
          continue

  return rf, pbl, invtop

if __name__ == "__main__":

  a=np.array([[1,2,3],[40,50,60],[100,200,400]])
  b=np.array([[1,2,3],[40,np.nan,60],[100,200,400]])
  z=np.array([10,200,400])
#  print color_ratio(a,b)
  rain(a,a,z,z,a)
