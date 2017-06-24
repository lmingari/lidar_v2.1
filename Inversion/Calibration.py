#!/usr/bin/python
import numpy as np
import Plots
from datetime import datetime

################# Parameters ####################
DoPlot = False
Debug  = True
#################################################

def vis_factor(vis,z):
  DZ    = z[1]-z[0]
  DZinv = 1.0/DZ

  iz_1200  = int(round(1.2 * DZinv))-1
  iz_6000  = int(round(6.0 * DZinv))-1

  with np.errstate(divide='ignore', invalid='ignore'):
    a = np.where(vis[iz_1200:iz_6000,:] > 0, np.log10(vis[iz_1200:iz_6000,:]), np.nan)
  hist, bin_edges = np.histogram(a[~np.isnan(a)], bins=60, range=[-3,2])
  maxsub          = np.nanargmax(hist)
  factor          = 1E4 / 10.**(bin_edges[maxsub])
  if Debug: 
    print "**** Function vis_factor ****"
    print "Index of maximum in histogram: {}".format(maxsub)
    print "Maximum histogram: {}".format(bin_edges[maxsub])
  if DoPlot: Plots.histo(hist, bin_edges)
  return factor
  
def color_ratio(vis, ir, z, maxint):
  output  = 1.0
  DZinv   = 1.0/(z[1]-z[0])
  iz_600  = int(round(0.6 * DZinv))-1
  iz_3000 = int(round(3.0 * DZinv))-1

  sub_vis = vis[iz_600:iz_3000,:]
  sub_ir  =  ir[iz_600:iz_3000,:]

  threshold = (DZinv * 6E-3) * vis.shape[-1]/96.0
  with np.errstate(invalid='ignore'): 
    cl  = (sub_vis>=maxint).sum()
  if cl > threshold:
    with np.errstate(invalid='ignore'): 
      x = sub_vis[sub_vis>=maxint]
      y = sub_ir [sub_vis>=maxint]
    params = np.polyfit(x,y,1)
    output = params[0]
    if Debug: 
      print "**** Function color_ratio ****"
      print "Number of Cloud points: {} - Color Ratio: {}".format(cl, output)
      print "Plotting channel relation..."
    if DoPlot: Plots.xy(x, y, params, "calibration.png")
  return output
  
def depol(par, per, z):
  output  = 1.0
  DZinv   = 1.0/(z[1]-z[0])
  iz_12000 = int(round(10.0 * DZinv))-1

  x = par[iz_12000,:]
  y = per[iz_12000,:]

  params = np.polyfit(x[~np.isnan(x)],y[~np.isnan(x)],1)
  output = params[0]
  if Debug: 
    print "**** Function depol ****"
    print "Linear fit: slope={}".format(output)
  if DoPlot: Plots.xy(x, y, params, "depol.png")
  return output
  
def power(x, filename):
  print "**** Function power ****"
  NX     = len(x)
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

def overlap(filename,z):
  yr     = np.ones_like(z, dtype=np.float)
  dtype1 = np.dtype([('height', 'f8'), ('yr', 'f8')])
  if Debug:
    print "**** Function overlap ****"
    print "Overlap File: {}".format(filename)
    print "Checking file consistence..."
  try:
    data   = np.loadtxt(filename, dtype=dtype1)
    if len(data)==len(z):
      if Debug: print "OK"
      for i in range(len(data)): yr[i]=data['yr'][i]
    else:
      if Debug: print "Incorrect number of levels"
  except IOError:
    if Debug: print "Overlap File not found"
  return yr

if __name__ == "__main__":
  vis = np.array([[1E5,2E5,1.2E5,1.45E5,np.nan,1.45E5],[1E5,2E5,1.2E5,1.45E5,1.2E5,1.45E5]])
  z   = [2,3, 120, 125, 142.5]
#  print histo_param(vis.T,z)
  filename = "YR.txt"
  yr = overlap(filename,z)
  print 1.0 / yr
