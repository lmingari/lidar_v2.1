#!/usr/bin/python

import numpy as np
from parse import parse
import struct

################# Parameters ####################
Debug = True
#################################################

class GlobalParameters(object):
  def __init__(self,line):
    self.Name = line.strip()

  def line2(self,line):
    s=parse("{} {:tg} {:tg} {:g} {:f} {:f} {:g}",line.strip())
    self.Station     = s[0]
    self.Start       = s[1]
    self.End         = s[2]
    self.HeightASL   = s[3]
    self.Longitude   = s[4]
    self.Latitude    = s[5] # OJO CON ESTO
    self.ZenithAngle = s[6]

  def line3(self,line):
    s=parse("{:d} {:d} {:d} {:d} {:d}",line.strip())
    self.Laser1Shots     = s[0]
    self.Laser1Frequency = s[1]
    self.Laser2Shots     = s[2]
    self.Laser2Frequency = s[3]
    self.Channels        = s[4]

class Channel(object):
  def __init__(self,line):
    s=parse("{:d} {:d} {:d} {:d} {:d} {:g} {:f} {:d}.{:w} 0 0 00 000 {:d} {:d} {:g} {:w}{:d}",line.strip())
    self.isActive         = s[0]==1
    self.isPhotonCounting = s[1]==1
    self.LaserNumber      = s[2]
    self.Bins             = s[3]
    self.isHV             = s[4]
    self.Votage           = s[5]
    self.BinSize          = s[6]
    self.Wavelength       = s[7]
    self.isPolarized      = s[8]=='p'
    self.ADC              = s[9]
    self.Shots            = s[10]
    if self.isPhotonCounting: 
      self.Scale          = s[11]
    else: 
      self.Scale          = s[11] * 1000.0 
    self.Transien         = s[13]

  def set_data(self, data):
    if self.isPhotonCounting:
      ScaleFactor = 150/self.BinSize                         # Signal [MCPS] (Mega Counts Per Second)
    else:
      ScaleFactor = self.Scale / 2**self.ADC                 # Signal [mV].
    if self.Shots==0:
      self.Shots=301
      if Debug: print "Asumiendo {} shots".format(self.Shots)
    self.Signal = np.array(data) * ScaleFactor / self.Shots  # Signal [MCPS] or [mV]
    self.Range  = ( np.arange(self.Bins)+1 ) * self.BinSize  # Range [m]

class LoadLicel(object):
  """
  Description:
  Function that Loads Licel data recorded by Licel VI's
 
  Input Parameters:
  FileName: The LICEL File Name
 
  Output Structure:
 
  data
      |
      |
      |_ GlobalParameters
      |                 |_ HeightASL       : Station Height             [m]
      |                 |_ Latitude        : Station Latitude         [deg]
      |                 |_ Longitude       : Station Longitude        [deg]
      |                 |_ ZenithAngle     : Laser Zenith Angle       [deg]
      |                 |_ Laser1Shots     : Number of Acquired Shots   [-]
      |                 |_ Laser1Frequency : Laser Repetition Rate     [Hz]
      |                 |_ Laser2Shots     : Number of Acquired Shots   [-]
      |                 |_ Laser2Frequency : Laser Repetition Rate     [Hz]
      |                 |_ Channels        : Active Channels            [-]
      |
      |_ Channel
                |_ isActive         : Is it active? (T/F)               [-] 
                |_ isPhotonCounting : Is PhCount (T) or Analog (F)?     [-]
                |_ LaserNumber      : To which Laser is it associated?  [-]
                |_ Bins             : Number of acquired bins           [-]
                |_ isHV             : (T) for HV on                     [-]
                |_ Votage           : Voltage of the PMT / APD          [V]
                |_ BinSize          : Size of the bin                   [m]
                |_ Wavelength       : Detection wavelength             [nm]
                |_ ADC              : Number of eq. bits of the ADC     [-]
                |_ Shots            : Number of acquired shots          [-]
                |_ Scale            : Voltage scale or Threshold level [mV]
                |_ Transient        : Associated transient recorder     [-]
                |_ Signal           : Signal                 [mV] or [MCPS]
                |_ Range            : Altitude Scale                [m AGL]
                |_ Time             : Time Scale from first record  [hours]
 
  or data = -1 if there is an error
  """
  def __init__(self, filename):
    self.channel          = []
    try:
        with open(filename, mode='rb') as file:
          self.globalparameters = GlobalParameters(file.readline())
          self.globalparameters.line2(file.readline())
          self.globalparameters.line3(file.readline())
          for i in range(self.globalparameters.Channels):
            new_channel = Channel(file.readline())
            self.channel.append(new_channel)
          file.seek(80*(self.globalparameters.Channels+3))
          for i in range(self.globalparameters.Channels):
            file.seek(2,1) # Skip CR/LF
            n    = self.channel[i].Bins
            data = struct.unpack('i'*n, file.read(4*n))
            self.channel[i].set_data(data)
    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror)

if __name__ == "__main__":
  import matplotlib.pyplot as plt
  path="/home/leonardo/Downloads/"
  filename="b1732900.011289"
  a = LoadLicel(path + filename)
  print a.globalparameters.Channels
  print np.max(a.channel[0].Signal)
  print (a.channel[0].Signal).shape
