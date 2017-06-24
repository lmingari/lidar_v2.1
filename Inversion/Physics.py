#!/usr/bin/python
# Extended to 18km
# Altitude / CIRA86
import numpy as np

################# Parameters ####################
s2 = 8.0 * np.pi / 3.0                          #
Na = 6.022E23                                   #(/mol)
m0 = 28.97                                      #(g/mol)
m0 = m0 / 1e3                                   #(kg/mol)
stdatm = [1.2250, 1.1117, 1.0066, 9.0925e-1, \
          8.1935e-1, 7.3643e-1, 6.6011e-1,   \
          5.9002e-1, 5.2579e-1, 4.6706e-1,   \
          4.1351e-1, 3.6480e-1, 3.1194e-1,   \
          2.6660e-1, 2.2786e-1, 1.9476e-1,   \
          1.6647e-1, 1.4230e-1, 1.2165e-1]      #US standard atmosphere (kg/m3)
#################################################
             

def rayleigh(z, wavelength, alt):
### z, alt     = meters
### month      = 1-12
### lat        = degree N
### wavelength = nm
  rho        = np.interp(z + alt, np.linspace(0,18, num=19) * 1e3, stdatm)
  Ng         = Na * rho / m0                                                    #(/m^3)
  crosssec_r = 5.45 * (550. / wavelength) ** 4 * 1e-28 * 1e-4
  beta_r     = crosssec_r * Ng
  alpha_r    = beta_r * s2
  return alpha_r, beta_r

def fernald(x, z, alpha_r, beta_r, s1, top, init_beta):
  """Fernald Algorithm
  Keyword arguments:
  x: range-corrected intensity
  alpha_r = extinction coefficient of molecules
  beta_r = backscattering coefficient of molecules
  s1 = lidar ratio
  ndata = number of data points
  resol = vertical resolution
  top = starting point
  init_beta = initial value of beta at TOP (optional)
  """
  NZ    = len(z)
  DZ    = z[1]-z[0]
  DZinv = 1.0/DZ
  DZ    = 1000.0 * DZ # in meters
  
  x = x * 1E-4

  iz_top    = int(round(top * DZinv))-1
  ndata     = iz_top

  alpha     = np.zeros(NZ)
  beta      = np.zeros(NZ)
  beta[iz_top-1]  = init_beta
  alpha[iz_top-1] = init_beta * s1

  for iz in reversed(range(1,iz_top)):
    a = (s1-s2)*(beta_r[iz-1]+beta_r[iz])*DZ
    num = x[iz-1] * np.exp(a)
    den = x[iz]/(beta[iz]+beta_r[iz]) 
    den = den + s1*(x[iz]+x[iz-1]*np.exp(a))*DZ
    if den>0:
      beta[iz-1] = num/den - beta_r[iz-1]
    else:
      beta[iz-1] = 0.0

  alpha = s1 * beta
  return alpha, beta

if __name__ == "__main__":
  z = np.array([0, 100, 2000, 5500])
  a,b = rayleigh(z, 1064.0, 0.0)
  print a
  print b
