#!/usr/bin/python

import ReadRaw
from datetime  import datetime, timedelta
import os.path
from ConfigParser import SafeConfigParser

################# Parameters ###############################################
cfgpath       = "../"                                                      # Path with cfg file
cfgfile       = "parameters.cfg"                                           # CFG file
station_block = 'Comodoro'                                                 # Configuration block
path          = "/home/leonardo/Doctorado/Japan2017/desarrollos/Data/"     # Path with input binary data
t1            = datetime(year=2015, month=2, day=15, hour=13, minute=0)    # Start date and time
t2            = datetime(year=2017, month=2, day=17, hour=18, minute=0)    # Final date and time
Debug         = True
############################################################################

config = SafeConfigParser()
config.read(cfgpath+cfgfile)

#Read parameters from an external file
block      = station_block
sampling   = config.getint(block, "sampling")
prefix     = config.get(block, "prefix")
ncpath_raw = config.get("Paths", "ncpath_raw")
ncfile_raw = config.get("Paths", "ncfile_raw")

print "Program for reading binary data from LICEL"
print "Output file (NetCDF format): {}".format(prefix+ncfile_raw)
t1_data, t2_data = ReadRaw.findtimes(t1,t2,sampling,path,prefix,ncpath_raw+prefix+ncfile_raw)
print "Time range: {}-{}".format(t1_data, t2_data)
if t2_data>t1_data:
  x, y1, y2, y3, z, height = ReadRaw.get_data(t1_data,t2_data,sampling,path,prefix)
  if os.path.isfile(ncpath_raw+prefix+ncfile_raw):
    print "Updating file: ", prefix+ncfile_raw
    ReadRaw.updatencd(x, y1, y2, y3, z, ncpath_raw+prefix+ncfile_raw)
  else:
    print "Creating a new file: ", prefix+ncfile_raw
    ReadRaw.createncd(x, y1, y2, y3, z, height, ncpath_raw+prefix+ncfile_raw)
else:
  print "Nothing to do"
  print "Waiting for new data"
print "******** Done ! **********"


