# $Source: /home/rgm/soft/ukidss/listdriven/RCS/listdriven_bulk.py,v $
# $Id: listdriven_bulk.py,v 1.1 2010/10/21 14:28:26 rgm Exp rgm $    
""" 

copied to casu_listdriven_batch.py

except IOError as (errno, strerror):
    print "I/O error({0}): {1}".format(errno, strerror)


  WISHLIST: 
   could split mef pawprints into single extensions to speed up by
   reducing unecessary processing.

  FIX/trap some of the crash conditions; maybe email alert option

  verbose levels; turn off/reduce the ssh debug levels

  HISTORY:

  Original: Richard McMahon, October, 2010

  2010-Oct-01(rgm): added trap for corrupted input files and errlog
    using solution from Eduardo.
  2010-Oct-11(rgm): added table_stats
  2010-Oct-21(rgm): added retry loop for scp
  2010-Oct-23(rgm): added a lockfile to allow multiple versions to run
  2012-Feb-20(rgm): added a try/except/pass to avoid some race conditions
                    that occur when multiple instances are running.

  $Id: listdriven_bulk.py,v 1.1 2010/10/21 14:28:26 rgm Exp rgm $     

  $Log: listdriven_bulk.py,v $
  Revision 1.1  2010/10/21 14:28:26  rgm
  Initial revision

  Login name of author of last revision:   $Author: rgm $ 
  Date and time (UTC) of revision:         $Date: 2010/10/21 14:28:26 $
  Login name of user locking the revision: $Locker: rgm $ 
  CVS revision number:                     $Revision: 1.1 $ 

"""

__version__ = "$Revision: 1.1 $"
__date__ = "$Date: 2010/10/21 14:28:26 $"

import os, sys, time, socket
import re, string
import cookielib, urllib, urllib2
import shlex, subprocess

from time import strftime, gmtime, sleep

from glob import glob
from optparse import OptionParser

import MultipartPostHandler

import pyfits
import splitfile 
import AstroUtils
from logger import *
from flock import flock

from table_stats import *
from pause import *

chunksize = 50000
nfilesmax = 1000
nearest=1

parser = OptionParser()

from time import strftime, gmtime, sleep

pid=os.getpid()
hostname=socket.gethostname()

time_str = strftime("%Y-%m-%dT%H-%M-%S", gmtime())
flog = open("Log_wsa_crossid_bulk_%s.txt" % time_str, 'wt')
flogerr = open("Logerr_wsa_crossid_bulk_%s.txt" % time_str, 'wt')

logdata='Start logfile'
logger(flog, logdata)

logdata='Start error logfile'
logger(flogerr, logdata)

logdata='hostname: ' + hostname
logger(flog, logdata)

logdata='pid: ' + str(pid)
logger(flog, logdata)

logdata="__name__ = " + __name__ 
logger(flog, logdata)

logdata='Current working directory: %s' % (os.getcwd())
logger(flog, logdata)

logdata='Executing: %s' % sys.argv[0]
logger(flog, logdata)

logdata='Version: %s %s %s ' % (sys.argv[0], __version__, __date__)

parser.add_option("-i", "--input", dest = "infile", 
 help = "Input filename")

parser.add_option("--listpath", dest = "listpath", 
 default = './',
 help = "Field list file path")

parser.add_option("-o", "--outpath", dest = "outpath", 
 default = './test/',
 help = "Output path")

parser.add_option("--debug", action = "store_true", dest = "debug", 
 default = False, 
 help = "print status messages to stdout")

parser.add_option("--selftest", action = "store_true", dest = "selftest", 
 default = False, 
 help = "run listphot using imcore catalogue from image")

parser.add_option("-v", "--verbose", 
 action = "store_true", dest = "verbose", 
 default = False,
 help = "print status messages to stdout")

parser.add_option("-w", "--wait", 
 dest = "wait", default = "0",
 help = "wait time between queries in seconds (0 = wait until query finishes) [0]")

parser.add_option("--clobber", dest = "clobber", 
 help = "Clobber or overwrite existing results")

#parser.add_option("-l", "--logfile", dest = "logfile", 
# help = "Logfile to log progress and summary")

(options, args) = parser.parse_args()

debug=options.debug

logfile = 'wsa_crossid_bulk.logfile'

logdata = 'options: '
logger(flog, logdata)

logdata = options
logger(flog, logdata)

logdata = 'args: '
logger(flog, logdata)

logdata = args
logger(flog, logdata)

infile = options.infile
logdata='infile: %s ' % infile
logger(flog, logdata)

outpath = options.outpath + '/'
logdata='outpath: %s ' % outpath
logger(flog, logdata)

# count the number lines in the input file
numLines = len(open(infile).readlines())
logdata = 'Input file has length: %d' % numLines
logger(flog, logdata)


# Read field list file skipping comment lines
records = [item for item in open(infile) if item[0]<>'#']
numFields = len(records)
logdata = 'Total number of fields: %d' % numFields
logger(flog, logdata)

starttime=time.time()

SelfTest = options.selftest

listpath=options.listpath 
ipos=listpath.rfind('/')
print ipos, len(listpath)
if ipos+1 < len: listpath=options.listpath + '/'
logdata='listpath: %s' % listpath 
logger(flog, logdata)

ifile=0
DryRun = False
nfiles=len(records)
#help(records)
#key=raw_input("Enter any key to continue: ")

logdata = "Total elapsed time %.3f seconds" % (time.time() - starttime)
logger(flog, logdata)

# loop through line by line
step=1
start=step-1
records=records[start:len(records):step]
verbose=False
n_errors=0
for line in records:
    ifile = ifile + 1
    line_save = line
    line = line.split()
    files = [line[i] for i in [3] if line[i]<>'null']
    filename = files[0]

    # find first fil
    if SelfTest:
      catfile = files[0][0:-4] + '_cat.fits'
      ipos1=catfile.rfind('/')
      ipos2=catfile.rfind('.fits')
      outfile=outpath+catfile[ipos1+1:ipos2] + '_listdriven.fits'
      if os.path.exists(outfile):
        logdata='List-driven data already exists for %s.' % (outfile)
        logger(flog, logdata)
        continue

    print 
    logdata= "Processing field: ", ifile, ' of ',nfiles
    logger(flog, logdata)

    logdata= 'Image filename: ' + filename
    logger(flog, logdata)

    deltatime=time.time()

    confname = files[0][:-4]+'_conf.fit'
    if verbose: 
      logdata='Confidence map: %s' % confname
      logger(flog, logdata)

    # strip of the path from the filename
    ipos=files[0].rfind('/')
    listfile = files[0][ipos+1:] + '.radec'
    if verbose: logdata='List file: %s' % listfile
    if verbose: logger(flog, logdata)

    if not SelfTest:
      ipos=files[0].rfind('/')
      listfile = files[0][ipos+1:] + '.radec'
      outfile=outpath + listfile + '_listdriven.fits'
      logdata='Output file: %s' % outfile
      logger(flog, logdata)
      listfile = listpath + listfile
      if os.path.exists(outfile):
        print 'Skipping since outfile already exists for %s.' % (outfile)
        continue

    if SelfTest:
      logdata= 'Running self test regression using imcore catalogue'   
      logger(flog, logdata)
      catfile = files[0][0:-4] + '_cat.fits'
      logdata='catalogue file: %s' % catfile
      logger(flog, logdata)
      ipos1=catfile.rfind('/')
      ipos2=catfile.rfind('.fits')
      outfile=outpath+catfile[ipos1+1:ipos2] + '_listdriven.fits'
      print 'outfile: ', outfile
      if os.path.exists(outfile):
        print 'List-driven data already exists for %s.' % (outfile)
        continue

    lockfile=outfile + '.lock' 
    if os.path.exists(lockfile):
      logdata='Skipping since lockfile exists: ' + lockfile
      logger(flog, logdata) 
      logdata = "Total elapsed time %.3f seconds" % (time.time() - starttime)
      logger(flog, logdata)
      continue

    # create lockfile
    lkfile = open(lockfile, 'wt')
    logdata = "Create lockfile %s" % lockfile
    logger(flog, logdata)

    # write pid and hostname into file
    #lockfile=listfile + '.lock.' + str(pid)
    #logdata = "Lockfile: " + lockfile
    #logger(flog, logdata)
    #lock = flock(lockfile, True).acquire()

    if not SelfTest and os.access(listfile, os.F_OK): 
      # Read list file skipping comment lines
      records = [item for item in open(listfile) if item[0]<>'#']
      numSources = len(records)
      logdata = 'Number of sources in listfile: %d' % numSources
      logger(flog, logdata)

    if not SelfTest and not os.access(listfile, os.F_OK): 
      logdata='List file problem'
      logger(flog, logdata)
      key=raw_input("Enter any key to continue: ")

    # build the scp command; note spaces between each parameter
    # use -v to debug the scp/ssh 
    #host='rgm@apm14.ast.cam.ac.uk:'
    #host='rgm@apm26.ast.cam.ac.uk:'
    host='rgm@apm25.ast.cam.ac.uk:'
    cmd ='scp -v -i ~/.ssh/id_dsa_nk ' \
      + ' '+ host + filename \
      + ' ' + outpath +'/. '

    logdata=cmd
    logger(flog, logdata)

    itry=0
    itry_max=10
    Transfered=False
    while (itry < itry_max) and not Transfered:
      itry=itry+1
      if debug: print cmd
      result=os.popen(cmd)
      while True:
        line = result.readline()
        if line == "": break
        logdata = line
        logger(flog, logdata)
        print line,
        if debug: key=raw_input("Enter any key to continue: ")

      # check the files was transfered
      ipos=filename.rfind('/')
      imagefile = outpath + filename[ipos+1:] 

      if os.access(imagefile, os.F_OK):
        Transfered=True

      if not os.access(imagefile, os.F_OK):
        logdata='WARNING: image file NOT transfered: %s' % imagefile
        logger(flogerr, logdata)
        Transfered=False
        delay=10*(itry*itry)
        logdata='WAITING: %s seconds' % delay
        logger(flog, logdata)
        time.sleep(delay)
        continue

      cmd ='scp -i ~/.ssh/id_dsa_nk ' \
       + ' '+ host + confname \
       + ' ' + outpath +'/. '

      logdata=cmd
      logger(flog, logdata)

      if debug: print cmd
      result=os.popen(cmd)
      while True:
        line = result.readline()
        if line == "": break
        logdata = line
        logger(flog, logdata)
        print line,
        if debug: key=raw_input("Enter any key to continue: ")

    # check the files was transfered
    ipos=confname.rfind('/')
    confmapfile=outpath+confname[ipos+1:]
    if not os.access(confmapfile, os.F_OK):
      n_errors=n_errors+1
      logdata= 'WARNING: confmap file NOT transferred: ', confmapfile
      logger(flog, logdata)
      logger(flogerr, logdata)
      if os.path.exists(lockfile):
        logdata = "Delete lockfile %s " % lockfile
        logger(flog, logdata)
        os.remove(lockfile)       
      continue

    if SelfTest:
      print 'SelfTest catalogue file: ' + catfile
      cmd ='scp -i ~/.ssh/id_dsa_nk ' \
        + ' '+ host + catfile \
        + ' ' + outpath +'/. '

      logdata=cmd
      logger(flog, logdata)

      if debug: print cmd
      result=os.popen(cmd)
      while True:
        line = result.readline()
        if line == "": break
        logdata = line
        logger(flog, logdata)
        print line,
        if debug: key=raw_input("Enter any key to continue: ")
      print 'SelfTest File transferred: ' + catfile

      ipos=catfile.rfind('/')
      catfile=outpath+catfile[ipos:]
      if os.access(catfile, os.F_OK):
        print 'catfile transfered OK: ',catfile
      if not os.access(catfile, os.F_OK):
        print 'WARNING: catfile NOT transfered: ',catfile
        key=raw_input("Enter any key to continue: ")

      listfile=catfile


    logdata = "Delta elapsed time %.3f seconds" % (time.time() - deltatime)
    logger(flog, logdata)

    if not DryRun:
      logdata= 'Start processing the image data'
      logger(flog, logdata)

      if debug: key=raw_input("Enter any key to continue: ")

      nustep=pyfits.getval(imagefile,'NUSTEP',0)
      print 'nustep = ',nustep
  
      if (nustep==1):
        rcore = "2.5"
        nbsize = "64"
        threshold = "1.5"
      elif (nustep==4): # i.e., 2x2
        rcore = "5.0"
        nbsize = "128"
        threshold = "1.25"
      elif (nustep==9): # i.e., 3x3
        rcore = "7.5"
        nbsize = "192"
        threshold = "1.25"

      print 'rcore, nbsize, threshold: ', rcore, nbsize, threshold

      cmd ='imcore_list ' \
        + ' ' + imagefile \
        + ' ' + confmapfile \
        + ' ' + listfile \
        + ' ' + outfile \
        + ' ' + threshold  \
        + ' --nbsize=' + nbsize \
        + ' --rcore=' + rcore \
        + ' --cattype=6 ' 

      # + ' --verbose '

      stdoutlog = open('logfile_stdout', 'w+')
      stderrlog = open('logfile_stderr', 'w+')

      command=cmd
      logdata=command
      logger(flog, logdata)

      args = shlex.split(command)
      print args
      p = subprocess.call(args, stderr=stderrlog, stdout=stdoutlog)
      logdata='subprocess error status: ' + str(p)
      logger(flog, logdata)
      if p is not 0:
	print 'Something went wrong: ', args
        logdata = 'Something went wrong: ' + command
        logger(flogerr, logdata)
        logger(flogerr, line_save)
        logger(flogerr, command)
        # delete the outfile if created
        if os.path.exists(outfile):
          logdata = "Delete %s " % outfile
          logger(flog, logdata)
          logger(flogerr, logdata)  
          os.remove(outfile)       
        if os.path.exists(lockfile):
          logdata = "Delete lockfile %s " % lockfile
          logger(flog, logdata)
          logger(flogerr, logdata)
          os.remove(lockfile)       
        continue
      else:
	logdata='imcore_list Finished'
        logger(flog, logdata)

      hdulist = pyfits.open(outfile)
      logdata= 'Number of extensions: %d ' % len(hdulist)
      logger(flog, logdata)

      n_ext=len(hdulist)
      for ext in range(1, n_ext):
        table_stats(outfile, ext=ext)

      print 'listdriven photometry completed'
      #key=raw_input("Enter any key to continue: ")

    print 'Deleting data files used'
    if os.path.exists(imagefile):
      print 'Remove the image file:' + imagefile
      cmd ='rm -v ' + imagefile
      try:
        result=os.popen(cmd)
      except OSError as (errno, strerror):
        logdata ="OS error({0}): {1}".format(errno, strerror)
        logger(flogerr, logdata)
        logdata = "error removing imagefile %s " % imagefile        
        logger(flogerr, logdata)
        pass

    if os.path.exists(confmapfile):
      print 'Remove the confidence map:' + confmapfile
      cmd ='rm -v ' + confmapfile
      try:
        result=os.popen(cmd)
      except:
        logdata = "error removing confmapfile %s " % confmapfile        
        logger(flogerr, logdata)
        pass

    if SelfTest:
      if os.path.exists(catfile):
        print 'Remove the catalogue fits file:' + catfile
        cmd ='rm -v ' + catfile
        try:
          result=os.popen(cmd)
        except:
          logdata = "error removing cataloge file%s " % catfile        
          logger(flogerr, logdata)
          pass

    if os.path.exists(lockfile):
      logdata = "Delete lockfile %s " % lockfile
      logger(flog, logdata)
      cmd ='rm -v ' + lockfile
      try:
        os.remove(lockfile)          
      except OSError as (errno, strerror):
        logdata ="OS error({0}): {1}".format(errno, strerror)
        logger(flogerr, logdata)
        logdata = "error removing lockfile%s " % lockfile        
        logger(flogerr, logdata)
        pass

    logdata = "Delta elapsed time %.3f seconds" % (time.time() - deltatime)
    logger(flog, logdata)

    logdata = "Total elapsed time %.3f seconds" % (time.time() - starttime)
    logger(flog, logdata)

    if debug: key=raw_input("Enter any key to continue: ")
  
if debug: print 'debug on'

elapsedtime = time.time() - starttime
logdata = "Elapsed time: %.3f seconds" % (time.time() - starttime)
logger(flog, logdata)

debug=False

logdata = "Query completed"
flog.close()
flogerr.close()
