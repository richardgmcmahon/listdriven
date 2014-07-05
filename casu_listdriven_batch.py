# $Source: /home/rgm/soft/ukidss/listdriven/RCS/listdriven_bulk.py,v $
# $Id: listdriven_bulk.py,v 1.1 2010/10/21 14:28:26 rgm Exp rgm $    
""" 



  WISHLIST: 

   could look at option to cache the image files on /tmp if this
   would give speed improvement; call stagepath 

   Need to trap the pyfits errors:
   except IOError as (errno, strerror):
   print "I/O error({0}): {1}".format(errno, strerror)

   could split mef pawprints into single extensions to speed up by
   reducing unecessary processing.

  FIX/trap some of the crash conditions; maybe email alert option

  verbose levels; turn off/reduce the ssh debug levels


  subprocess.call; subprocess.Popen; os.popen

  Both subprocess.call and subprocess.Popen are used for
  historical reasons due to previous use of os.popen

  google: python subprocess call popen

  See:
http://stackoverflow.com/questions/7681715/whats-the-difference-between-subprocess-popen-and-call-how-do-you-use-them-to

http://stackoverflow.com/questions/6934696/python-subprocess-call-and-subprocess-popen-giving-me-different-outputs

subprocess.Popen returns a Popen object, which you can use to communicate with the process and get output, however subprocess.call will only return the return code of the process:

So subprocess.call is basically equivalent to the following code, and only exists for convenience:



def call(*popenargs, timeout=None, **kwargs):
    
    Run command with arguments.  Wait for command to complete or
    timeout, then return the returncode attribute.

    The arguments are the same as for the Popen constructor.  Example:

    retcode = call(["ls", "-l"])


    with Popen(*popenargs, **kwargs) as p:
        try:
            return p.wait(timeout=timeout)
        except:
            p.kill()
            p.wait()
            raise


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

from __future__ import print_function, division

import os
import errno

def make_sure_path_exists(path):
    """
    see:
http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary

    this should catch a race condition where the directory gets created 
    between the check and the creation
    see the lengthy discussion

    """
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def imcore_list_run(options=None, imagefile=None,
  confmapfile=None, listfile=None, outfile=None):

  # set the imcorelist parameters; maybe need to check these againsts
  # values in the header of the catalogue file
  nustep=-1
  print('options.vhs: ',options.vhs)
  print('options.cache: ',options.cache)

  if not options.vhs:
    nustep=pyfits.getval(imagefile,'NUSTEP',0)
    print('nustep = ', nustep)
  
  if options.vhs: 
    rcore = "3.0"
    nbsize = "64"
    threshold = "1.5"

  # UKIDSS values
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

  print('rcore, nbsize, threshold: ', rcore, nbsize, threshold)

  # needs the binary location to be a parameter for portability
  command ='nice -n19 /home/rgm/bin/imcore_list ' \
        + ' ' + imagefile \
        + ' ' + confmapfile \
        + ' ' + listfile \
        + ' ' + outfile \
        + ' ' + threshold  \
        + ' --nbsize=' + nbsize \
        + ' --rcore=' + rcore \
        + ' --cattype=6 ' 

      # + ' --verbose '

  stdoutlog = open(logpath+'Logfile_stdout', 'w+')
  stderrlog = open(logpath+'Logfile_stderr', 'w+')

  logdata=command
  logger(flog, logdata)

  args = shlex.split(command)
  print('Print command as args: ')
  print(args)
  result = subprocess.call(args, \
       stderr=stderrlog, stdout=stdoutlog)
  logdata='subprocess error status: ' + str(result)
  logger(flog, logdata)
  trace = traceback.extract_stack()[-1]
  print(os.path.basename(trace[0]), ' line :', str(trace[1]))
  if result is not 0: print(len(result))
  print(result)

  if result is not 0:
    print(result)
    trace = traceback.extract_stack()[-1]
    print(os.path.basename(trace[0]), ' line :', str(trace[1]))
    print('Something went wrong: ', args)
    logger(flogerr, logdata)
    logdata = 'Something went wrong: ' + command
    logger(flogerr, logdata)
    logger(flogerr, line_save)
    logger(flogerr, command)
   
    while True:
      line = result.readline()
      if line == "": break
      logdata = line
      logger(flog, logdata)
      print(line,)
      if debug: key=raw_input("Enter any key to continue: ")

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
    #continue
  
    usage = resource.getrusage(resource.RUSAGE_CHILD)
    for name, desc in [
        ('ru_utime', 'User time'),
        ('ru_stime', 'System time'),
        ]:
        print('%-25s (%-10s) = %s' % (desc, name, getattr(usage, name)))


  else:
    logdata='imcore_list Finished'
    logger(flog, logdata)

  print('Read back the results and check integrity: ', outfile)
  if not os.path.exists(outfile):       
    print('File does not exist: ', outfile)

  try:
    hdulist = pyfits.open(outfile)
  except:
    traceback.print_exc(file=sys.stdout)

  logdata= 'Number of extensions: %d ' % len(hdulist)
  logger(flog, logdata)

  n_ext=len(hdulist)
  for ext in range(1, n_ext):
    table_stats(outfile, ext=ext)

  print('listdriven photometry completed')
  #key=raw_input("Enter any key to continue: ")


def get_vista_pawprints(imagefile=None, filename=None, 
 stagepath=None, SelfTest=False):
  """
  Under developement option to allow pawprint processing
  """
  # basic scp with out retry

  print('filename: ', filename)

  pathname = os.path.dirname(filename) 

  hdulist = pyfits.open(imagefile)
  header=hdulist[1].header
  print()
  print('PROV files:')
  #print(header['PROV0000'])
  list=[header['PROV0001']]
  list.append(header['PROV0002'])
  list.append(header['PROV0003'])
  list.append(header['PROV0004'])
  list.append(header['PROV0005'])
  list.append(header['PROV0006'])
  #print(header['PROV0001'])
  #print(header['PROV0002'])
  #print(header['PROV0003'])
  #print(header['PROV0004'])
  #print(header['PROV0005'])
  #print(header['PROV0006'])
  #print()
  print(list)

  #key=raw_input("Enter any key to continue: ")

  for image in list:

    filename = pathname + '/' + image
    command ='scp  -i ~/.ssh/id_dsa_nk ' \
      + ' '+ host + filename \
      + ' ' + stagepath +'/. '

    logdata=command
    logger(flog, logdata)

    #result=os.popen(command)
    result=subprocess.Popen(command)  

    # create confidence map filename from the image filename
    confname = filename[:-4]+'_conf.fit'

    command ='scp  -i ~/.ssh/id_dsa_nk ' \
      + ' '+ host + confname \
      + ' ' + stagepath +'/. '

    logdata=command
    logger(flog, logdata)

    #result=os.popen(command)
    result=subprocess.Popen(command)  

    if SelfTest:
      catfile = files[0][0:-4] + '_cat.fits'

      print('SelfTest catalogue file: ' + catfile)
      command ='scp -i ~/.ssh/id_dsa_nk ' \
        + ' '+ host + catfile \
        + ' ' + stagepath +'/. '

      logdata=command
      logger(flog, logdata)

      if debug: print(command)

      #result=os.popen(command)
      result=subprocess.Popen(command)  
      while True:
        line = result.readline()
        if line == "": break
        logdata = line
        logger(flog, logdata)
        print(line,)
        if debug: key=raw_input("Enter any key to continue: ")
      print('SelfTest File transferred: ' + catfile)

    #key=raw_input("Enter any key to continue: ")

  return list


def process_image(files=None, outpath=None, Search=False):

    catfile = files[0][0:-4] + '_cat.fits'
    logdata='process_image: Number of files = %s' % (str(len(files)))
    logger(flog, logdata)

    # find first file; needs to be refactored since only one file is input
    if SelfTest:
      ipos1=catfile.rfind('/')
      ipos2=catfile.rfind('.fits')
      outfile=outpath+catfile[ipos1+1:ipos2] + '_listdriven.fits'
      if os.path.exists(outfile):
        logdata='List-driven data already exists for %s.' % (outfile)
        logger(flog, logdata)
        return

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
        print('Skipping since outfile already exists for %s.' % (outfile))
        return
        

    if SelfTest:
      logdata= 'Running self test regression using imcore catalogue'   
      logger(flog, logdata)
      catfile = files[0][0:-4] + '_cat.fits'
      logdata='catalogue file: %s' % catfile
      logger(flog, logdata)
      ipos1=catfile.rfind('/')
      ipos2=catfile.rfind('.fits')
      outfile=outpath+catfile[ipos1+1:ipos2] + '_listdriven.fits'
      print('outfile: ', outfile)
      if os.path.exists(outfile):
        print('List-driven data already exists for %s.' % (outfile))
        return

    lockfile=outfile + '.lock' 
    print('lockfile: ', lockfile)
    if os.path.exists(lockfile):
      logdata='Skipping since lockfile exists: ' + lockfile
      logger(flog, logdata) 
      logdata = "Total elapsed time %.3f seconds" % (time.time() - starttime)
      logger(flog, logdata)
      return

    # create lockfile
    lkfile = open(lockfile, 'wt')
    logdata = "Create lockfile %s" % lockfile
    logger(flog, logdata)

    # write pid and hostname into file
    lkfile.write(strftime("%Y-%m-%dT%H-%M-%S", gmtime()))
    lkfile.write(':pid: '+ str(pid) + '\n')
    lkfile.flush()

    #lockfile=listfile + '.lock.' + str(pid)
    #logdata = "Lockfile: " + lockfile
    #logger(flog, logdata)
    # use flock module
    #lock = flock(lockfile, True).acquire()

    if not SelfTest and os.access(listfile, os.F_OK): 
      # Read list file skipping comment lines
      records = [item for item in open(listfile) if item[0]<>'#']
      numSources = len(records)
      logdata = 'Number of sources in listfile: %d' % numSources
      logger(flog, logdata)

    if not SelfTest and not os.access(listfile, os.F_OK): 
      logdata='List file: %s' % listfile
      logger(flog, logdata)
      logdata='List file problem'
      logger(flog, logdata)
      key=raw_input("Enter any key to continue: ")

    # build the scp command; note spaces between each parameter
    # use -v to debug the scp/ssh 
    #cmd ='scp -v -i ~/.ssh/id_dsa_nk ' \
    command ='scp  -i ~/.ssh/id_dsa_nk ' \
      + ' '+ host + filename \
      + ' ' + stagepath +'/. '

    logdata=command
    logger(flog, logdata)

    itry=0
    Transfered=False

    while (itry < iretry_max) and not Transfered:
      itry=itry+1
      if debug: 
        print(command)
      #result=os.popen(command)
      #help(result)
      result=subprocess.Popen(command, 
       stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
      output, errors = result.communicate()
      #result=subprocess.check_output(command)
      #help(result)
      print('subprocesss output: ', output)
      print('subprocess errors: ', errors)
      #while True:
      #  line = errors.readline()
      #  if line == "": break
      #  logdata = line
      #  logger(flog, logdata)
      #  print(line,)
      #  if debug: key=raw_input("Enter any key to continue: ")

      # check the files was transfered
      ipos=filename.rfind('/')
      imagefile = stagepath + filename[ipos+1:] 

      if os.access(imagefile, os.F_OK):
        Transfered=True

      if not os.access(imagefile, os.F_OK):
        scpfailurefile=outfile + '.scpfailure' 
        scpfailurefileh = open(scpfailurefile, 'wt')

        logdata='WARNING: image file NOT transfered: %s' % imagefile
        logger(flogerr, logdata)
        Transfered=False
        delay=delaytime*(itry*itry)
        logdata='WAITING: %s seconds' % delay
        logger(flog, logdata)
        time.sleep(delay)
        continue

      command ='scp -i ~/.ssh/id_dsa_nk ' \
       + ' '+ host + confname \
       + ' ' + stagepath +'/. '

      logdata=command
      logger(flog, logdata)

      if debug: print(command)
      usage = resource.getrusage(resource.RUSAGE_SELF)
      for name, desc in [
        ('ru_utime', 'User time'),
        ('ru_stime', 'System time'),
        ]:
        print('%-25s (%-10s) = %s' % (desc, name, getattr(usage, name)))

      #result=os.popen(command)
      result=subprocess.Popen(command, 
       stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
      output, errors = result.communicate()
      #help(result)
      print('subprocesss output: ', output)
      print('subprocess errors: ', errors)
      usage = resource.getrusage(resource.RUSAGE_SELF)  
      for name, desc in [
        ('ru_utime', 'User time'),
        ('ru_stime', 'System time'),
        ]:
        print('%-25s (%-10s) = %s' % (desc, name, getattr(usage, name)))

      #while True:
      #  line = result.readline()
      #  if line == "": break
      #  logdata = line
      #  logger(flog, logdata)
      #  print(line,)
      #  if debug: key=raw_input("Enter any key to continue: ")

    # check the files was transfered
    ipos=confname.rfind('/')
    confmapfile=stagepath + confname[ipos+1:]
    if not os.access(confmapfile, os.F_OK):
      n_errors=n_errors+1
      logdata= 'WARNING: confmap file NOT transferred: ', confmapfile
      logger(flog, logdata)
      logger(flogerr, logdata)
      if os.path.exists(lockfile):
        logdata = "Delete lockfile %s " % lockfile
        logger(flog, logdata)
        os.remove(lockfile)       


    if options.pawprints:
      # read the file to determine the constituent pawprints
      get_vista_pawprints(imagefile=imagefile, filename=filename, 
       stagepath=stagepath)

    if SelfTest or Search:
      print('SelfTest/Search catalogue file: ' + catfile)
      if os.path.exists(catfile):
        logdata='catalogue already exists %s.' % (catfile)
        logger(flog, logdata)

      if not os.path.exists(catfile):
        command ='scp -i ~/.ssh/id_dsa_nk ' \
         + ' '+ host + catfile \
         + ' ' + stagepath +'/. '

        logdata=command
        logger(flog, logdata)

        if debug: print(command)
        result=subprocess.Popen(command, 
         stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
        output, errors = result.communicate()
        #result=os.popen(command)
        #result=subprocess.Popen(command)  
        #while True:
        #  line = result.readline()
        #  if line == "": break
        #  logdata = line
        #  logger(flog, logdata)
        #  print(line,)
        #  if debug: key=raw_input("Enter any key to continue: ")
      
      print('SelfTest cataloge file transferred: ' + catfile)

      ipos=catfile.rfind('/')
      catfile= stagepath + catfile[ipos:]
      if os.access(catfile, os.F_OK):
        print('catfile transfered OK: ',catfile)
      if not os.access(catfile, os.F_OK):
        print('WARNING: catfile NOT transfered: ',catfile)
        key=raw_input("Enter any key to continue: ")

      listfile=catfile

    logdata = "Delta elapsed time %.3f seconds" % (time.time() - deltatime)
    logger(flog, logdata)

    if not DryRun:
      logdata= 'Start processing the image data'
      logger(flog, logdata)

      if debug: 
        trace = traceback.extract_stack()[-1]
        print(os.path.basename(trace[0]), ' line :', str(trace[1]))
        key=raw_input("Debug: Enter any key to continue: ")

      imcore_list_run(options=options, imagefile=imagefile,
        confmapfile=confmapfile, listfile=listfile, outfile=outfile)

      print('listdriven photometry completed')
      #key=raw_input("Enter any key to continue: ")

    print('options.cache: ',options.cache)
    if os.path.exists(imagefile) and not options.cache:
      print('Deleting data files used')
      print('Remove the image file:' + imagefile)
      try:
        os.remove(imagefile)       
      except OSError as (errno, strerror):
        logdata ="OS error({0}): {1}".format(errno, strerror)
        logger(flogerr, logdata)
        logdata = "error removing imagefile %s " % imagefile        
        logger(flogerr, logdata)
        pass

    if os.path.exists(confmapfile) and not options.cache:
      print('Remove the confidence map:' + confmapfile)
      try:
        os.remove(confmapfile)       
      except:
        logdata = "error removing confmapfile %s " % confmapfile        
        logger(flogerr, logdata)
        pass

    if SelfTest:
      if os.path.exists(catfile)  and not options.cache:
        print('Remove the catalogue fits file:' + catfile)
        try:
          os.remove(catfile)       
        except:
          logdata = "error removing cataloge file%s " % catfile        
          logger(flogerr, logdata)
          pass

    if os.path.exists(lockfile):
      logdata = "Delete lockfile %s " % lockfile
      logger(flog, logdata)
      command ='rm -v ' + lockfile
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

    



__version__ = "$Revision: 1.1 $"
__date__ = "$Date: 2010/10/21 14:28:26 $"


import os, sys, time, socket
import subprocess
import re, string
import resource

import cookielib, urllib, urllib2
import shlex
import traceback

from time import strftime, gmtime, sleep

from glob import glob
from optparse import OptionParser

import MultipartPostHandler

import pyfits
import splitfile 
import AstroUtils
from logger import *
from flock import flock

import multiprocessing

# use cat_cal(cat, ra, dec, chipno)
# returns 1 line table
# 
#sys.path.append('/home/sr525/Python_Code/')
sys.path.append('/home/rgm/soft/sreed/')
#import cat_cal as srlib 

from table_stats import *
from pause import *

chunksize = 50000
nfilesmax = 1000
nearest=1

ncores=multiprocessing.cpu_count()
print('Number of cores: ', ncores)
nworkers = ncores*4
print('Number of workers: ', nworkers)

# create Queue in module multiprocessing.queues object
work_queue = multiprocessing.Queue()
done_queue = multiprocessing.Queue()

# create list for the processes
processes = []

starttime=time.time() 
computetime=0.0


iretry_max=1
nskip_files=0
delaytime=1.0

import ConfigParser
config = ConfigParser.RawConfigParser()

# read config file
config.read('casu_listdriven_batch.cfg')
host = config.get('casu','host')
user = config.get('casu','user')
#host='rgm@apm21.ast.cam.ac.uk:'
host=user+'@'+host+':'



# OptionParser is deprecated so need to replace with argparse
# see section on Upgrading optparse code in argparse
parser = OptionParser()

from time import strftime, gmtime, sleep

pid=os.getpid()
hostname=socket.gethostname()

time_str = strftime("%Y-%m-%dT%H-%M-%S", gmtime())

parser.add_option("-i", "--input", dest = "infile", 
 help = "Input filename")

parser.add_option("--listpath", dest = "listpath", 
 default = './',
 help = "Field list file path")

parser.add_option("-o", "--outpath", dest = "outpath", 
 default = '.',
 help = "Output path")

parser.add_option("--cache", action = "store_true", dest = "cache", 
 default = False, 
 help = "cache casu files")

parser.add_option("--stagepath", dest = "stagepath", 
 help = "Data staging area  path [Default outpath]")

parser.add_option("--logpath", dest = "logpath", 
 default = './',
 help = "Output path")

parser.add_option("--skip", dest = "nskip_files", 
 default = '0',
 help = "Number of files to skip")

parser.add_option("--nfiles", dest = "nfiles", 
 default = '1',
 help = "Number of files to process")

parser.add_option("--debug", action = "store_true", dest = "debug", 
 default = False, 
 help = "print status messages to stdout")

parser.add_option("--selftest", action = "store_true", dest = "selftest", 
 default = False, 
 help = "run listphot using imcore catalogue from image")

parser.add_option("--pawprints", action = "store_true", dest = "pawprints", 
 default = False, 
 help = "process constituent pawprints")

parser.add_option("--getcat", dest = "getcat", 
 default = False, 
 help = "copy catalogue from archive location")

parser.add_option("--search", action = "store_true", dest = "search", 
 default = False, 
 help = "Search catalogue")

parser.add_option("--radius", dest = "search_radius", 
 default = '4',
 help = "Catalogue search radius (arc seconds)")

parser.add_option("--vhs", action = "store_true", dest = "vhs", 
 default = False, 
 help = "run using vhs parameters")

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


logpath=options.logpath


debug=options.debug

logfile = 'casu_listdriven_batch.logfile'


make_sure_path_exists(logpath)
flog = open(logpath+"/Log_casu_listdriven_batch_%s.txt" % time_str, 'wt')
flogerr = open(logpath+"/Logerr_casu_listdriven_batch_%s.txt" % time_str, 'wt')

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
logger(flog, logdata)

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

if options.stagepath:
  stagepath = options.stagepath + '/'
if not options.stagepath: stagepath=outpath

logdata='stagepath: %s ' % stagepath
logger(flog, logdata)

#key=raw_input("Enter any key to continue: ")

# count the number lines in the input file
logdata = 'Opening file: %s' % infile
logger(flog, logdata)

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
print(ipos, len(listpath))
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


# This is the main work loop
# loop through line by line; there is an off-by-one issue
step=1
start=step-1
if options.nskip_files > 0: start=int(options.nskip_files)
print("start, step: ", start, step, len(records))

end=len(records)
if options.nfiles > 0: end=start+int(options.nfiles)
print("start, end, step: ", start, end, step)

#records=records[start:end:step]

nfiles=len(records)
verbose=False
n_errors=0
print('Number of files to process: ', len(records))
for line in records:
    ifile = ifile + 1
    line_save = line
    line = line.split()
    files = [line[i] for i in [3] if line[i]<>'null']
    filename = files[0]

    print() 
    logdata= "Processing field: ", ifile, ' of ',nfiles
    logger(flog, logdata)

    logdata= "Processing: ", filename
    logger(flog, logdata)

    process_image(files=files, outpath=outpath, Search=options.search)

    usage = resource.getrusage(resource.RUSAGE_SELF)
    for name, desc in [
      ('ru_utime', 'User time'),
      ('ru_stime', 'System time'),
      ]:
      print('%-25s (%-10s) = %s' % (desc, name, getattr(usage, name)))

    if debug: key=raw_input("Enter any key to continue: ")

if debug: print('debug on')

elapsedtime = time.time() - starttime
logdata = "Elapsed time: %.3f seconds" % (time.time() - starttime)
logger(flog, logdata)

debug=False

logdata = "Query completed"
flog.close()
flogerr.close()
