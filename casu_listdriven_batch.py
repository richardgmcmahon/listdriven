# $Source: /home/rgm/soft/ukidss/listdriven/RCS/listdriven_bulk.py,v $
# $Id: listdriven_bulk.py,v 1.1 2010/10/21 14:28:26 rgm Exp rgm $
"""

  WISHLIST:

   check the lockfile approach for race conditions and cases
   where program aborts; investigate lockfile methodology

   could look at option to cache the image files on /tmp if this
   would give speed improvement; call stagepath

   Need to trap the pyfits errors:
   except IOError as (errno, strerror):
   print "I/O error({0}): {1}".format(errno, strerror)

   could split mef pawprints into single extensions to speed up by
   reducing unecessary processing.

  FIX/trap some of the crash conditions; maybe email alert option

  Also, need to log errors to a file and print at the end
  an error repot summary

  verbose levels; turn off/reduce the ssh debug levels

  could add time to subprocess call to get some profiling info

  subprocess.call; subprocess.Popen; os.popen

  Both subprocess.call and subprocess.Popen are used for
  historical reasons due to previous use of os.popen

  subprocess.Popen is the main class with a number of wrapper
  convenience fucntions e.g. call, check_call, check_output

  see google: python subprocess call popen

  See also:

http://www.pythonforbeginners.com/os/subprocess-for-system-administrators

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


def get_remotefile(filename=None):
    # scp the science image
    # build the scp command; note spaces between each parameter
    # use -v to debug the scp/ssh
    #cmd ='scp -v -i ~/.ssh/id_dsa_nk ' \
    command ='time scp  -i ~/.ssh/id_dsa_nk ' \
      + ' ' + host + filename \
      + ' ' + stagepath +'/. '

    logdata=command
    logger(flog, logdata)

    itry=0
    Transfered=False

    trace = traceback.extract_stack()[-1]
    print(os.path.basename(trace[0]), ' line :', str(trace[1]))

    while (itry < iretry_max) and not Transfered:
      itry=itry+1
      if debug or verbose:
        print(command)
      #result=os.popen(command)
      #help(result)
      result=subprocess.Popen(command,
       stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
      print('Using subprocess.Popen and communicate')
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



def imcore_list_run(args=None, imagefile=None,
  confmapfile=None, listfile=None, outfile=None):
  """
  runs imcore_list using subprocess




  """

  # set the imcorelist parameters; maybe need to check these againsts
  # values in the header of the catalogue file
  nustep=-1
  print('args.vhs: ', args.vhs)
  print('args.cache: ', args.cache)

  if not args.vhs:
    nustep=pyfits.getval(imagefile,'NUSTEP',0)
    print('nustep = ', nustep)

  if args.vhs:
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
  #command = 'time nice -n19 ' + IMCORE_LIST  \
  command = 'time nice -n19 /home/rgm/bin/imcore_list ' \
        + ' ' + imagefile \
        + ' ' + confmapfile \
        + ' ' + listfile \
        + ' ' + outfile \
        + ' ' + threshold  \
        + ' --nbsize=' + nbsize \
        + ' --rcore=' + rcore \
        + ' --cattype=6 '

      # + ' --verbose '

  # save stdout and stderr to a logfile
  stdoutlog = open(logpath+'Logfile_stdout', 'w+')
  stderrlog = open(logpath+'Logfile_stderr', 'w+')

  logdata=command
  logger(flog, logdata)

  # convert comamd to string list for subprocess
  command = shlex.split(command)
  print('Print command as args: ')
  print(command)

  result = subprocess.call(command, \
       stderr=stderrlog, stdout=stdoutlog)

  Popen=False
  if Popen:
    result = subprocess.Popen(command, shell=True,
      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print('Using subprocess.Popen and communicate')
    output, errors = result.communicate()
    result=subprocess.check_output(command)
    #help(result)
    print('subprocesss output: ', output)
    print('subprocess errors: ', errors)

  #logdata='subprocess error status: ' + str(result)
  #logger(flog, logdata)
  trace = traceback.extract_stack()[-1]
  print(os.path.basename(trace[0]), ' line :', str(trace[1]))
  #if result is not 0: print(len(result))
  #print(result)

  result=0
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

    trace = traceback.extract_stack()[-1]
    print(os.path.basename(trace[0]), ' line :', str(trace[1]))

    #result=os.popen(command)
    result=subprocess.Popen(command)

    # create confidence map filename from the image filename
    confname = filename[:-4]+'_conf.fit'

    command ='scp -i ~/.ssh/id_dsa_nk ' \
      + ' '+ host + confname \
      + ' ' + stagepath +'/. '

    logdata=command
    logger(flog, logdata)

    trace = traceback.extract_stack()[-1]
    print(os.path.basename(trace[0]), ' line :', str(trace[1]))

    #result=os.popen(command)
    result=subprocess.Popen(command)

    if SelfTest:
      catfile = files[0][0:-4] + '_cat.fits'

      print('SelfTest catalogue file: ' + catfile)
      command ='time scp -i ~/.ssh/id_dsa_nk ' \
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


def process_image(filename=None, files=None, outpath=None):
    """
    refactoring in progress

    """

    deltatime=time.time()

    logdata= 'Image filename: ' + filename
    logger(flog, logdata)

    catfile = filename[0:-4] + '_cat.fits'
    logdata='Catalogue file: %s' % catfile
    logger(flog, logdata)

    confname = filename[:-4]+'_conf.fit'
    if verbose:
      logdata='Confidence map: %s' % confname
      logger(flog, logdata)


    # this is duplicated below!
    # find first file; needs to be refactored since only one file is input
    #if SelfTest:
    #  ipos1=catfile.rfind('/')
    #  ipos2=catfile.rfind('.fits')
    #  outfile=outpath+catfile[ipos1+1:ipos2] + '_listdriven.fits'
    #  if os.path.exists(outfile):
    #    logdata='List-driven data already exists for %s.' % (outfile)
    #    logger(flog, logdata)
    #    return


    if not SelfTest:
      # strip of the path from the filename
      ipos=filename.rfind('/')
      listfile = filename[ipos+1:] + '.radec'

      outfile = outpath + listfile + '_listdriven.fits'
      logdata ='Output file: %s' % outfile
      logger(flog, logdata)

      listfile = listpath + listfile
      logdata='Listfile: %s' % listfile


    if SelfTest:
      logdata= 'Running self test regression using imcore catalogue'
      logger(flog, logdata)
      catfile = filename[0:-4] + '_cat.fits'
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

    #another method
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

    get_remotefile(filename=filename)

    get_remotefile(filename=confname)

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

    if args.pawprints:
      # read the file to determine the constituent pawprints
      get_vista_pawprints(imagefile=imagefile, filename=filename,
       stagepath=stagepath)

    if SelfTest:
      print('SelfTest/Search catalogue file: ' + catfile)
      if os.path.exists(catfile):
        logdata='catalogue already exists %s.' % (catfile)
        logger(flog, logdata)

      if not os.path.exists(catfile):
        command ='time scp -i ~/.ssh/id_dsa_nk ' \
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
        print(os.path.basename(trace[0]), ' line:', str(trace[1]))
        key=raw_input("Debug: Enter any key to continue: ")

      imcore_list_run(args=args, imagefile=imagefile,
        confmapfile=confmapfile, listfile=listfile, outfile=outfile)

      print('listdriven photometry completed')
      #key=raw_input("Enter any key to continue: ")

    print('args.cache: ', args.cache)
    if os.path.exists(imagefile) and not args.cache:
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

    if os.path.exists(confmapfile) and not args.cache:
      print('Remove the confidence map:' + confmapfile)
      try:
        os.remove(confmapfile)
      except:
        logdata = "error removing confmapfile %s " % confmapfile
        logger(flogerr, logdata)
        pass

    if SelfTest:
      if os.path.exists(catfile)  and not args.cache:
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

def get_file(host=None, infile=None, transport='scp'):
    """
    scp file from remote host
    this could also be a http etc
    """

    # build the scp command; note spaces between each parameter
    # use -v to debug the scp/ssh
    #cmd ='scp -v -i ~/.ssh/id_dsa_nk '
    scp_verbose='-v'
    scp_verbose=''

    # ssh key location could be in cfg
    command ='time scp ' + scp_verbose + ' -i ~/.ssh/id_dsa_nk ' \
      + ' '+ host + infile \
      + ' ' + stagepath +'/. '

    logdata=command
    logger(flog, logdata)

    itry=0
    Transfered=False

    trace = traceback.extract_stack()[-1]
    print(os.path.basename(trace[0]), ' line :', str(trace[1]))

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
      ipos=infile.rfind('/')
      outfile = stagepath + infile[ipos+1:]

      if os.access(outfile, os.F_OK):
        Transfered=True

      if not os.access(outfile, os.F_OK):
        scpfailurefile=outfile + '.scpfailure'
        scpfailurefileh = open(scpfailurefile, 'wt')

        logdata='WARNING: file NOT transfered: %s' % infile
        logger(flogerr, logdata)
        Transfered=False
        delay=delaytime*(itry*itry)
        logdata='WAITING: %s seconds' % delay
        logger(flog, logdata)
        time.sleep(delay)
        continue

def search_catalogue(filename=None, listfile=None, outpath=None, radius=2.0):
    """

    """

    print('listfile: ', listfile)
    table=Table.read(listfile, format='ascii')
    table.pprint()
    ralist=table['col1']
    declist=table['col2']

    catfile = filename[0:-4] + '_cat.fits'
    logdata='Process catalogue file: %s' % (catfile)
    logger(flog, logdata)

    deltatime=time.time()

    ipos=filename.rfind('/')
    outfile=outpath + filename[ipos+1:-4] + '_search.fits'
    logdata='Output file: %s' % outfile
    logger(flog, logdata)
    if os.path.exists(outfile):
      print('Skipping since outfile already exists for %s.' % (outfile))
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

    # write pid and hostname into lockfile
    lkfile.write(strftime("%Y-%m-%dT%H-%M-%S", gmtime()))
    lkfile.write(':pid: '+ str(pid) + '\n')
    lkfile.flush()

    #lockfile=listfile + '.lock.' + str(pid)
    #logdata = "Lockfile: " + lockfile
    #logger(flog, logdata)
    # use flock module
    #lock = flock(lockfile, True).acquire()

    if os.access(listfile, os.F_OK):
      # Read list file skipping comment lines
      records = [item for item in open(listfile) if item[0]<>'#']
      numSources = len(records)
      logdata = 'Number of sources in listfile: %d' % numSources
      logger(flog, logdata)

    if not os.access(listfile, os.F_OK):
      logdata='List file: %s' % listfile
      logger(flog, logdata)
      logdata='List file problem'
      logger(flog, logdata)
      key=raw_input("Enter any key to continue: ")

    get_file(host=host, infile=catfile, transport='scp')

    logdata = "Delta elapsed time %.3f seconds" % (time.time() - deltatime)
    logger(flog, logdata)

    if not DryRun:
      logdata= 'Start processing the data'
      logger(flog, logdata)

      if debug:
        trace = traceback.extract_stack()[-1]
        print(os.path.basename(trace[0]), ' line :', str(trace[1]))
        key=raw_input("Debug: Enter any key to continue: ")

      ipos=filename.rfind('/')
      catfile = stagepath + catfile[ipos+1:]
      result=srlib.cat_cal(catfile, ralist, declist, 1, radius = radius)
      if result != None: result.write(outfile)

      print('catalogue search completed')
      #key=raw_input("Enter any key to continue: ")

    print('args.cache: ',args.cache)
    if os.path.exists(catfile) and not args.cache:
      print('Deleting data files used')
      print('Remove the cat file:' + catfile)
      try:
        os.remove(catfile)
      except OSError as (errno, strerror):
        logdata ="OS error({0}): {1}".format(errno, strerror)
        logger(flogerr, logdata)
        logdata = "error removing  catfile %s " % catfile
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

from astropy.table import Table

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
import cat_cal as srlib

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

# read config file; need to decide capatilisation rules
# for variable names
config.read('casu_listdriven_batch.cfg')
host = config.get('casu','host')
user = config.get('casu','user')
IMCORE_LIST = config.get('casu','imcore_list')
host = user + '@' + host + ':'

from time import strftime, gmtime, sleep

pid=os.getpid()
hostname=socket.gethostname()

time_str = strftime("%Y-%m-%dT%H-%M-%S", gmtime())

import argparse
parser = argparse.ArgumentParser(description='Listdriven photometry worker script; still under development')
"""
Input file format

# Sun Jun 22 15:00:56 2014
# /home/rgm/soft/idl/librgm/survey_radec_field.pro 1065
# /home/rgm/Projects/PS1/20140612/ps1zdrops_vhsarea1.fits
# /data/vhs/dqc/2014/vistaqc_20140430_tiles_vhs.fits
# J
       1      13      13  /data/apm29_a/vista/processed/20091103/v20091103_00093_st_tl.fit
       2      14      27  /data/apm29_a/vista/processed/20091103/v20091103_00129_st_tl.fit
       3       6      33  /data/apm29_a/vista/processed/20091103/v20091103_00165_st_tl.fit
"""

parser.add_argument("-i", "--input", dest = "infile",
 help = "Input filename")

parser.add_argument("--listpath", dest = "listpath",
 default = './',
 help = "Field list file path")

parser.add_argument("-o", "--outpath", dest = "outpath",
 default = '.',
 help = "Output path")

parser.add_argument("--cache", action = "store_true", dest = "cache",
 default = False,
 help = "cache processed image files in stagepath or outpath")

parser.add_argument("--stagepath", dest = "stagepath",
 help = "Data staging area  path [Default outpath]")

parser.add_argument("--logpath", dest = "logpath",
 default = './',
 help = "logpath")

parser.add_argument("--skip", dest = "nskip_files",
 default = '0',
 help = "Number of list files to skip")

parser.add_argument("--nfiles", dest = "nfiles",
 default='0',
 help = "Number of list files to process")

parser.add_argument("--debug", action = "store_true", dest = "debug",
 default = False,
 help = "print status messages to stdout")

parser.add_argument("--selftest", action = "store_true", dest = "selftest",
 default = False,
 help = "run listphot using imcore catalogue from image")

parser.add_argument("--dryrun", dest = "dryrun",
 default = False,
 help = "run Dry Run mode that test the input file handling")

parser.add_argument("--pawprints", action = "store_true", dest = "pawprints",
 default = False,
 help = "process constituent pawprints")

parser.add_argument("--getcat", dest = "getcat",
 default = False,
 help = "copy catalogue from archive location")

parser.add_argument("--search", dest = "search",
 default = False,
 help = "Search catalogue")

parser.add_argument("--radius", dest = "search_radius",
 default = '4',
 help = "Catalogue search radius (arc seconds)")

parser.add_argument("--vhs", dest = "vhs",
 default = False,
 help = "run using vhs parameters")

parser.add_argument("-v", "--verbose",
 action = "store_true", dest = "verbose",
 default = False,
 help = "print status messages to stdout")

parser.add_argument("-w", "--wait",
 dest = "wait", default = "0",
 help = "wait time between queries in seconds (0 = wait until query finishes) [0]")

parser.add_argument("--clobber", dest = "clobber",
 help = "Clobber or overwrite existing results")

#parser.add_argument("-l", "--logfile", dest = "logfile",
# help = "Logfile to log progress and summary")
# (options, args) = parser.parse_args()

args = parser.parse_args()
if args.verbose:
    print("verbose turned on")

logpath=args.logpath

debug=args.debug

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

logdata = 'args: '
logger(flog, logdata)

logdata = args
logger(flog, logdata)

logdata = 'args: '
logger(flog, logdata)

logdata = args
logger(flog, logdata)

infile = args.infile
logdata='infile: %s ' % infile
logger(flog, logdata)

outpath = args.outpath + '/'
logdata='outpath: %s ' % outpath
logger(flog, logdata)

if args.stagepath:
  stagepath = args.stagepath + '/'
if not args.stagepath: stagepath=outpath

logdata='stagepath: %s ' % stagepath
logger(flog, logdata)

#key=raw_input("Enter any key to continue: ")

# count the number lines in the input file
logdata = 'Opening file: %s' % infile
logger(flog, logdata)

try:
  numLines = len(open(infile).readlines())
except:
  traceback.print_exc(file=sys.stdout)
  print('input file cannot be opened; exiting.... ', None)
  sys.exit(1)

logdata = 'Input file has length: %d' % numLines
logger(flog, logdata)


# Read field list file skipping comment lines
records = [item for item in open(infile) if item[0]<>'#']
numFields = len(records)

logdata = 'Total number of fields: %d' % numFields
logger(flog, logdata)

starttime=time.time()

SelfTest = args.selftest

listpath=args.listpath
ipos=listpath.rfind('/')
print(ipos, len(listpath))
if ipos+1 < len: listpath=args.listpath + '/'
logdata='listpath: %s' % listpath
logger(flog, logdata)

ifile=0
DryRun = args.dryrun
if DryRun: print('Running in Dry Run mode')
nfiles=len(records)
#help(records)
#key=raw_input("Enter any key to continue: ")

logdata = "Total elapsed time %.3f seconds" % (time.time() - starttime)
logger(flog, logdata)


# This is the main work loop
# loop through line by line; there is an off-by-one issue
step=1
start=step-1
if args.nskip_files > 0: start=int(args.nskip_files)
print("start, step: ", start, step, len(records))

end=len(records)
if args.nfiles > 0: end=start+int(args.nfiles)
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

    process_image(filename=filename, files=files, outpath=outpath)

    if args.search:
      # strip of the path from the filename
      ipos=filename.rfind('/')
      listfile = listpath + filename[ipos+1:] + '.radec'
      logdata='Process search with list file: %s' % (listfile)
      search_catalogue(filename=filename, listfile=listfile,
        outpath=outpath, radius=4.0)

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
