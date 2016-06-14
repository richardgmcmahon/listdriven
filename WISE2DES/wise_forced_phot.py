"""

multiprocessing version based on

(i) multisite website crawler example:
http://toastdriven.com/blog/2008/nov/11/brief-introduction-multiprocessing/

(ii) listdriven_bulk.py written for SDSS BOSS

WARNINGs:

Be careful where you write the log file to:

Check for stagnant lockfile and also maybe completed files

Problems:

Could have a Deadlock
It was not caused by logging;
It occurs with files are being skipped; not sure why


TODO:
(-1) limit the wavebands e.g. we just need g and i for moment for Lensed quasars
     but not worth the effort restricting it maybe?
(0)    add command line input for nworkers, sleeptime, log location, overwrite
(i)    need to add DES data staging
(ii)   could include running imcore on unWISE-R images
(iii)  VISTA
(iv)   VST
(v)    SDSS etc
(vi)   Add into LSST framework
(viii) Add dry run to check all input files exist
(ix)   Run the WISE list part only to check the WSDB part
(x)    Continue option

It takes about 4 minutes to process 5 bands for a single tile.
The i/o is 100MB per file which is 1 second accross the
network i.e. 5seconds in 200.
Therefore if we can in principle get 40 threads running without
saturating the disk i/o.

i.e. 50 tiles in 4 minutes; or 500-1000 tiles per hour
=> 300deg^2 per hour or all DES in less than 24 hours
over 40 cores; i.e. 1000 core hours

4 mins per tile; 40,000 mins gives same.

calx024 tests: 12 tiles in 12mins

HISTORY

20160411: rgm- parallelising
20160411: rgm- made into function
20160411: rgm- making pep8 and future compliant
20160323: slr- original version

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

# standard libs
import os
import sys
import time

import logging
import traceback
import configparser
import warnings

import subprocess

# 3rd party libs
import numpy as np

import astropy.io.fits as fits
from astropy import wcs
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning


# local libs
sys.path.append('/home/rgm/soft/python/lib/')

import sqlutil

from librgm.rd_config_wsdb import rd_config_wsdb

from parse_config import *

warnings.filterwarnings('ignore', category=AstropyWarning, append=True)

def create_logger(logdir=None, loglevel=logging.INFO):
    """
    Need to add traceback to the error log output via a function

    TODO:
    (1) add a lockfile like listdriven
    (2) Quieter skipping

    """

    logfile = '/tmp/logfile.log'
    if logdir is not None:
        logfile = logdir + '/' + logfile

    # handler = logging.handlers.RotatingFileHandler(
    #    logfile,
    #    backupCount=5
    # )

    # logger = logging.getLogger("Rotating Log")
    # handler = logging.handlers.TimedRotatingFileHandler(logfile,
    #                                                    when="h",
    #                                                    interval=1,
    #                                                    backupCount=24)

    format = ('%(levelname)s:'
              '%(asctime)s:'
              # '%(threadName)s:'
              # '%(thread)d:'
              '%(processName)s:'
              '%(process)d:'
              '%(filename)s:'
              '%(funcName)s():'
              '%(lineno)d:'
              '%(message)s')

    datefmt='%Y-%m-%dT%H:%M:%S'
    # DEBUG: Detailed information, typically only for diagnosing problems.
    # INFO:  A message confirming that things are working as expected
    # WARNING: An indication that something unexpected happened
    # ERROR: indicates a more serious issue, including exceptions

    # level = logging.DEBUG
    # level = level
    # level = logging.WARNING
    # level = logging.ERROR
    # level = logging.CRITICAL
    # Note basicConfig can only be called once
    logging.basicConfig(level=loglevel,
                        filename=logfile,
                        format=format,
                        datefmt=datefmt)

    # logging.getLogger('').addHandler(handler)
    # logger.addHandler(handler)
    error_logfile = 'logfile.error'
    if logdir is not None:
        error_logfile = logdir + '/' + error_logfile

    # logging.basicConfig(level=logging.ERROR,
    #                    filename=error_logfile,
    #                    format='%(levelname)s:'
    #                           '%(asctime)s:'
    #                           '%(filename)s:'
    #                           '%(funcName)s():'
    #                           '%(lineno)d:'
    #                           '%(message)s',
    #                    datefmt='%Y-%m-%dT%H:%M:%S')

    # create another logger from console
    console = logging.StreamHandler()
    logging.getLogger('').addHandler(console)
    console.setLevel(loglevel)
    # console.setLevel(logging.ERROR)

    # add another formatted record
    # formatter = logging.Formatter(
    #    '%(levelname)s - %(asctime)s - %(name)s - %(message)s'
    # )
    #

    formatter = logging.Formatter(format, datefmt)
    console.setFormatter(formatter)

    # logging.basicConfig(level=logging.INFO)
    # logger = logging.getLogger(__name__)
    # formatter = logging.Formatter(
    #    '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    # logger.setFormatter(formatter)

    # console = logging.StreamHandler()
    # console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    # formatter = logging.Formatter('%(levelname)-8s:%(asctime)s: %(name)-12s: %(message)s', "%Y-%m-%dT%H:%M:%S")

    # formatter = logging.Formatter('%(levelname)-8s:%(asctime)s: %(name)-12s: %(message)s')

    # tell the handler to use this format
    # console.setFormatter(formatter)
    # add the handler to the root logger
    # logging.getLogger('').addHandler(console)

    # formatter = logging.Formatter(
    #    '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    # )

    return


def get_cats_status(tilename=None,
                    outdir=None,
                    config_file=None):
    """Determine if outfile already exists


    returns 1 if data already exists and thus next stage can be skipped

    """

    import os.path

    cats_status = 0

    # get the default logger
    logger = logging.getLogger()

    logger.info('Processing Tile: %s', tilename)

    logger.debug('Reading config file: %s', config_file)
    config = configparser.RawConfigParser()
    config.read(config_file)

    # could check these exist here and stop if missing
    IMCORE_LIST = config.get("casu", "imcore_list")
    CLASSIFY = config.get("casu", "classify")

    datpath = config.get("des", "datapath")
    # should add trailing slash if not present

    OUTDIR = config.get("des", "outpath")

    # this is a bit of a mess
    outpath = OUTDIR + '/' + tilename + '/'

    logger.debug('Validating outpath: %s', outpath)
    if os.path.exists(outpath):
        logger.debug('Outpath exists: %s', outpath)

    if not os.path.exists(outpath):
        logger.debug('Creating: %s', outpath)
        try:
            os.makedirs(outpath)
        except Exception as error:
            logger.exception(error)
            print("Unexpected error:", sys.exc_info()[0])
            traceback.print_exc(file=sys.stdout)
            raise

    logger.debug('Outpath exists: %s', outpath)

    coord_file = outpath + '/' + tilename + '_' + "test_input.txt"
    # could skip getting the coord file
    # currently skipped within mk_cats if it already exists
    # os.path.exists(coord_file)

    # check if output file exists for files that have images
    # need to get the missing images e.g. get_images function

    logger.debug('Processing wavebands')
    wavebands = ["g", "r", "i", "z", "Y"]
    for band in wavebands:

        logger.debug("Processing band: %s", band)
        im_file = datapath + tilename + "/" + tilename + \
             "_" + band + ".fits.fz"
        logger.debug('Processing file: %s %s', im_file, band)

        if os.path.isfile(im_file):
            outfile = outpath + tilename + "_WISEfp_" + band + ".fits"

            if os.path.isfile(outfile):
                cats_status = cats_status + 1

    # could give status per waveband via a dict
    logger.debug('cats_status: %s', cats_status)
    return cats_status


def get_coords(infile=None, coord_file=None, tilename=None):
    """
    # Take DES tile, find all WISE objects on that tile.
    # infile = "/data/desardata/Y1A1/" + tile + "/" + tile + "_z.fits.fz"
    # DATAPATH = '/data/desardata/Y1A1/'

    """
    import math

    logger = logging.getLogger()

    outpath = os.path.dirname(coord_file)
    if not os.path.exists(outpath):
        try:
            os.makedirs(outpath)
        except Exception as error:
            logger.exception(error)
            print("Unexpected error:", sys.exc_info()[0])
            traceback.print_exc(file=sys.stdout)
            raise

    logger.debug('Reading tile image to get corners: %s', infile)
    with fits.open(infile) as fhlist:
        hdr = fhlist[1].header
        w = wcs.WCS(hdr, naxis=2)

    # If this is ever used for non square inmages this needs checking
    pix_corners = [[0, 0], [0, hdr["NAXIS1"]], [
        hdr["NAXIS2"], 0], [hdr["NAXIS1"], hdr["NAXIS2"]]]
    w_corners = w.wcs_pix2world(pix_corners, 1)

    [corner_ra, corner_dec] = zip(*w_corners)

    """
    corner_ra = sorted(corner_ra)
    corner_dec = sorted(corner_dec)

    #Take the slightly smaller square inside the not square image
    ra_min = corner_ra[1]
    ra_max = corner_ra[2]
    dec_min = corner_dec[1]
    dec_max = corner_dec[2]
    """

    delta_ra = max(corner_ra) - min(corner_ra)
    delta_dec = max(corner_dec) - min(corner_dec)

    dec_mean = (min(corner_dec) + min(corner_dec)) / 2.0
    delta_ra = delta_ra * math.cos(math.radians(dec_mean))

    logger.debug('dec_mean: %s', dec_mean)
    logger.debug('delta_ra, delta_dec: %s; %s', delta_ra, delta_dec)

    # raw_input("Enter any key to continue: ")
    # Some of these objects won't be on the DES tile so will return infs
    ra_min = str(min(corner_ra))
    ra_max = str(max(corner_ra))
    dec_min = str(min(corner_dec))
    dec_max = str(max(corner_dec))

    query_wise = "select ra, dec FROM allwise.main WHERE ra > " + \
        ra_min + \
        " and ra < " + ra_max + " and dec > " + \
        dec_min + " and dec < " + dec_max

    # special case where crossing  24hrs ra > ra_max or ra < ra_min
    if delta_ra > 180.0:

        query_wise = "select ra, dec FROM allwise.main WHERE (ra > " + \
            ra_max + \
            " or ra < " + ra_min + ") and dec > " + \
            dec_min + " and dec < " + dec_max

    db, host, user, password, db_table = rd_config_wsdb(table='wise')

    RA_wise, DEC_wise = sqlutil.get(
        query_wise, db=db, host=host,
        user=user, password=password
    )
    logger.info('%s: Number of wise sources: %s', tilename, len(RA_wise))

    info = ""
    for (ra, dec) in zip(RA_wise, DEC_wise):
        info += (str(ra) + " " + str(dec) + "\n")

    with open(coord_file, "w") as f:
        logger.debug('Writing coords: %s', coord_file)
        try:
            f.write(info)
            logger.info('Successfully written coordinate file: %s', coord_file)
        except IOError as e:
            print("I/O error({0}): {1}".format(e.errno, e.strerror))
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise

        # except IOError as (errno, strerror):
        #    print("I/O error({0}): {1}".format(errno, strerror))
        #    print('Failed to write coordinate file: ' + str(e))
        #    logger.error('Failed to write coordinate file: ' + str(e))
        #    sys.exit()
        #    pass

        # except OSError as (errno, strerror):
        #    print("Hello World:OS error({0}): {1}".format(errno, strerror))
        #    pass

        # except:
        #    print("Unexpected error:", sys.exc_info()[0])
        #    # traceback.print_exc(file=sys.stdout)
        #    raise

    return coord_file


def mk_cats(tile, rcore, outdir, config_file, overwrite=False):
    """

    """

    import math

    # get the default logger
    logger = logging.getLogger()
    logger.info('Processing Tile: %s', tile)

    logger.debug('Reading config file: %s', config_file)

    config = configparser.RawConfigParser()
    config.read(config_file)

    IMCORE_LIST = config.get("casu", "imcore_list")
    CLASSIFY = config.get("casu", "classify")

    datapath = config.get("des", "datapath")
    OUTPATH = config.get("des", "outpath")
    outdir = OUTPATH
    outpath = outdir + '/' + tile + '/'

    rcore = config.get("des", "rcore")

    logger.debug('Checking outpath: %s', outpath)

    if not os.path.exists(outpath):
        logger.debug('Creating: %s', outpath)
        try:
            os.makedirs(outpath)
        except Exception as e:
            print("Unexpected error:", sys.exc_info()[0])
            traceback.print_exc(file=sys.stdout)
            raise

    coord_file = outpath + '/' + tile + '_' + "test_input.txt"

    logger.debug('Checking coord_file: %s', coord_file)
    if os.path.isfile(coord_file):
         logger.debug('coord_file: %s  Exists', coord_file)
    # could skip getting the coord file
    if not os.path.isfile(coord_file):
        logger.debug('coord_file: %s  Does not exist', coord_file)
        infile = datapath + tile + "/" + tile + "_z.fits.fz"

        if os.path.isfile(infile):
            logger.debug('footprint file : %s  exists', infile)

        if not os.path.isfile(infile):
            logger.error('File does not exist: %s', infile)
            logger.error('Exiting....')
            #sys.exit()
            logger.error('Should have exited....')
            return 'Hello_World'

        logger.debug('Get the coordinate list for: %s', infile)

        get_coords(infile=infile, coord_file=coord_file, tilename=tile)

    # Now process the image data using the list of
    for band in ["g", "r", "i", "z", "Y"]:

        logger.debug('Processing band: %s rcore: %s', band, str(rcore))

        im_file = datapath + tile + "/" + tile + "_" + band + ".fits.fz"
        logger.info('Processing file: %s %s', im_file, band)

        outfile = outpath + tile + "_WISEfp_" + band + ".fits"
        logger.info("Output file: %s", outfile)

        if not os.path.exists(outpath):
            logger.debug('Creating: %s', outpath)
            try:
                os.makedirs(outpath)
            except Exception as e:
                print("Unexpected error:", sys.exc_info()[0])
                traceback.print_exc(file=sys.stdout)
                raise

        if not os.path.exists(im_file):
            logger.warning('Input image file does not exist: %s', im_file)

        if not os.path.exists(coord_file):
            logger.warning('Input coordinate file does not exist: %s', coord_file)

        if not os.path.exists(IMCORE_LIST):
            print('File does not exist:', IMCORE_LIST)

        if os.path.exists(outfile) and not overwrite:
            print('Output file already exists:', outfile)
            print('Skipping')
            # continue

        if not os.path.exists(outfile) or overwrite:

            logger.info('Starting IMCORE_LIST on %s: ', im_file)
            subprocess.call([IMCORE_LIST, im_file, "noconf",
                            coord_file, outfile,
                            "1.5", "--rcore=" + rcore,
                            "--cattype=6", "--noell"])

            logger.info('Starting CLASSIFY on %s: ', outfile)
            subprocess.call([CLASSIFY, outfile])

            logger.info('Completed processing: %s', im_file)

    return outpath

def calibrate(tile, outdir, config_file=None, overwrite=False):

    # Calibrate the cats, just doing aper 3 for now
    # Negative flux -> negative magnitude

    from astropy import wcs
    from astropy.table import Table
    import numpy as np

    logger = logging.getLogger()

    config = configparser.RawConfigParser()
    config.read(config_file)

    release = config.get("des", "release")

    for band in ["g", "r", "i", "z", "Y"]:
        """
        im_file = "/data/desardata/Y1A1/" + tile + "/" +
            tile + "_" + band + ".fits.fz"
        with fits.open(im_file) as fhlist:
            des_hdr = fhlist[1].header
        zpt = des_hdr["SEXMGZPT"]
        expt = des_hdr["EXPTIME"]
        """

        fp_file = outdir + tile + "_WISEfp_" + band + ".fits"
        logger.info('Reading: %s', fp_file)
        t = Table.read(fp_file)

        # Remove Blanks
        for col in t.columns:
            if "Blank" in col:
                del t[col]

        # Add coords
        im_file = datapath + tile + "/" + tile + "_" + band + ".fits.fz"
        logger.info('Reading: %s', im_file)
        with fits.open(im_file) as fhlist:
            des_hdr = fhlist[1].header
        w = wcs.WCS(des_hdr, naxis=2)
        pix_coords = zip(t["X_coordinate"], t["Y_coordinate"])
        w_coords = w.wcs_pix2world(pix_coords, 1)
        w_coords = zip(*w_coords)
        t["RA_CALC"] = w_coords[0]
        t["DEC_CALC"] = w_coords[1]

        with fits.open(fp_file) as fhlist:
            fp_hdr = fhlist[1].header

        rcore = fp_hdr["RCORE"]
        skyn = fp_hdr["SKYNOISE"]
        zpt = fp_hdr["SEXMGZPT"]
        expt = fp_hdr["EXPTIME"]

        aperlist = ['1','2','3','3','4','5','6','7']

        for aper in aperlist:

            apcor = 0.0
            if aper <= '7':
                apcor = fp_hdr["APCOR" + aper]
            logger.debug('APCOR%s: %s', aper, apcor)

            try:
                apcor = fp_hdr["APCOR" + aper]
            except Exception as e:
                print("Unexpected error:", sys.exc_info()[0], fp_file)
                traceback.print_exc(file=sys.stdout)
                raise

            # extinction and airmass -
            #    in theory included in the DES zeropoints

            t["MAG_" + aper + "_CAL"] = [0.0] * len(t)
            t["MAG_ERR_" + aper + "_CAL"] = [0.0] * len(t)

            # deal with -ve flux; = 0 not covered
            ids = np.where((t["Aper_flux_" + aper] < 0.0))[0]

            # t["MAG_3_CAL"][ids] = \
            #   (-1 * (zpt - 2.5 * np.log10(np.fabs(t["Aper_flux_3"][ids])) +
            #   2.5 * np.log10(expt) - apcor))

            # log10 (abs(flux))
            t["MAG_" + aper + "_CAL"][ids] = \
               -1 * \
               (zpt - 2.5*np.log10(np.fabs(t["Aper_flux_" + aper][ids])) \
              - apcor)

            # add flux error so that magerr can be calculated
            fluxes = t["Aper_flux_" + aper][ids] + \
                t["Aper_flux_" + aper + "_err"][ids]

            # t["MAG_ERR_3_CAL"][ids] = \
            #    (zpt - 2.5 * np.log10(np.fabs(fluxes)) +
            #    2.5 * np.log10(expt) - apcor) - np.fabs(t["MAG_3_CAL"][ids])

            t["MAG_ERR_" + aper + "_CAL"][ids] = \
                (zpt - 2.5*np.log10(np.fabs(fluxes)) - apcor) \
                - np.fabs(t["MAG_" + aper + "_CAL"][ids])

            ids = np.where((t["Aper_flux_" + aper] > 0.0))[0]

            t["MAG_" + aper + "_CAL"][ids] = \
                (zpt - \
                 2.5 * np.log10(np.fabs(t["Aper_flux_" + aper][ids])) - apcor)

            # t["MAG_3_CAL"][ids] = \
            #    (zpt - 2.5 * np.log10(np.fabs(t["Aper_flux_3"][ids])) +
            #     2.5 * np.log10(expt) - apcor)

            fluxes = t["Aper_flux_" + aper][ids] + \
                t["Aper_flux_" + aper + "_err"][ids]

            # t["MAG_ERR_3_CAL"][ids] = \
            #    (zpt - 2.5 * np.log10(np.fabs(fluxes)) +
            #     2.5 * np.log10(expt) - apcor) - np.fabs(t["MAG_3_CAL"][ids])
            t["MAG_ERR_" + aper + "_CAL"][ids] = \
               (zpt - 2.5*np.log10(np.fabs(fluxes)) - apcor) \
               - np.fabs(t["MAG_" + aper + "_CAL"][ids])

        # compute the Kron magnitude

        t["MAG_KRON_CAL"] = [0.0] * len(t)
        t["MAG_ERR_KRON_CAL"] = [0.0] * len(t)

        ids = np.where((t['Kron_flux'] < 0.0))[0]

        t["MAG_KRON_CAL"][ids] = \
               -1 * \
               (zpt - 2.5*np.log10(np.fabs(t["Kron_flux"][ids])) \
              - apcor)

        # add flux error so that magerr can be calculated
        fluxes = t["Kron_flux"][ids] + \
                 t["Kron_flux_err"][ids]

        t["MAG_ERR_KRON_CAL"][ids] = \
                (zpt - 2.5*np.log10(np.fabs(fluxes)) - apcor) \
                - np.fabs(t["MAG_KRON_CAL"][ids])

        ids = np.where((t["Kron_flux"] > 0.0))[0]

        t["MAG_KRON_CAL"][ids] = \
                (zpt - \
                 2.5 * np.log10(np.fabs(t["Kron_flux"][ids])) - apcor)

        fluxes = t["Kron_flux"][ids] + \
                 t["Kron_flux_err"][ids]

        t["MAG_ERR_KRON_CAL"][ids] = \
               (zpt - 2.5*np.log10(np.fabs(fluxes)) - apcor) \
               - np.fabs(t["MAG_KRON_CAL"][ids])

        calfile = fp_file[:-5] + "_cal.fits"
        t.write(calfile, overwrite=overwrite)

        logger.info('Calibration complete: %s', calfile)


def join_cats(tile, outdir, overwrite=False):
    """

    """
    from astropy.table import Table, hstack

    logger = logging.getLogger()

    t_out = Table.read(outdir + tile + "_WISEfp_g_cal.fits")

    for col in t_out.columns:
        t_out.rename_column(col, col + "_G")

    for band in ["r", "i", "z", "Y"]:
        fp_file = outdir + tile + "_WISEfp_" + band + "_cal.fits"
        t = Table.read(fp_file)

        for col in t.columns:
            t.rename_column(col, col + "_" + band.upper())

        # stack onto end of the table
        t_out = hstack([t_out, t])

    outfile = outdir + tile + "_WISEfp.fits"
    # add a counter to the WISE list
    t_out["WISE_NUM"] = np.arange(1, len(t_out)+1)
    t_out.write(outfile, overwrite=overwrite)

    logger.info('Join complete: %s', outfile)


def add_DEScat(tile, outdir, filename_suffix='_WISEfp',
               overwrite=False):


    from astropy.table import Table, hstack
    import numpy as np
    import match_lists

    config = configparser.RawConfigParser()
    config.read(config_file)

    release = config.get("des", "release")

    infile = datapath + tile + "/" + tile + ".fits"

    logger.info('Read: %s', infile)

    t_cat = Table.read(infile)
    t_fp = Table.read(outdir + tile + "_WISEfp.fits")

    dists, inds = match_lists.match_lists(t_cat["ALPHAWIN_J2000_G"].data, t_cat["DELTAWIN_J2000_G"].data, t_fp["RA_CALC_G"].data, t_fp["DEC_CALC_G"].data, 5.0/3600.0, 1)

    ids = np.where( (inds <> len(t_fp)) )[0]
    t_cat = t_cat[ids]
    t_fp = t_fp[inds[ids]]
    t_out = hstack([t_fp, t_cat])
    t_out.write(outdir + tile + "_WISEfp_DEScat.fits",
        overwrite=overwrite)



def phot_check(tile, outdir):

    import matplotlib.pyplot as plt
    from astropy.table import Table, hstack
    import numpy as np

    t = Table.read(outdir + tile + "_WISEfp_DEScat.fits")

    for band in ["G", "R", "I", "Z", "Y"]:
        print('band: ', band)
        ids = np.where( (t["MAG_APER_4_" + band] < 90.0) & (t["MAG_3_CAL_" + band] > 0.0) )[0]
        plt.plot(t["MAG_APER_4_" + band][ids], t["MAG_3_CAL_" + band][ids]-t["MAG_APER_4_" + band][ids], "k.")
        plt.ylabel("Calculated Magnitude - Catalogue Aperture 4 Magnitude")
        plt.xlabel("Catalogue Aperture 4 Magnitude")
        plt.title(tile + " " + band)
        plt.axhline(0.0)
        plt.show()


        ids = np.where( (t["MAG_PSF_" + band] < 90.0) & (t["MAG_3_CAL_" + band] > 0.0) )[0]
        plt.plot(t["MAG_PSF_" + band][ids], t["MAG_3_CAL_" + band][ids]-t["MAG_PSF_" + band][ids], "k.")
        plt.ylabel("Calculated Magnitude - Catalogue PSF Magnitude")
        plt.xlabel("Catalogue PSF Magnitude")
        plt.title(tile + " " + band)
        plt.axhline(0.0)
        plt.show()

def WISE_match(tile, outdir, width_arcsecs = 5.0,
               checkplots=False, overwrite=False):

    import numpy as np
    from astropy.table import Table, hstack
    import srpylib

    infile = outdir + tile + "_WISEfp_DEScat.fits"
    print('Reading: ', infile)
    t = Table.read(outdir + tile + "_WISEfp_DEScat.fits")

    t = srpylib.WISE_match(t, "RA_CALC_G", "DEC_CALC_G",
         width = width_arcsecs/3600.0,
         c_graph = checkplots, DES_box = True)

    t.write(outdir + tile + "_WISEfp_DEScat_WISE_match.fits",
        overwrite=overwrite)

def nearest_neighbour(tile, outdir, checkplots=False):

    from astropy.table import Table, hstack
    import srpylib
    import match_lists
    import matplotlib.pyplot as plt
    import numpy as np

    db, host, user, password, db_table = rd_config_wsdb(table='wise')

    t = Table.read(outdir + tile + "_WISEfp_DEScat.fits")

    ra_min, ra_max, dec_min, dec_max = srpylib.min_max_match(t["RA_CALC_G"], t["DEC_CALC_G"], 1.0, DES_box = True)

    query_wise = "select ra, dec FROM allwise.main WHERE ra > " + ra_min + " and ra < " + ra_max + " and dec > " + dec_min + " and dec < " + dec_max
    RA, DEC = sqlutil.get(query_wise, db=db, host=host, user=user, password=password)

    dists, inds = match_lists.match_lists(t["RA_CALC_G"].data, t["DEC_CALC_G"].data, RA, DEC, 1.0, 2)

    match_distance = dists[:,1]*3600.0

    if checkplots:
        plt.axvline(6.1)
        plt.axvline(6.4)
        plt.axvline(6.5)
        plt.axvline(12.0)
        plt.hist(match_distance, bins = 100)
        plt.xlabel('Match Distance, "')
        plt.ylabel("Number")
        plt.title("Nearest Neighbours, " + tile)
        plt.show()

    print("Median Distance:", np.median(match_distance))


def wise_forced_phot(tilename=None, overwrite=False, dryrun=False,
                     config_file=None, start_point=None):
    """

    need to review the path munging and flow

    """

    #print('wise_forced_phot:', tilename)

    # get the default logger
    logger = logging.getLogger()
    logger.info('Starting Tile: %s', tilename)

    if start_point is not None:
        logger.info('Processing start point: %s', start_point)

    logger.debug('overwrite: %s', overwrite)
    logger.debug('dryrun: %s', dryrun)

    logger.debug('Reading config file: %s', config_file)
    config = configparser.RawConfigParser()
    config.read(config_file)

    datapath = config.get("des", "datapath")
    OUTDIR = config.get("des", "outpath")

    logger.info('OUTDIR: %s', OUTDIR)
    logger.info('tilename: %s', tilename)

    # this is a bit of a mess
    outpath = OUTDIR + '/' + tilename + '/'

    logger.info('outpath: %s', outpath)

    lockfile = outpath + tilename + '.lock'

    logger.info('lockfile: %s', lockfile)

    if os.path.exists(lockfile):
        logger.info('lockfile exists: ' + lockfile)

    cats_status = 0
    if not os.path.exists(lockfile):
        cats_status = get_cats_status(tilename=tilename,
                                      outdir=outdir,
                                      config_file=config_file)
        logger.info('cats_status: %s', cats_status)

    if cats_status != 0:
        logger.info('Catalogue already completed for Tile:%s', tilename)

    if (cats_status == 0 or
         overwrite is True) and not os.path.exists(lockfile):

        # create lockfile
        try:
            logger.info("Create lockfile: %s", lockfile)
            open(lockfile, 'wt')
        except Exception as e:
            print("Unexpected error:", sys.exc_info()[0])
            traceback.print_exc(file=sys.stdout)
            raise

        start_point = 'mk_cats'
        if start_point == 'mk_cats':
            logger.info('Make the list driven catalogues: %s', tilename)
            outpath = mk_cats(tilename, rcore, outdir, config_file,
                              overwrite=overwrite)

            logger.info('Calibrate Tile: %s', tilename)
            calibrate(tilename, outpath, config_file=config_file,
                      overwrite=overwrite)

            logger.info('Join Tile catalogues: %s', tilename)
            join_cats(tilename, outpath, overwrite=overwrite)


        logger.info('Pairwise match to DES catalogues: %s', tilename)
        add_DEScat(tilename, outpath, overwrite=overwrite)

        logger.info('Calibration check plots: %s', tilename)
        if checkplots: phot_check(tilename, outpath)

        logger.info('Pairwise match to WISE: %s', tilename)
        WISE_match(tilename, outpath, checkplots=checkplots,
                   overwrite=overwrite)

        logger.info('WISE pairwise self-neighbour plot: %s', tilename)
        nearest_neighbour(tilename, outpath, overwrite=overwrite)

        logger.info('Completed Tile:%s', tilename)

        if os.path.exists(lockfile):
            # remove lockfile
            try:
                logger.info("Delete lockfile: %s", lockfile)
                os.remove(lockfile)
            except Exception as e:
                print("Unexpected error:", sys.exc_info()[0])
                traceback.print_exc(file=sys.stdout)
                raise

        completedfile = outpath + tilename + '.completed'
        open(completedfile, 'a')

    # without this things can hang maybe; not sure it helps
    sleep_time = 1
    time.sleep(sleep_time)

    return


# This is the worker function
def worker_tile(work_queue, done_queue):
    """

    based on website crawler example def worker

    """

    logger = logging.getLogger()

    logger.info('module name: %s', __name__)

    # print process pid and and parent process pid: ppid
    if hasattr(os, 'getppid'):  # only available on Unix
        logger.info('Parent process: %s', os.getppid())
    logger.info("Worker process: {}".format(os.getpid()))

    try:
        # interate over work_queue until you reach 'STOP'
        iprocess = 0
        for item in iter(work_queue.get, 'STOP'):
            iprocess = iprocess + 1
            starttime = time.time()

            logger.info('%s  Processing: %s' % (str(iprocess), item))

            # logger.info(str(iprocess)+': Processing: %s', item)

            # I am not sure where tilename, overwrite and config_file come from
            # maybe use inspect to get the values
            status_code = wise_forced_phot(tilename=item,
                                           overwrite=overwrite,
                                           config_file=config_file)

            logger.info("Status Code: %s %s", status_code, item)

            logger.info(str(iprocess) + ': Completed: %s', item)

            done_queue.put(
                "%s - %s got %s." %
                (current_process().name, item, status_code)
            )

            name = multiprocessing.current_process().name
            logger.info("Process name: %s %s", name, item)
            logger.info("Worker process: {}".format(os.getpid()))

            # add a short sleep to avoid some race conditions
            sleep_time = 10
            sleep_time = 1
            if sleep_time > 0:
                logger.info('Sleeping...: %s %s %s', sleep_time, str(iprocess), item)
            time.sleep(sleep_time)
            logger.info('Time elapsed: %s %s' % (time.time() - starttime, item))

    except Exception, e:
        done_queue.put(
            "%s failed on %s with: %s" %
            (current_process().name, item, e.message)
        )

    return True


def mp_forcephot(tilelist=None,
                 overwrite=False,
                 dryrun=False,
                 nworkers=None,
                 loglevel=logging.INFO):
    """

    calls worker function: worker_tile

    """

    logger = logging.getLogger()

    ncores = multiprocessing.cpu_count()
    logger.info('Number of cores: %s', ncores)
    if nworkers is None:
        nworkers = ncores * 1

    logger.info('Number of workers: %s', nworkers)

    # create Queue in module multiprocessing.queues object
    logger.info('Create empty queues')
    work_queue = multiprocessing.Queue()
    done_queue = multiprocessing.Queue()

    # create list for the processes
    processes = []

    starttime = time.time()
    computetime = 0.0

    logger.info('Fill work_queue')
    index = 0
    for item in tilelist:
        index = index + 1
        logger.info('Putting item into queue %s :%s', str(index), item)
        work_queue.put(item)

    # raw_input("Enter any key to continue: ")
    logger.info('Process work_queue')
    for w in xrange(nworkers):
        multiprocessing.log_to_stderr(loglevel)

        # create progress to run work_tile
        p = Process(target=worker_tile, args=(work_queue, done_queue))

        logger.info('Worker: %s', w + 1)

        name = multiprocessing.current_process().name
        logger.info('Start process: %s', name)
        # raw_input("Enter any key to continue: ")
        p.start()
        logger.info('processes.append')
        # raw_input("Enter any key to continue: ")
        processes.append(p)
        work_queue.put('STOP')
        # raw_input("Enter any key to continue: ")

    # from http://stackoverflow.com/questions/20167735/python-multiprocessing-pipe-deadlock
    # for p in processes:
    for workerno, process in enumerate(processes):
        logger.info(
            'Attempting process.join() of workerno {0}'.format(workerno)
        )
        p.join()

    done_queue.put('STOP')

    logger.info('process results')
    for status in iter(done_queue.get, 'STOP'):
        logger.info('status: %s', status)

    logger.info('total time elapsed: %s', time.time() - starttime)


def parse_args(version=None):
    """Parse the command line arguments

    Returns the args as an argparse.Namespace object

    """
    import sys
    import argparse


    description = "Parse command line arguements"

    # use formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # so that --help lists the defaults
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument(
    #    'config_file', type=str, nargs='?',
    #    help='the configuration file'
    # )


    parser.add_argument(
        '--nworkers', type=int, action='store', default=0,
        help='integer number of worker processes')

    parser.add_argument(
        '--skip', type=int, action='store', default=0,
        help='skip n tiles')

    parser.add_argument(
        '--count', type=int, action='store', default=0,
        help='count tiles and exit')

    parser.add_argument(
        '--modulo', type=int, action='store', default=0,
        help='process modulo n tiles')

    parser.add_argument(
        '--remainder', type=int, action='store', default=0,
        help='modulus remainder to us: [maximum = MODULO -1 ]')

    parser.add_argument(
        '-l', '--log_file', type=str, action='store', default=None,
        help='filename for storing logging output ' +
        '[default is to stream to stdout]'
    )

    parser.add_argument(
        '--overwrite', action='store_true', default=False,
        help='overwrite existing output data files')

    parser.add_argument(
        '--checkplots', action='store_true', default=False,
        help='create check plots')

    parser.add_argument(
        '--version', action='store_const', default=False, const=True,
        help='show the version')

    parser.add_argument(
        '--tilename', dest='tilename', action='store',
        help='tilename overrides tilename in config file')

    parser.add_argument(
        '--tilelistfile', action='store', default=None,
        help='file contain list of tilenames; astropy table formats')

    parser.add_argument("--debug", action='store_true', default=False,
                        dest='debug', help="debug option; logging level DEBUG")

    parser.add_argument(
        '--decmin', action='store', default=-90.0,
        help='Tile centre minimum  Declination')

    parser.add_argument(
        '--decmax', action='store', default=90.0,
        help='Tile centre maximum  Declination')

    parser.add_argument(
        '--ramin', action='store', default=0.0,
        help='Tile centre minimum  RA')

    parser.add_argument(
        '--ramax', action='store', default=24.0,
        help='Tile centre maximum  RA')

    parser.add_argument(
        '--start_point', type=str, action='store', default='mk_cats',
        choices = ('mk_cats', 'add_DEScat'),
        help='Starting point for processing'
    )


    # WARNING could be a proxy for level 1 verbosity
    parser.add_argument(
        "--quiet", action='store_true', default=False,
        dest='quiet',
        help="NOT IMPLEMENTED YET: quiet option; logging level WARNING")

    parser.add_argument(
        '-v', '--verbosity', type=int, action='store', default=0,
        choices=(0, 1, 2, 3),
        help='NOT IMPLEMENTED YET: integer verbosity level')



    parser.add_argument(
        "--force", action='store_true', default=False,
        help="force processing ignoring lockfile or existing data")

    args = parser.parse_args()

    # if args.config_file is None:
    #    if args.version:
    #        print(version)
    #    else:
    #        parser.print_help()
    #    sys.exit()
    # elif args.version:
    #    print(version)

    return args

if __name__ == '__main__':
    """

    """

    import os
    import time
    import glob

    import multiprocessing
    from multiprocessing import Lock, Process, Queue
    from multiprocessing import current_process, cpu_count

    t0 = time.time()

    nice_level = os.nice(0)
    print('Current level of nice:', nice_level)
    nice_increment = 19 - nice_level
    nice_level = os.nice(nice_increment)
    print('Current level of nice:', nice_level)

    tile = "DES0453-4457"
    tile = "DES2359+0043"
    # File where you want the input coords saved
    coord_file = tile + '_' + "test_input.txt"
    rcore = "5"

    # Folder where you want the output cats saved
    outdir = ""
    config_file = 'wise2des.cfg'
    cfg = parse_config(config_file=config_file, debug=True)
    release = cfg.get("des", "release")
    datapath = cfg.get("des", "datapath")


    tile = None
    single = False
    args = parse_args(version=None)

    checkplots = args.checkplots
    nworkers = args.nworkers
    nskip = args.skip
    modulo = args.modulo
    overwrite = args.overwrite
    remainder = args.remainder
    remainder = min(remainder, modulo - 1)

    ramin = float(args.ramin)
    ramax = float(args.ramax)

    decmin = float(args.decmin)
    decmax = float(args.decmax)

    loglevel = logging.INFO
    if args.debug: loglevel=logging.DEBUG

    create_logger(loglevel=loglevel)
    logger = logging.getLogger()

    if args.tilename is not None:
        single = True
        tilename = args.tilename

    if args.tilelistfile is not None:
        tilelistfile = args.tilelistfile
        print('Reading:', tilelistfile)
        tilelist = Table.read(tilelistfile)
        tilelist.info()

        raw_input("Enter any key to continue: ")



    if single:
        # 24hr edge test
        # if tile is not None:
            # tile = 'DES2327-5248'
            # tile = 'DES2359+0043'
        wise_forced_phot(tilename=tilename, overwrite=overwrite,
                         config_file=config_file)
        raw_input("Enter any key to continue: ")


    # multisite()

    Y1A1_GravLense_Tiles = ['DES2327-5248', 'DES0406-5414']

    SVA1_WL_TestBed_Tiles = ['DES0449-4706', 'DES0453-4706', 'DES0457-4706',
                             'DES0449-4748', 'DES0453-4748', 'DES0457-4748',
                             'DES0449-4831', 'DES0453-4831', 'DES0458-4831',
                             'DES0449-4914', 'DES0454-4914', 'DES0458-4914']

    tilelist = SVA1_WL_TestBed_Tiles
    # tilelist = Y1A1_GravLense_Tiles

    logger.info('Skip n tiles: %s', str(nskip))


    # glob on DESHHMMsDDMM to get a list of tilenames
    tilepath = datapath
    print('glob tilepath: ', tilepath)
    tilepaths = glob.glob(tilepath + '/DES?????????')
    i = -1
    itile = -1
    ntiles = len(tilepaths)
    logger.info('Number of input tiles: %s', str(ntiles))
    tilelist = []
    # apply modulo remainder sceduling for separate processing instances
    for tile in tilepaths[nskip:]:
        itile = itile + 1
        if (modulo == 0) or (modulo > 0 and itile % modulo == remainder):
            i = i + 1
            tilename = tile.split('/')[-1]
            tilelist.append(tilename)
            logger.info('Append to tilelist: %d %s', itile, tilename)

    tilelist = np.sort(tilelist)

    # limit tiles by ramin and ramax from the tilename based directory
    # e.g. DESHHMM+DDSS, DES2348-3056
    #      012345678901
    itile = -1
    ikeep = -1
    tilelist_keep = []
    print('ramin, ramax: ', ramin, ramax)
    print('decmin, decmax: ', decmin, decmax)
    for tile in tilelist:
      itile = itile + 1
      print('ramin, ramax: ', ramin, ramax)
      print('decmin, decmax: ', decmin, decmax)
      print('tilename:', tile)
      ratile = float(tile[3:3+2]) + float(tile[5:5+2])/60.0
      dectile_sign = tile[7:7+1]
      dectile = float(tile[8:8+2]) + float(tile[10:10+2])/60.0
      print('dectile_sign:', dectile_sign)
      if dectile_sign == '-':
          dectile = -1.0 * dectile

      print(itile, tile, ratile, dectile)
      if (ratile <= ramax and ratile >= ramin and
          dectile <= decmax and dectile >= decmin):
          ikeep = ikeep + 1
          tilelist_keep.append(tile)

    nkeep = ikeep + 1
    tilelist = tilelist_keep
    print('Number of tiles in range: ', nkeep)
    print(tilelist)
    print('ramin: ', ramin)
    print('ramax: ', ramax)
    print('Number of tiles in range: ', nkeep)

    raw_input("Enter any key to continue: ")

    itile = 0
    for tile in tilelist:
        itile = itile + 1
        logger.info('tilename: %s, %s', str(itile), tile)

    ntiles = len(tilelist)
    logger.info('Number of tiles to process: %s', str(ntiles))
    raw_input("Enter any key to continue: ")
    # process the data
    if nworkers == 0:
        nworkers = None

    mp_forcephot(tilelist=tilelist, overwrite=overwrite, nworkers=nworkers)

    print('Total time elapsed:', time.time() - t0)
