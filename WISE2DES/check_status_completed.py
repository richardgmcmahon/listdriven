from __future__ import print_function, division

from astropy.io import fits


def check_status_completed(tilename=None, datapath=None):
    """
    /data/desardata2/Y1A1/WISE2DES/R4/*/DES*.completed

    returns a list of tiles that are not completed

    could also verify that DES2327-5248_WISEfp_DEScat_WISE_match.fits
    is valid FITS file


    http://docs.astropy.org/en/stable/io/fits/usage/verification.html#verification-using-the-fits-checksum-keyword-convention


    see also: https://github.com/astropy/astropy/tree/master/astropy/io/fits/scripts

    fitscheck.py




    """

    return status


def verify_compliance(filename):
    """Check for FITS standard compliance."""

    hdulist = fits.open(filename)
    try:
        hdulist.verify('exception')
    except fits.VerifyError as exc:
        log.warning('NONCOMPLIANT %r .. %s' %
                 (filename), str(exc).replace('\n', ' '))
        return 1
    return 0

def setup_logging():

    if OPTIONS.verbose:
        log.setLevel(logging.INFO)
    else:
        log.setLevel(logging.WARNING)

    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter('%(message)s'))
    log.addHandler(handler)


def process_file(filename):
    """
    Handle a single .fits file,  returning the count of checksum and compliance
    errors.
    """

    try:
        checksum_errors = verify_checksums(filename)
        if OPTIONS.compliance:
            compliance_errors = verify_compliance(filename)
        else:
            compliance_errors = 0
        if OPTIONS.write_file and checksum_errors == 0 or OPTIONS.force:
            update(filename)
        return checksum_errors + compliance_errors
    except Exception as e:
        log.error('EXCEPTION %r .. %s' % (filename, e))
        return 1


def handle_options(args):
    if not len(args):
        args = ['-h']

    parser = optparse.OptionParser(usage=textwrap.dedent("""
        fitscheck [options] <.fits files...>
        .e.g. fitscheck example.fits
        Verifies and optionally re-writes the CHECKSUM and DATASUM keywords
        for a .fits file.
        Optionally detects and fixes FITS standard compliance problems.
        """.strip()))

    parser.add_option(
        '-k', '--checksum', dest='checksum_kind',
        type='choice', choices=['standard', 'nonstandard', 'either', 'none'],
        help='Choose FITS checksum mode or none.  Defaults standard.',
        default='standard', metavar='[standard | nonstandard | either | none]')

    parser.add_option(
        '-w', '--write', dest='write_file',
        help='Write out file checksums and/or FITS compliance fixes.',
        default=False, action='store_true')

    parser.add_option(
        '-f', '--force', dest='force',
        help='Do file update even if original checksum was bad.',
        default=False, action='store_true')

    parser.add_option(
        '-c', '--compliance', dest='compliance',
        help='Do FITS compliance checking; fix if possible.',
        default=False, action='store_true')

    parser.add_option(
        '-i', '--ignore-missing', dest='ignore_missing',
        help='Ignore missing checksums.',
        default=False, action='store_true')

    parser.add_option(
        '-v', '--verbose', dest='verbose', help='Generate extra output.',
        default=False, action='store_true')

    global OPTIONS
    OPTIONS, fits_files = parser.parse_args(args)

    if OPTIONS.checksum_kind == 'none':
        OPTIONS.checksum_kind = False

    return fits_files


def check_tile(inpath=None, debug=False, fast=False,
    ncheckfiles=None, clean=False):
    """

    """

    import shutil

    if debug: print('Processing:', inpath)
    debug = False
    filelist = glob.glob(inpath + '/*fit*')
    nfiles = len(filelist)
    print('Number of fits files:', nfiles)

    if fast and clean and nfiles != ncheckfiles and nfiles != 0:
        reply = raw_input("Delete directory "+ inpath +"? [y/[n]] ")
        if reply == 'y' or reply == 'Y':
            shutil.rmtree(inpath)

        filelist = glob.glob(inpath + '/*fit*')
        nfiles = len(filelist)
        print('Number of fits files:', nfiles)

    if not fast:
        for infile in filelist:
            if debug: print(infile)

            hdulist = fits.open(infile)
            nhdu = len(hdulist)
            if debug: print('Number of extension:', nhdu)

            for ihdu in xrange(0, nhdu):
                if hdulist[ihdu].verify() is not None:
                   print(ihdu, hdulist[ihdu].verify())

    if fast and nfiles != ncheckfiles:
        for infile in filelist:
            if debug: print(infile)

            hdulist = fits.open(infile)
            nhdu = len(hdulist)
            if debug: print('Number of extension:', nhdu)

            for ihdu in xrange(0, nhdu):
                if hdulist[ihdu].verify() is not None:
                   print(ihdu, hdulist[ihdu].verify())

    return nfiles


if __name__ == '__main__':
    """

    """
    import os
    import time
    import glob

    import argparse

    description = "Parse command line arguements"
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '--skip', action='store', default=None,
        help='Starting point for processing')

    parser.add_argument(
        '--fast', action='store_true', default=False,
        help='Fast option than checks number of files and only verifies the final output file')

    parser.add_argument(
        '--nfiles', action='store', default=13,
        help='Expected number of files')


    parser.add_argument("--debug", action='store_true', default=False,
        dest='debug', help="debug option; logging level DEBUG")

    parser.add_argument(
        '--clean', action='store_true', default=False,
        help='Clean/Delete incomplete directories')

    args = parser.parse_args()

    clean = args.clean
    debug = args.debug
    fast = args.fast
    ncheckfiles = args.nfiles
    skip = args.skip

    inpath = "/data/desardata2/Y1A1/WISE2DES/R4/"

    pathlist = sorted(glob.glob(inpath + 'DES?????????'))
    nlist = len(pathlist)
    print('Number of tiles:', nlist)
    itile = 0
    if skip is not None:
        iskip = int(skip)
        pathlist = pathlist[iskip:]
        itile = iskip

    for path in pathlist:
        itile = itile + 1
        print('Processing:', itile, path)
        nfiles = check_tile(inpath=path, fast=fast,
            ncheckfiles=ncheckfiles, clean=clean)

        if debug:
            print(path.split('/'))
        tile = path.split('/')[-1]
        if debug:
            print('tile:', tile)

        if nfiles > 0:
            filename_suffix = '_WISEfp_DEScat_WISE_match'

            infile = inpath + tile + '/' + tile + filename_suffix + '.fits'
            print('Read:', infile)

            hdulist = fits.open(infile)
            nhdu = len(hdulist)
            if debug: print('Number of extension:', nhdu)

            for ihdu in xrange(0, nhdu):
                if hdulist[ihdu].verify() is not None:
                    print(ihdu, hdulist[ihdu].verify())

    filename = 'DES2327-5248.completed'
