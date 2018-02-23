def add_kernel_DEScat(tile, outdir, kernel_path, suffix_for_kernelcats):

	from astropy.table import Table, hstack
	import numpy as np
	import match_lists
	import astropy.io.fits as fits
	import subprocess
	import os
	config = ConfigParser.RawConfigParser()
	config.read(config_file)

	tile_header = fits.open('/data/desardata/Y1A1/'+str(tile)+'/'+str(tile)+'_g.fits.fz')[1].header
	RA_centre = tile_header['CRVAL1']
	DEC_centre = tile_header['CRVAL2']
	colours = ["g","r","i","z","Y"]
	for colour in colours:
		if not os.path.exists('/data/desardata/Y1A1/'+str(tile)+'/'+str(tile)+'_'+str(colour)+'.fits'):
			subprocess.call(['bash','/data/cl522/WISE2DES/kernelised_extraction/funpack.sh','/data/desardata/Y1A1/'+str(tile)+'/'+str(tile)+'_'+str(colour)+'.fits.fz'])

	subprocess.call(['bash','/data/cl522/WISE2DES/kernelised_extraction/sextractor_script.sh',str(tile),str(RA_centre),str(DEC_centre),str(kernel_path),str(suffix_for_kernelcats),str(outdir)])

	for (n, band) in enumerate(["G", "R", "I", "Z"]): #, "Y"]): removed because final fits file must have < 1000 columns :/
		if n!=4:
			t = Table.read(outdir+"/"+tile+"_"+band.lower()+"_cat"+suffix_for_kernelcats+".fits")
		if n==4:
			t = Table.read(outdir+"/"+tile+"_"+band+"_cat"+suffix_for_kernelcats+".fits")
		for col in t.columns:
			t.rename_column(col, col + "_" + band)
		if n == 0:
			t_cat = t
		else:
			t_cat = hstack([t_cat, t])
	t_cat.write(outdir+"/"+tile+suffix_for_kernelcats+".fits", overwrite = True)

	release = config.get("des", "release")
	t_cat = Table.read(outdir+"/"+tile+suffix_for_kernelcats+".fits")
	t_fp = Table.read(outdir + tile + "_WISEfp.fits")

	dists, inds = match_lists.match_lists(t_cat["ALPHAWIN_J2000_G"].data, t_cat["DELTAWIN_J2000_G"].data, t_fp["RA_CALC_G"].data, t_fp["DEC_CALC_G"].data, 5.0/3600.0, 1)

	ids = np.where( (inds <> len(t_fp)) )[0]
	t_cat = t_cat[ids]
	t_fp = t_fp[inds[ids]]
	t_out = hstack([t_fp, t_cat])
	t_out.write(outdir + "/" + tile + "_WISEfp_DEScat"+suffix_for_kernelcats+".fits", overwrite = True)