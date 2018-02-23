
First Header | Second Header
------------ | -------------
Content from cell 1 | Content from cell 2
Content in the first column | Content in the second column
 - | test 


WISE FWHM:


http://wise2.ipac.caltech.edu/docs/release/prelim/expsup/sec4_5c.html#psf


     (1) | (2) | (3) | (4) | (5) 
   ----- | --- | --- | --- | ----
-  |  Airy  | Airy     |        FWHM | 
-  | radius | diameter |       |
   |  (")   |  |   |
W1 |  2.11  |  4.22 |  6.1 |  6.1 |  5.8
W2 |  2.89  |  5.78 |  6.4 |  6.8 |  6.4
W3 |  7.27  | 14.54 |  6.5 |  7.4 |  6.6
W4 | 13.90  | 27.80 | 12.0 |   -  | 11.9


Optimal detection that maximise Signal/Noise for point sources
Irwin(1985, CASU) radius = FWHM (diameter = 2 * FWHM) is recommended
Masci(2008)       radius = 0.67 FWHM

(1,2)) Rayleigh diffraction limit (Airy radius) 1.22 lamdba/telescope diameter
(3) Wright+2010 http://adsabs.harvard.edu/abs/2010AJ....140.1868W 
(4) ALLWISE: http://wise2.ipac.caltech.edu/docs/release/allwise/expsup/sec1_2.html
(5) Table 1 in Aniano+20


CASU imcore aperture definitions:
http://casu.ast.cam.ac.uk/surveys-projects/vista/technical/catalogue-generation

The radii are:
1:  1/2 rcore
2:  1/sqrt(2) rcore
3:  rcore
4:  sqrt(2)  rcore
5:  2 rcore
6:  2 sqrt(2) rcore
7:  4 rcore
8:  5  rcore
9:  6  rcore
10: 7  rcore
11: 8  rcore
12: 10 rcore
13: 12 rcore




DES Y1A1+ pixel size is 0.263 arcsec per pixel 

     (1) | (2) | (3) | (4) | (5) |
	 ---
DES FWHM(1)     rcore        rcore 
     "       pixels arcsec
g   0.9       2.0   0.526    3.0  0.789
r   0.9       2.0   
i   0.9       2.0
z   0.9       2.0
Y   0.9       2.0


Note CASU Classify only computes APCOR upto APER 7 ie. APCOR7; it assumes
that this recovers 100% of flux and APER8+ are for large extended objects like
galaxies.

rcore = 2
Aper2(radius) = 
Aper3(radius) = 0.53"
Aper4(radius) = 
Aper5(radius) = 1.05    (2 * rcore) [should be best for DES point sources]
Aper6(radius) = 1.49    (2 * 1.4 rcore)
Aper7(radius) = 2.14"   (4 * rcore)
Aper8(radius) = 2.63"   (5 * rcore)
Aper9(radius) = 3.16"   (6 * rcore)



rcore = 3
Aper2(radius) = 
Aper3(radius) = 0.79"
Aper4(radius) = 
Aper5(radius) = 1.59    (2 * rcore) 
Aper6(radius) =         (2 * 1.4 rcore)
Aper7(radius) = 3.16"   (4 * rcore)
Aper8(radius) = 3.95"   (5 * rcore)

