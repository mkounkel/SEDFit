Python module to downloads photometry and to performs a multi-component SED fitting. All the templates are already packaged, no additional download is necessary.

Example usage:

```
from SEDFit.sed import SEDFit
x=SEDFit(86.4236567891,-0.06733864568,1)
x.addguesses(dist=166,av=0.,r=[2],teff=[10000],logg=4)
x.addrange(dist=[100,200],r=[0.5,5],logg=[3,5])
x.fullfit()
x.makeplot()
print("Distance: {} pc".format(x.getdist()))
print("AV: {} mag".format(x.getav()))
print("Radius: {} Rsun".format(x.getr()))
print("Teff: {} K".format(x.getteff()))
print("Log g: {} ".format(x.getlogg()))
print("Fe/H: {}".format(x.getfeh()))
```

<img width="810" alt="image" src="https://user-images.githubusercontent.com/6738789/232641720-92d7f35f-3295-4c7e-b68b-5fbf72e16f77.png">
Distance: 195.97488165738196 pc <br>
AV: 0.26569459140300744 mag  <br>
Radius: [10.298017329573494] Rsun <br>  
Teff: [4523.642995967771] K  <br>
Log g: [3.510695283221381]  <br>
Fe/H: -0.4072943677937134  


----------------------------------

Detailed description

Initialization:
```
x=SEDFit(ra,dec,radius,**kwargs)
```
Required keywords:<br>
ra: Right Ascention of the target, can be in decimal or hexadecimal format, float or a string.<br>
dec: Declination of the target, can be in decimal or hexadecimal format, float or a string.<br>
radius: Radius in arcseconds around which to search Vizier for photometry.<br>

Optional keywords:<br>
grid_type: Which theoretical SEDs should be used for fitting. Options: 'kurucz' (default, Kurucz 1992), 'coelho' (Coelho 2014)<br>
use_gaia: Whether Gaia XP spectrum should be used or not, true by default, boolean<br>
flux_filename: FITS filename where photometry is/should be saved, string<br>
gaia_filename: FITS filename where Gaia spectrum is/should be saved, string<br>
download_flux: Whether flux should be (re-)downloaded. False by default if flux_filename file already exists, boolean<br>
download_gaia: Whether Gaia spectrum should be (re-)downloaded. False by default if gaia_filename file already exists, boolean<br>
gaia: Whether Gaia fluxes (G, BP, RP) should be downloaded, True by default, boolean<br>
sdss: Whether SDSS fluxes (u,g,r,i,z) should be downloaded, True by default, boolean<br>
panstarrs: Whether PanSTARRS fluxes (g,r,i,z,y) should be downloaded, True by default, boolean<br>

johnson: Whether Johnson fluxes (U,B,V,R,I) should be downloaded, True by default, boolean<br>

cousins: Whether Cousins fluxes (U,B,V,R,I) should be downloaded, True by default, boolean<br>

tmass: Whether WISE fluxes (W1, W2, W3, W4) should be downloaded, True by default, boolean<br>

galex: Whether GALEX fluxes (FUV, NUV) should be downloaded, True by default, boolean<br>

```
x.addguesses()
```

