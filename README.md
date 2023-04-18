Python module that queries Vizier to download photometry and to perform a multi-component SED fitting. All the templates are already packaged, no additional download is necessary.

To install:
```
pip install git+https://github.com/mkounkel/SEDFit.git
```

Example usage:

```
from SEDFit.sed import SEDFit
x=SEDFit('02:03:47.1141597864','+35:35:28.665702672',1)
x.addguesses(r=[2],teff=[10000],logg=3)
x.addrange(logg=[1,4])
x.fullfit()
x.makeplot()
print("Distance: {} pc".format(x.getdist()))
print("AV: {} mag".format(x.getav()))
print("Radius: {} Rsun".format(x.getr()))
print("Teff: {} K".format(x.getteff()))
print("Log g: {} ".format(x.getlogg()))
print("Fe/H: {}".format(x.getfeh()))
print("Chi squared: {}".format(x.getchisq()))
x.sed['model']=x.mags
x.sed
```

<img width="810" alt="image" src="https://user-images.githubusercontent.com/6738789/232641720-92d7f35f-3295-4c7e-b68b-5fbf72e16f77.png">
<img width="685" alt="image" src="https://user-images.githubusercontent.com/6738789/232651588-d3b2bad2-c6b9-499c-9f48-3d9e1f71c670.png">
Distance: 195.97488165738196 pc <br>
AV: 0.26569459140300744 mag  <br>
Radius: [10.298017329573494] Rsun <br>  
Teff: [4523.642995967771] K  <br>
Log g: [3.510695283221381]  <br>
Fe/H: -0.4072943677937134  <br>
Chi squared: 31.84604213491017

----------------------------------

Detailed description


```
x=SEDFit(ra,dec,radius,**kwargs)
```
Initializes SED fitting

Required keywords:<br>
- ra: Right Ascention of the target, can be in decimal or hexadecimal format, float or a string.<br>
- dec: Declination of the target, can be in decimal or hexadecimal format, float or a string.<br>
- radius: Radius in arcseconds around which to search Vizier for photometry.<br>

Optional keywords:<br>
- grid_type: Which theoretical SEDs should be used for fitting. Options: 'kurucz' (default, Kurucz 1992), 'coelho' (Coelho 2014)<br>
- use_gaia_params: Whether Gaia astrometry should be downloaded to estimate distances, true by default, boolean<br>
- use_gaia_xp: Whether Gaia XP spectrum should be used or not, true by default, boolean<br>
- gaia_params: FITS filename where Gaia astrometry is/should be saved, uses coordinates as default, string<br>
- flux_filename: FITS filename where photometry is/should be saved, uses coordinates as default, string<br>
- gaia_filename: FITS filename where Gaia spectrum is/should be saved, uses coordinates as default, string<br>
- download_flux: Whether flux should be (re-)downloaded. False by default if flux_filename file already exists, boolean<br>
- download_gaia: Whether Gaia XP spectrum should be (re-)downloaded. False by default if gaia_filename file already exists, boolean<br>
- gaia: Whether Gaia fluxes (G, BP, RP) should be downloaded, True by default, boolean<br>
- sdss: Whether SDSS fluxes (u,g,r,i,z) should be downloaded, True by default, boolean<br>
- panstarrs: Whether PanSTARRS fluxes (g,r,i,z,y) should be downloaded, True by default, boolean<br>
- johnson: Whether Johnson fluxes (U,B,V,R,I) should be downloaded, True by default, boolean<br>
- cousins: Whether Cousins fluxes (U,B,V,R,I) should be downloaded, True by default, boolean<br>
- tmass: Whether WISE fluxes (W1, W2, W3, W4) should be downloaded, True by default, boolean<br>
- galex: Whether GALEX fluxes (FUV, NUV) should be downloaded, True by default, boolean<br>
- parallax_sigma: Range of distances that should be used in fitting based on Gaia parallax and corresponding uncertainty, default is 3 sigma.

```
x.addguesses()
```
Provides initial guesses to the fitting process

Optional keywords:<br>
- dist: Estimated distance to the system in parsec, 1000 pc by default <br>
- av: Estimated extinction to the system in mag, 0 mag by default <br>
- r: Estimated radii of all of the stars in Rsun, stored as list for all of the stars in the system, [1] Rsun by default<br>
- teff: Estimated Teff of all of the stars in K, stored as list for all of the stars in the system, [5000] K by default<br>
- logg: Estimated log g of all of the stars, stored as list for all of the stars in the system, [4] dex by default<br>
- feh: Estimated Fe/H of the system, 0 dex by default<br>
- alpha: Estimated alpha/H of the system, 0 dex by default<br>
- area: Whether surface area is given instead of radius in r, False by default, boolean<br>

Stipulations:<br>
- Lists for all of the variables should match each other in length.<br>
- If there are multiple stars, but only one supplied value, for, e.g., teff, it is assumed that both stars will be initialized with the same teff<br>
- At least one of the lists needs to explicitly have the length of the desired number of elements<br>
- Nested lists are allowed (e.g., if the arrays of areas, teff, logg are produced by the PHOEBE mesh, they can go inside the list as one of the elements, e.g, r=[[meshes for star 1], [meshes for star2], radius3]). In this case, r is treated as the surface area for the stars where mesh is provided.
- Sources where meshes are provided should be passed as is, without additional fitting.
- Initial guesses for teff, logg, feh, and alpha need to be within the bounds of the grid of synthetic SEDs, otherwise it would not be able to function. Refer to the papers of the included models for description.

```
x.addrange()
```
Provides the allowed range of values that can be returned by the fitting, specifying lower and upper bound

Optional parameters<br>
- dist: Range of distances in pc, default [0, 1e5]<br>
- av: Range of AVs in mag, default [0., los_max.], where los_max is the maximum extinction along the line of sight from Fitzpatrick (1999) map <br>
- r: Range of radii in Rsun, default [0,2000.]<br>
- teff: Range of Teff in K, default [0, 1e6]<br>
- logg: Range of log g, default [-2,10]<br>
- feh: Range of Fe/H, default [-5.,2.]<br>

Stipulations:<br>
- Alpha is included as is, cannot be changed in the model in the fitting process, thus no range for it is provided<br>
- If fitting multiple stars, lists should be nested, e.g, r=[[1,5],[3,4]], this will set the allowed range of radii to be 1-5 Rsun for the first star, and 3-4 Rsun for the second star<br>
- The number of elements for the nested lists should be consistent with the number of stars provided in the addguess() function<br>
- If fitting multiple stars, but providing only two numbers, e.g., r=[1,5], then this range is assumed as valid for all stars<br>
- The ranges can be arbitrary, but keep in mind the bounds of the grid that is being used<br>

```
x.fullfit()
```
Performs the full SED fitting, including r, teff, logg, dist, av, feh

Optional parameters<br>
- fitstar: list of stars that should be used in the fitted process, vs those the parameters of which should be passed directly from the initial guesses unaltered. 1=include, 0=exclude, e.g., [0,0,1] would preserve the parameters of the first two stars and would try to fit the residual fluxes using only the third one. System-wide parameters such as feh, av, dist would be fit in all cases. By default, fits all stars. If meshes are used, fitstar for the stars where this is the case must be set to 0.
- use_gaia: whether to use Gaia XP spectrum in the fitting process, True by default, boolean
- use_mag: list of indices of x.sed table that are used for fitting, as some fluxes could be of poor quality, such as due to saturation, source mismatch, or IR/UV excess. Uses all fluxes by default. E.g, in the example above, Fluxes from SDSS and Pan-STARRS i band are highly discrepant from the model. To improve the quality of the fit and to exclude them,<br> use_mag=np.where((x.sed['sed_filter']!='PAN-STARRS.PS1.i') & (x.sed['sed_filter']!='SDSS.i'))[0]

```
x.fit()
```
An alternative to fullfit, fitting only dist, av, and radii of the stars, the remainig parameters are preserved. 

Optional parameters<br>
- fitstar: same as above
- use_gaia: same as above
- use_mag: same as above
- ratio: Whether the flux ratio between the stars in a binary at one of the passbands is known, to provide a constraint on the radii. False by default, boolean. If True, requires both ratiolambda and fluxratio.
- ratiolambda: if ratio=True, reference wavelength at which the flux ratios are known
- fluxratio: if ratio=True, known flux ratio.

```
x.makeplot()
```
Produces a plot of the SED fit.

Optional parameters<br>
- file: string containing a filename (including the extension) to a file where the figure should be saved.
- getplot: whether to return the editable matplotlib ax variable instead of displaying it. False by default.If True,
ax=x.makeplot(getplot=False)<br>
....<br>
plt.show()

```
x.getchisq()
```
Calculates Chi squared of the fit.

Optional parameters<br>
- idx: list of indices of fluxes that should be used in the chisq calculation (similar format to use_mag in fullfit)
- gaia: whether to use Gaia XP spectrum in chisq calculation, boolean.

```
x.getdist(),x.getav(),x.getr(),x.getteff(),x.getlogg(),x.getfeh()
```
Returns fitted variables. Same as<br>
x.dist, x.av, x.r, x.teff, x.logg, x.feh<br>

Other usefulf variables:<br>
- x.sed: table consisting of the downloaded fluxes from Vizier
- x.gaia: table consisting of the Gaia XP spectrum
- x.mags: array of the model-predicted fluxes for all of the elements in the x.sed table
- x.spec: array of the model-predicted fluxes for all of the elements in the x.gaia table
- x.la: array of wavelength for the model
- x.flux: nested array of the raw fluxes for all of the stars provided by the model (i.e., not scaled by distance, radius, av, etc)
- x.fx: nested array of the fitted model fluxes for all of the stars
- x.f: array of the co-added fluxes for all of the model fitted fluxes
