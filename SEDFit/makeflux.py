from astropy.table import Table
import numpy as np
import astropy.units as u
import warnings
warnings.simplefilter("ignore")

empty=Table.read('emptyflux.fits')
def makeflux(tmass_j=np.nan,tmass_h=np.nan,tmass_k=np.nan,cousins_u=np.nan,cousins_b=np.nan,cousins_v=np.nan,
             cousins_r=np.nan,cousins_i=np.nan,gaia_g=np.nan,gaia_bp=np.nan,gaia_rp=np.nan,galex_fuv=np.nan,galex_nuv=np.nan,
             johnson_u=np.nan,johnson_b=np.nan,johnson_v=np.nan,johnson_r=np.nan,johnson_i=np.nan,panstarrs_g=np.nan,
             panstarrs_r=np.nan,panstarrs_i=np.nan,panstarrs_z=np.nan,panstarrs_y=np.nan,sdss_u=np.nan,sdss_g=np.nan,sdss_r=np.nan,
             sdss_i=np.nan,sdss_z=np.nan,wise_w1=np.nan,wise_w2=np.nan,wise_w3=np.nan,xmm_v=np.nan,xmm_b=np.nan,xmm_u=np.nan,
             xmm_uvw1=np.nan,xmm_uvm2=np.nan,xmm_uvw2=np.nan,spitzer_36=np.nan,spitzer_45=np.nan,spitzer_58=np.nan,spitzer_8=np.nan,
             spitzer_24=np.nan,xpflux=None,xpflux_error=None):
    t=empty.copy()

    def assign(mag,filter,zeropoint=1):
        if np.isfinite(mag).any():
            a=np.where(t['sed_filter']==filter)[0][0]
            if type(mag)==list:
                flux=(10 ** (-mag[0] / 2.5)) * zeropoint
                if not np.isfinite(mag[1]): mag[1]=0.
                eflux=(10 ** (-mag[0] / 2.5)) * zeropoint*mag[1]*np.log(10)/2.5
            else:
                flux=(10 ** (-mag / 2.5)) * zeropoint
                eflux=0.
            if eflux/flux<0.02: eflux=flux*0.02 
            t['eflux'][a]=np.round(eflux/flux/np.log(10),4)
            t['flux'][a]=np.round(np.log10(flux*t['la'][a]),4)
        return

    #vega
    assign(tmass_j,'2MASS.J',3.13e-10)
    assign(tmass_h,'2MASS.H',1.13e-10)
    assign(tmass_k,'2MASS.Ks',4.28e-11)

    #vega
    assign(cousins_u,'Cousins.U',3.49719e-9)
    assign(cousins_b,'Cousins.B',6.72553e-9)
    assign(cousins_v,'Cousins.V',3.5833e-9)
    assign(cousins_r,'Cousins.R',2.23895e-9)
    assign(cousins_i,'Cousins.I',1.19038e-9)

    #vega
    assign(gaia_g,'GAIA.GAIA3.G',2.50386e-9)
    assign(gaia_bp,'GAIA.GAIA3.Gbp',4.07852e-9)
    assign(gaia_rp,'GAIA.GAIA3.Grp',1.26902e-9)

    #ab
    assign(galex_fuv,'GALEX.FUV',4.6194e-8)
    assign(galex_nuv,'GALEX.NUV',2.05634e-8)

    #vega
    assign(johnson_u,'Johnson.U',3.49719e-9)
    assign(johnson_b,'Johnson.B',6.72553e-9)
    assign(johnson_v,'Johnson.V',3.5833e-9)
    assign(johnson_r,'Johnson.R',1.87529e-9)
    assign(johnson_i,'Johnson.I',9.23651e-10)

    #ab
    assign(panstarrs_g,'PAN-STARRS.PS1.g',4.62937e-9)
    assign(panstarrs_r,'PAN-STARRS.PS1.r',2.83071e-9)
    assign(panstarrs_i,'PAN-STARRS.PS1.i',1.91728e-9)
    assign(panstarrs_z,'PAN-STARRS.PS1.z',1.44673e-9)
    assign(panstarrs_y,'PAN-STARRS.PS1.y',1.17434e-9)

    #ab
    assign(sdss_u,'SDSS.u',8.60588e-9)
    assign(sdss_g,'SDSS.g',4.92255e-9)
    assign(sdss_r,'SDSS.r',2.85425e-9)
    assign(sdss_i,'SDSS.i',1.94038e-9)
    assign(sdss_z,'SDSS.z',1.35994e-9)

    #vega
    assign(wise_w1,'WISE.W1',8.02178e-12)
    assign(wise_w2,'WISE.W2',2.36824e-12)
    assign(wise_w3,'WISE.W3',6.39263e-14)

    #ab
    assign(xmm_v,'XMM-OT.V',3.69915e-9)
    assign(xmm_b,'XMM-OT.B',5.78686e-9)
    assign(xmm_u,'XMM-OT.U',9.004e-9)
    assign(xmm_uvw1,'XMM-OT.UVW1',1.26684e-8)
    assign(xmm_uvm2,'XMM-OT.UVM2',2.0104e-8)
    assign(xmm_uvw2,'XMM-OT.UVW2',2.38676e-8)

    #vega
    assign(spitzer_36,'Spitzer.IRAC.3.6',6.57558e-12)
    assign(spitzer_45,'Spitzer.IRAC.4.5',2.65608e-12)
    assign(spitzer_58,'Spitzer.IRAC.5.8',1.04921e-12)
    assign(spitzer_8,'Spitzer.IRAC.8.0',3.14081e-13)
    assign(spitzer_24,'Spitzer.MIPS.24',3.82549e-15)

    a=np.where(np.isfinite(t['flux']))
    t=t[a]

    d=Table()

    if xpflux is not None:
        d['la']=np.linspace(4000,8000,41)*u.AA
        d['flux']=(xpflux*u.W/u.m**2/u.nm).to((u.erg/u.s/(u.cm**2)/u.AA))
        try:
            d['eflux']=(xpflux_error*u.W/u.m**2/u.nm).to((u.erg/u.s/(u.cm**2)/u.AA))
        except:
            d['eflux']=0*(u.erg/u.s/(u.cm**2)/u.AA)
        d['width']=20*u.AA
        d['sed_filter']='XP'
        a=np.where(d['flux']>0)[0]
        d=d[a]
        d['eflux']=d["eflux"]/d["flux"]/np.log(10)
        d['flux'] =np.log10(d["flux"].value*d['la'])
    
    return t,d
