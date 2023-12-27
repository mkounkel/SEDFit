import astropy.units as u
import astropy.coordinates as coord
from astropy.io import fits
from astropy.table import Table,vstack
from astropy.modeling import models
from dust_extinction.averages import GCC09_MWAvg,G21_MWAvg
from dust_extinction.parameter_averages import F99
from scipy.optimize import curve_fit
from scipy.interpolate import RegularGridInterpolator
from astroquery.gaia import Gaia

import dustmaps
from dustmaps.sfd import SFDQuery

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import glob
import pickle
from tqdm import tqdm
import pkg_resources
import warnings
import bz2
import tarfile

from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', AstropyWarning)

class SEDFit:
    def __init__(self,ra,dec,radius,frame='icrs',flux_filename='',gaia_filename='',gaia_params='',parallax_sigma=3,
                 download_flux=False,download_gaia=False,use_gaia_params=True,use_gaia_xp=True,nstar=1,**kwargs):
        
        crd=(str(ra)+'_'+str(dec)+'_'+str(radius)).replace(':','_')
        if flux_filename=='':
            flux_filename=crd+'_sed.fits'
        if gaia_filename=='':
            gaia_filename=crd+'_gaiaxp.fits'
        if gaia_params=='':
            gaia_params=crd+'_gaiaparams.fits'
        
        if (' ' in str(ra)) or (':' in str(ra)): raunit=u.hourangle
        else: raunit=u.deg
        c=coord.SkyCoord(ra=ra, dec=dec,unit=(raunit, u.deg),frame='icrs')
        self.c=c
        self.ra=c.ra.value
        self.dec=c.dec.value
        self.radius=radius
        self.area=False
        
        self.getmaxreddening(c)
        
        if (not download_flux) & (len(glob.glob(flux_filename))>0):
            self.sed=Table.read(flux_filename)
        else:
            self.downloadflux(**kwargs)
            if len(self.sed)==0:
                print("Can't download SED for this source. Try again later - if the issue persists no star may be found at this position, or the radius is too large.")
                return
            else:
                self.sed.write(flux_filename,overwrite=True)
                
        
        if use_gaia_params:
            if (len(glob.glob(gaia_params))>0):
                self.gaiaparams=Table.read(gaia_params)
            else:
                try:
                    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
                    j = Gaia.cone_search_async(c, self.radius*u.arcsec,verbose=False)
                    self.gaiaparams=Table(j.get_results())
                    for col in self.gaiaparams.columns:
                        if self.gaiaparams[col].dtype != object: continue
                        self.gaiaparams[col] = [str(n) for n in self.gaiaparams[col]]
                        
                    a=np.argsort(self.gaiaparams['dist'])
                    self.gaiaparams=self.gaiaparams[a][0:1]
                    
                    if len(self.gaiaparams)>0:
                        self.gaiaparams.write(gaia_params)
                    else:
                        self.gaiaparams = Table()
                        self.gaiaparams['parallax']=[np.nan]
                        print('Could not find Gaia astrometry for this source')
                except:
                    self.gaiaparams = Table()
                    self.gaiaparams['parallax']=[np.nan]
                    print('Could not find Gaia astrometry for this source')
        else:
            self.gaiaparams = Table()
            self.gaiaparams['parallax']=[np.nan]
        
        if (use_gaia_xp) & (use_gaia_params):
            if (not download_gaia) & (len(glob.glob(gaia_filename))>0):
                self.gaia=Table.read(gaia_filename)
            else:
                try:
                    self.gaia=self.getgaia(c)
                    self.gaia.write(gaia_filename,overwrite=True)
                    print('Gaia XP spectra available')
                except:
                    print('No Gaia XP spectra found')
                    self.gaia=[]
        else:
            self.gaia=[]
            
        
        self.add_new_grid(**kwargs)
        
        self.extuv=GCC09_MWAvg()
        self.extir=G21_MWAvg()
        
        self.refuv=np.where(self.la<1.1*u.micron)
        self.refir=np.where(self.la>=1.1*u.micron)
        
        self.nstar=int(nstar)
        
        if np.isfinite(self.gaiaparams['parallax'][0]):
            d=1000/self.gaiaparams['parallax'][0]
            ed=[(1000/(self.gaiaparams['parallax']+parallax_sigma*self.gaiaparams['parallax_error']))[0],(1000/(self.gaiaparams['parallax']-parallax_sigma*self.gaiaparams['parallax_error']))[0]]
            ruwe=self.gaiaparams['ruwe'][0]
            print('Gaia DR3 distance towards this source is ' + str(np.round(d,1)) +' pc')
            print(str(parallax_sigma)+' sigma uncertainty in distance is ' + str(np.round(ed[0],1)) +' - '+str(np.round(ed[1],1)) +' pc' )
            if ruwe<1.4:
                print('RUWE is ' + str(np.round(ruwe,3)) )
            else:
                print('RUWE is ' + str(np.round(ruwe,3)) +', distance measurement may be unreliable, proceed with caution')
        else:
            print('Gaia parallax is not available, distance is arbitrarily set by default and should be manually adjusted.')
            d=1000.
            ed=[0,1e5]
        print('Maximum Av along the line of sight is '+str(np.round(self.maxav,3)))
        self.addrange(dist=ed,av=[0.,self.maxav],r=[0.,2000.],teff=[0.,1e6],logg=[-2.,10.],feh=[-5.,2.])
        self.addguesses(dist=d,av=0.,r=[1.]*nstar,teff=[5000.]*nstar,logg=[4.]*nstar,feh=0.,alpha=0.)
        
        
    def getmaxreddening(self,coords):
        a=dustmaps.std_paths.data_dir()
        if len(glob.glob(a+'/sfd'))==0:
            dustmaps.sfd.fetch()
        sfd = SFDQuery()
        self.maxav = sfd(coords)*3.1
        return
        
        
    def downloadflux(self,**kwargs):

        target=str(self.ra)+'%20'+str(self.dec)
        try:
            self.sed=Table.read(f"https://vizier.cds.unistra.fr/viz-bin/sed?-c={target}&-c.rs={self.radius}")
        except:
            self.sed=[]
            return
        self.sed['la']=self.sed['sed_freq'].to(u.AA, equivalencies=u.spectral())
        a=np.where((self.sed['la']<15*u.micron) & (self.sed['la']>1000*u.AA))[0]
        self.sed=self.sed[a]
        
        self.sed["sed_flux"]=self.sed["sed_flux"].to((u.erg/u.s/(u.cm**2)/u.AA),equivalencies=u.spectral_density(self.sed['la'].data*u.AA))
        self.sed["sed_eflux"]=self.sed["sed_eflux"].to((u.erg/u.s/(u.cm**2)/u.AA),equivalencies=u.spectral_density(self.sed['la'].data*u.AA))
        
        a=np.where(self.sed['sed_flux']>0)[0]
        self.sed=self.sed[a]
        
        a=np.where((np.isnan(self.sed["sed_eflux"])==True) | (self.sed["sed_eflux"]/self.sed["sed_flux"]<0.02))[0]
        self.sed["sed_eflux"][a]=self.sed["sed_flux"][a]*0.02
        
        self.sed['eflux']=self.sed["sed_eflux"]/self.sed["sed_flux"]/np.log(10)
        self.sed['flux'] =np.log10(self.sed["sed_flux"]*self.sed['la'])
        
        self.definefilter(**kwargs)
        a=np.argsort(self.sed['la'])
        self.sed=self.sed[a]
        self.sed=self.sed[['sed_filter','la','width','flux','eflux']]
        
        
    def definefilter(self,tmass=True,cousins=True,gaia=True,galex=True,johnson=True,panstarrs=True,sdss=True,wise=True,**kwargs):
        idx=[]
        self.sed['width']=0*u.AA
        if tmass:
            filters=['2MASS:J','2MASS:H','2MASS:Ks']
            width=np.array([0.15,0.24,0.25])/2*u.micron
            idx.extend(self.selectflux(filters,width))
        if cousins:
            filters=['Cousins:U','Cousins:B','Cousins:V','Cousins:R','Cousins:I']
            width=np.array([0.0639,0.0928,0.0843,0.1297,0.095])/2*u.micron
            idx.extend(self.selectflux(filters,width))
        if gaia:
            filters=['GAIA/GAIA3:G','GAIA/GAIA3:Gbp','GAIA/GAIA3:Grp']
            width=np.array([0.4053,0.2158,0.2924])/2*u.micron
            idx.extend(self.selectflux(filters,width))
        if galex:
            filters=['GALEX:FUV','GALEX:NUV']
            width=np.array([0.0269,0.0616])/2*u.micron
            idx.extend(self.selectflux(filters,width,checkpos=False))
        if johnson:
            filters=['Johnson:U','Johnson:B','Johnson:V','Johnson:R','Johnson:I']
            width=np.array([0.0619,0.0891,0.0818,0.1943,0.2176])/2*u.micron
            idx.extend(self.selectflux(filters,width))
        if panstarrs:
            filters=['PAN-STARRS/PS1:g','PAN-STARRS/PS1:r','PAN-STARRS/PS1:i','PAN-STARRS/PS1:z','PAN-STARRS/PS1:y']
            width=np.array([0.1166,0.1318,0.1243,0.09658,0.06149])/2*u.micron
            idx.extend(self.selectflux(filters,width))
        if sdss:
            filters=['SDSS:u','SDSS:g','SDSS:r','SDSS:i','SDSS:z']
            width=np.array([0.0555,0.1245,0.1262,0.1291,0.1326])/2*u.micron
            idx.extend(self.selectflux(filters,width))
        if wise:
            filters=['WISE:W1','WISE:W2','WISE:W3']
            width=np.array([0.66,1.04,5.51])/2*u.micron
            idx.extend(self.selectflux(filters,width,checkpos=False))
            
        self.sed=self.sed[idx]
        for i in range(len(self.sed)):
            self.sed['sed_filter'][i]=self.sed['sed_filter'][i].replace(':','.').replace('/','.')
        return
            
    def selectflux(self,filters,width,checkpos=True):
        idx,dist=[],[]
        for i in range(len(filters)):
            a=np.where(self.sed['sed_filter']==filters[i])[0]
            self.sed['width'][a]=np.round(width[i].to(u.AA))
            if len(a)>0:
                d=np.sqrt((self.ra-self.sed['_RAJ2000'][a])**2+(self.dec-self.sed['_DEJ2000'][a])**2)
                b=np.argmin(d)
                idx.append(a[b])
                dist.append(d[b])
        idx=np.array(idx)
        if (checkpos) & (len(dist)>0):
            b=np.min(dist)
            a=np.where(dist==b)[0]
            idx=idx[a]
        return idx
    def getgaia(self,c):
        sid=self.gaiaparams['source_id']
        
        datalink = Gaia.load_data(ids=sid, data_release = 'Gaia DR3', retrieval_type='XP_SAMPLED',
                                  data_structure = 'INDIVIDUAL', verbose = False, output_file = None)
        dl_keys  = [inp for inp in datalink.keys()]
        dl_keys.sort()
        i=list(datalink.keys())[0]
        d=datalink[i][0].to_table() 
        d['wavelength']=np.round(d['wavelength'].to(u.AA))
        d['flux']=d['flux'].to((u.erg/u.s/(u.cm**2)/u.AA))
        d['flux_error']=d['flux_error'].to((u.erg/u.s/(u.cm**2)/u.AA))
        d['wavelength'].name='la'
        d['flux_error'].name='eflux'
        d['width']=20*u.AA
        d['sed_filter']='XP'
        
        d['eflux']=d["eflux"]/d["flux"]/np.log(10)
        d['flux'] =np.log10(d["flux"]*d['la'])
        return d

    def load_coelho_sed(self,name):
        hdul = fits.open(name)
        header=hdul[0].header
        la=10**(np.arange(header['NAXIS1'])*header['CDELT1']+header['CRVAL1'])*u.AA
        spec = hdul[0].data*np.pi*u.erg/u.s/(u.cm**2)/u.AA
        hdul.close()
        
        bb = models.BlackBody(temperature=header['teff']*u.K)
        l=np.arange(100000,320000,5000)*u.AA
        b=((bb(l)*u.sr).to(u.erg/u.s/(u.cm**2)/u.AA,equivalencies=u.spectral_density(l)))
        b=b/b[0]*spec[-1]
        la=np.append(la,l[1:])
        spec=np.append(spec,b[1:])
        
        a=np.where((la>1000*u.AA))[0]
        return [np.round(np.log10((spec[a]*la[a]).value),3)],la[a],[header['TEFF']],[header['LOG_G']],[header['FEH']],[header['AFE']]

    def load_phoenix_sed(self,name,la,splits):
        hdul = fits.open(name)
        spec = hdul[0].data*u.erg/u.s/(u.cm**3)
        spec=spec.to(u.erg/u.s/(u.cm**2)/u.AA)
        header=hdul[0].header
        hdul.close()
        
        bb = models.BlackBody(temperature=header['PHXTEFF']*u.K)
        ll=np.arange(55000,350000,10)*u.AA
        b=((bb(ll)*u.sr).to(u.erg/u.s/(u.cm**2)/u.AA,equivalencies=u.spectral_density(ll)))
        b=b/b[0]*spec[-1]
        spec=np.append(spec,b[1:])
        
        l=10**np.arange(np.log10(913),np.log10(320000),0.003)*u.AA
            
        s=[]
        for j in splits:
            s.append(np.sum(spec[j].value)/len(j))
        
        return [np.round(np.log10((np.array(s)*l).value),3)],l,[header['PHXTEFF']],[header['PHXLOGG']],[header['PHXM_H']],[header['PHXALPHA']]

    def load_btsettl_sed(self,name):
        n=name.split('/')[-1].replace('m','-').replace('p','+')
        t=Table.read(name)
        la=t['WAVELENGTH'].data*u.AA
        x=t.keys()
        x.pop(0)
        spec,logg,feh,teff,alpha=[],[],[],[],[]
        tt=int(n.split('_')[1].split('.')[0])
        ff=int(n[7:10])/10
        
        l=10**np.arange(np.log10(913),np.log10(320000),0.003)*u.AA
        
        ll=np.diff(l.value)
        d=[ll[0]]
        d.extend(ll)
        d.append(ll[-1])
        d=np.array(d)
        l1=l.value-d[:-1]*0.5
        l2=l.value+d[1:]*0.5

        splits=[]
        for j in range(len(l)):
            splits.append(np.where((la.value>l1[j]) & (la.value<l2[j]))[0])
        
        for i in x:
            if np.sum(t[i].data)>0:
                s=[]
                for j in splits:
                    s.append(np.sum(t[i].data[j])/len(j))
                s=np.array(s)*u.erg/u.s/(u.cm**2)/u.AA*l
                spec.append(np.round(np.log10(s.value),3))
                logg.append(int(i[1:])/10)
                teff.append(tt)
                feh.append(ff)
                alpha.append(0)
                
        return spec,l,teff,logg,feh,alpha

    def load_kurucz_sed(self,name):
        n=name.split('/')[-1].replace('m','-').replace('p','+')
        t=Table.read(name)
        la=t['WAVELENGTH'].data*u.AA
        x=t.keys()
        x.pop(0)
        spec,logg,feh,teff,alpha=[],[],[],[],[]
        tt=int(n.split('_')[1].split('.')[0])
        ff=int(n[1:4])/10
        
        a=np.where((la<10*u.micron) & (la>950*u.AA))[0]
        
        l=np.arange(100000,320000,1000)*u.AA
        bb = models.BlackBody(temperature=tt*u.K)
        b=((bb(l)*u.sr).to(u.erg/u.s/(u.cm**2)/u.AA,equivalencies=u.spectral_density(l)))
        
        
        for i in x:
            if np.sum(t[i].data)>0:
                s=t[i].data[a]*u.erg/u.s/(u.cm**2)/u.AA*la[a]
                s=np.append(s,b[1:]/b[0]*s[-1])
                spec.append(np.round(np.log10(s.value),3))
                logg.append(int(i[1:])/10)
                teff.append(tt)
                feh.append(ff)
                alpha.append(0)
                
        la=np.append(la[a],l[1:])
        return spec,la,teff,logg,feh,alpha

    def add_new_grid(self,grid_path='',grid_pickle='',grid_type='kurucz',loadgrid=True,laname=None,compressed=True,**kwargs):
        
        if grid_pickle=='':
            grid_pickle='sed_'+grid_type+'.p'
            fullpath=pkg_resources.resource_filename('SEDFit', grid_pickle)
        else:
            fullpath=grid_pickle
        if loadgrid & len(glob.glob(fullpath))>0:
            if compressed:
                flux,t,la = pickle.load(bz2.BZ2File(fullpath,'rb'))
            else:
                flux,t,la = pickle.load(open(fullpath,'rb'))
        else:
            path= glob.glob(grid_path)
            
            if len(path)>1:
                a=np.argsort(path)
                path=np.array(path)[a]
            elif len(path)==1:
                path=np.array(path)
            else:
                raise Exception("Can't find SED templates, check grid_path") 
                
            t = Table(names=('teff','logg','feh','alpha'))
            flux=[]
            
            if grid_type=='phoenix':
                    if laname is None:
                        raise Exception("Please provide wavelength array path (laname)")
                    hdul = fits.open(laname)
                    la=hdul[0].data*u.AA
                    hdul.close()
                    l=np.arange(55000,350000,10)*u.AA
                    la=np.append(la,l[1:])
                    
                    l=10**np.arange(np.log10(913),np.log10(320000),0.003)*u.AA
                    
                    ll=np.diff(l.value)
                    d=[ll[0]]
                    d.extend(ll)
                    d.append(ll[-1])
                    d=np.array(d)
                    l1=l.value-d[:-1]*0.5
                    l2=l.value+d[1:]*0.5
                    
                    splits=[]
                    for j in range(len(l)):
                        splits.append(np.where((la.value>l1[j]) & (la.value<l2[j]))[0])
                        
            for i in tqdm(path):
                if grid_type=='phoenix':
                    f,la,teff,logg,feh,alpha=self.load_sed(grid_type,i,la=la,splits=splits)
                else:
                    f,la,teff,logg,feh,alpha=self.load_sed(grid_type,i)
                for j in range(len(f)):
                    flux.append(f[j])
                    t.add_row((teff[j],logg[j],feh[j],alpha[j]))  
                    
            flux=np.array(flux)
            a=np.where(np.isfinite(flux)==False)
            b=np.where(np.isfinite(flux)==True)
            flux[a]=np.nanmin(flux[b])
                
            if compressed:
                pickle.dump( [flux,t,la], bz2.BZ2File(grid_pickle,'wb') )
            else:
                pickle.dump( [flux,t,la], open(grid_pickle,'wb') )
            
        self.interp=self.makeinterp(t,flux)
        self.la=la
        self.grid_table=t
        return
    def load_sed(self,grid_type,i,la=None,splits=None):
        if grid_type=='coelho':
            return self.load_coelho_sed(i)
        elif grid_type=='kurucz':
            return self.load_kurucz_sed(i)
        elif grid_type=='btsettl':
            return self.load_btsettl_sed(i)
        elif grid_type=='phoenix':
            return self.load_phoenix_sed(i,la,splits)
        else:
            raise Exception("Unsupported grid type, available options are 'btsettl', 'coelho', 'kurucz', and 'phoenix'")

    def makeinterp(self,t,flux):
        teff=np.unique(np.array(t['teff']))
        logg=np.unique(np.array(t['logg']))
        feh=np.unique(np.array(t['feh']))
        alpha=np.unique(np.array(t['alpha']))
        if len(alpha)==1:
            d=t.copy()
            t['alpha']=alpha-0.5
            d['alpha']=alpha+0.5
            t=vstack([t,d])
            alpha=np.unique(np.array(t['alpha']))
            flux=np.vstack([flux,flux])
            
        data=np.zeros((len(teff),len(logg),len(feh),len(alpha),len(flux[0])))
        for d in range(len(alpha)):
            q1=np.abs(t['alpha']-alpha[d])
            w1=np.where(q1==np.min(q1))[0]
            for c in range(len(feh)):
                q2=np.abs(t['feh'][w1]-feh[c])
                w2=w1[np.where(q2==np.min(q2))[0]]
                for b in range(len(teff)):
                    q3=np.abs(t['teff'][w2]-teff[b])
                    w3=w2[np.where(q3==np.min(q3))[0]]
                    for a in range(len(logg)):
                        q4=np.abs(t['logg'][w3]-logg[a])
                        w4=w3[np.argmin(q4)]
                        data[b,a,c,d,:]=flux[w4]
        interp = RegularGridInterpolator((teff,logg,feh,alpha), data,bounds_error=False,fill_value=0.)
        return interp
        
    
    def addguesses(self,dist=None,av=None,r=None,teff=None,logg=None,feh=None,alpha=None,area=None,nstar=1,verbose=True):
        if dist is None: dist=self.dist
        if av is None: av=self.av
        if r is None: r=self.r
        if teff is None: teff=self.teff
        if logg is None: logg=self.logg
        if feh is None: feh=self.feh
        if alpha is None: alpha=self.alpha
        if area is None: area=self.area
        
        if not type(r) in [list,tuple,np.ndarray]: r=[r]
        if not type(teff) in [list,tuple,np.ndarray]: teff=[teff]
        if not type(logg) in [list,tuple,np.ndarray]: logg=[logg]
            

        self.nstar=int(np.max([len(r),len(teff),len(logg),nstar]))
        
        if len(r)!=self.nstar: self.r=r*self.nstar
        else: self.r=r
        if len(teff)!=self.nstar: self.teff=teff*self.nstar
        else: self.teff=teff
        if len(logg)!=self.nstar: self.logg=logg*self.nstar
        else: self.logg=logg

        self.dist=dist
        self.av=av
        self.feh=feh
        self.alpha=alpha
        self.area=area
        
        
        self.flux=[]
        for i in range(len(self.r)):
            spec=self.interp((self.teff[i],self.logg[i],self.feh,self.alpha))
            
            if np.sum(spec)==0:
                if verbose:
                    print('Initial guesses for star '+str(i+1)+' are outside of the grid edges; using blackbody instead')
                bb = models.BlackBody(temperature=teff[i]*u.K)
                bx=((bb(self.la)*u.sr).to(u.erg/u.s/(u.cm**2)/u.AA,equivalencies=u.spectral_density(self.la)))
                spec=np.log10(bx.value*self.la.value*np.pi)
            self.flux.append(spec)
        self.f,self.fx,self.mags,self.spec=self.getfluxsystem(self.dist,self.av,self.r)
        self.addrange()
        
        
        return
    

    def addrange(self,dist=None,av=None,r=None,teff=None,logg=None,feh=None):
        if dist is None: dist=self.distrange
        if av is None: av=self.avrange
        if r is None: r=self.rrange
        if teff is None: teff=self.teffrange
        if logg is None: logg=self.loggrange
        if feh is None: feh=self.fehrange
        
        r=np.array(r)
        teff=np.array(teff)
        logg=np.array(logg)
        
        if r.shape[0]==r.size: r=[r]
        if teff.shape[0]==teff.size: teff=[teff]
        if logg.shape[0]==logg.size: logg=[logg]
        
        if len(r)<self.nstar: r=[r[0]]*self.nstar
        if len(teff)<self.nstar: teff=[teff[0]]*self.nstar
        if len(logg)<self.nstar: logg=[logg[0]]*self.nstar
        
        if len(r)>self.nstar: r=r[0:self.nstar]
        if len(teff)>self.nstar: teff=teff[0:self.nstar]
        if len(logg)>self.nstar: logg=logg[0:self.nstar]
        
        self.distrange=dist
        self.avrange=av
        self.rrange=r
        self.teffrange=teff
        self.loggrange=logg
        self.fehrange=feh
        return
        
    def getfluxsystem(self,dist,av,r,fullfit=False,teff=None,logg=None,feh=None):
        redden=self.la.value*0.
        redden[self.refuv] = self.extuv.extinguish(self.la[self.refuv], Av=av)
        redden[self.refir] = self.extir.extinguish(self.la[self.refir], Av=av)
        
        fx=[]
        for x in range(self.nstar):
            if fullfit:
                flux=self.interp((teff[x],logg[x],feh,self.alpha))
            else:
                flux=self.flux[x]
            if type(r[x]) in [list,tuple,np.ndarray]:
                if len(r[x])==len(self.flux[x]):
                    fxx=[]
                    for xx in range(len(r[x])):
                        fxx.append(((10**(flux[xx])*(r[x][xx]*(u.Rsun.to(u.cm))**2))/(4*np.pi*((dist*u.pc).to(u.cm))**2)*redden).value)
                    fx.append(np.sum(fxx,axis=0))
                else:
                    raise Exception('Number of facets areas has to equal the number of temperatures')
            elif self.area:
                fx.append(((10**(flux)*( r[x]*(u.Rsun.to(u.cm))**2))/(4*np.pi*((dist*u.pc).to(u.cm))**2)*redden).value)
            else:
                fx.append(((10**(flux)*((r[x]*u.Rsun).to(u.cm))**2)/(((dist*u.pc).to(u.cm))**2)*redden).value)
        f=np.sum(fx,axis=0)
        mags=self.fluxtomag(f)
        if len(self.gaia)>0:
            spec=np.log10(self.gaiatomag(f))
        else:
            spec=[]
        return np.log10(f),np.log10(fx),np.log10(mags),spec
    
    def fluxtomag(self,f):
        
        y=[]
        for i in range(len(self.sed)):
            t=Table.read(pkg_resources.resource_filename('SEDFit', self.sed['sed_filter'][i]+'.dat'),format='ascii')
            match=np.interp(t['col1']*u.AA, self.la,f/self.la)
            y.append((np.sum(match*t['col2']/np.sum(t['col2']))*self.sed['la'][i]).value)
        return np.array(y)
    
    def gaiatomag(self,f):
        match=np.interp(self.gaia['la'], self.la,f)
        return match
    
    def makeplot(self,file='',getplot=False):
        plt.rc('font', size=24) 
        fig = plt.figure(figsize=(12,12))
        gs = fig.add_gridspec(2, 1, hspace=0, wspace=0,height_ratios=[3, 1])
        ax=gs.subplots(sharex='col')
        ax[0].errorbar(self.sed["la"]/1e4, self.sed["flux"],
                     yerr=self.sed["eflux"],xerr=self.sed["width"]/1e4,
                     linestyle='',zorder=4,c='black')
        ax[0].set_xscale('log')
        ax[0].scatter(self.sed["la"]/1e4,self.mags,c='#004488',zorder=5,s=80)
        ylims=ax[0].get_ylim()
        
        for i in range(self.nstar):
            ax[0].plot(self.la/1e4,self.fx[i],zorder=0,c='#DDAA33')
        ax[0].plot(self.la/1e4,self.f,zorder=1,c='#BB5566')
        
        
        ylims1=ax[0].get_ylim()
        
        ax[0].set_ylim(ylims[0],ylims1[1])
        ax[0].set_xlim(0.1,20)
        
        ax[1].scatter(self.sed["la"]/1e4,self.sed["flux"]-self.mags,zorder=1,c='#004488',s=80)
        ax[1].errorbar(self.sed["la"]/1e4, self.sed["flux"]-self.mags,
                     yerr=self.sed["eflux"],
                     linestyle='',zorder=0,c='black')
        
        if len(self.gaia)>0:
            ax[0].plot(self.gaia['la']/1e4,self.gaia['flux'],c='black')
            ax[1].plot(self.gaia['la']/1e4,self.gaia['flux']-self.spec,c='black')
        
        ax[1].axhline(y = 0.0, color = 'r')
        ax[1].set_xscale('log')
        ax[1].set_xlim(0.1,20)
        
        ax[1].set_xlabel('$\lambda$ ($\mu$m)')
        ax[0].set_ylabel('log $\lambda F_\lambda$ (erg s$^{-1}$ cm$^{-2}$)')
        ax[1].set_ylabel('Residuals')
        
        ax[1].xaxis.set_major_locator(ticker.FixedLocator([0.1,0.2,0.5,1,2,5,10]))
        ax[1].set_xticklabels(['0.1','0.2','0.5','1','2','5','10'])
        
        if file != '':
            plt.savefig(file, bbox_inches="tight",dpi=150)
        if getplot:
            return ax
        else:
            plt.show()
        
    def getchisq(self,idx=None,gaia=True):
        if idx==None: idx=range(len(self.mags))
        if (gaia) & (len(self.gaia)>0):
            m=np.append(self.mags[idx],self.spec)
            r=np.append(self.sed['flux'][idx],self.gaia['flux'])
            e=np.append(self.sed['eflux'][idx],self.gaia['eflux'])
            l=len(idx)+len(self.gaia)
        else:
            m=self.mags[idx]
            r=self.sed['flux'][idx]
            e=self.sed['eflux'][idx]
            l=len(idx)
            
        return np.sum(((m-r) / e) ** 2)/(l-1)
    
    def getr(self):
        return self.r
    
    def getdist(self):
        return self.dist
    
    def getav(self):
        return self.av
        
    def getteff(self):
        return self.teff
        
    def getlogg(self):
        return self.logg
        
    def getfeh(self):
        return self.feh
    
    
    def fit(self,use_gaia=True,use_mag=[],fitstar=[],teffratio=None,teffratio1=None,teffratio2=None,teffratio_error=None,
                fluxratiolambda=None,fluxratio=None,fluxratio_error=None,
                radiusratio=None,radiusratio_error=None,full=True,fitteff=True,fitlogg=True,fitfeh=True):
        self.full=full
        self.fitteff=fitteff
        self.fitlogg=fitlogg
        self.fitfeh=fitfeh
        self.fluxratiolambda=fluxratiolambda
        self.teffratio=teffratio
        self.teffratio1=teffratio1
        self.teffratio2=teffratio2
        self.radiusratio=radiusratio
        if fluxratio is not None:
            if fluxratiolambda is None:
                raise Exception('Reference wavelength (fluxratiolambda) for the flux ratio is needed')
            try:
                a=fluxratiolambda.unit
            except:
                raise Exception('Reference wavelength for the flux ratio has to have appropriate units')
        
        if len(use_mag)==0: use_mag=range(len(self.sed))
        self.use_mag=use_mag
        if (len(self.gaia)>0) & (use_gaia):
            flux=np.append(self.sed['flux'][self.use_mag],self.gaia['flux'])
            eflux=np.append(self.sed['eflux'][self.use_mag],self.gaia['eflux'])
            self.use_gaia=True
        else:
            flux=self.sed['flux'][self.use_mag]
            eflux=self.sed['eflux'][self.use_mag]
            self.use_gaia=False
        if (fluxratio is not None) & (self.nstar>1):
            flux=np.append(flux,fluxratio)
            if fluxratio_error is not None:
                eflux=np.append(eflux,fluxratio_error)
            else:
                eflux=np.append(eflux,0.1)
        if (teffratio is not None) & (self.nstar>1) & (self.full):
            flux=np.append(flux,teffratio)
            if teffratio_error is not None:
                eflux=np.append(eflux,teffratio_error)
            else:
                eflux=np.append(eflux,0.01)
        if (teffratio1 is not None) & (self.nstar>2) & (self.full):
            flux=np.append(flux,teffratio1)
            if teffratio_error is not None:
                eflux=np.append(eflux,teffratio_error)
            else:
                eflux=np.append(eflux,0.01)
        if (teffratio2 is not None) & (self.nstar>3) & (self.full):
            flux=np.append(flux,teffratio2)
            if teffratio_error is not None:
                eflux=np.append(eflux,teffratio_error)
            else:
                eflux=np.append(eflux,0.01)
        if (radiusratio is not None) & (self.nstar>1):
            flux=np.append(flux,radiusratio)
            if radiusratio_error is not None:
                eflux=np.append(eflux,radiusratio_error)
            else:
                eflux=np.append(eflux,0.1)
            
        if fitstar==[]: fitstar=[1]*self.nstar
        self.fitstar=fitstar
        
        for i in range(len(self.r)):
            if (type(self.r[i]) in [list,tuple,np.ndarray]):
                self.fitstar[i]=0
            
        boundlower,boundupper=[self.distrange[0],self.avrange[0]],[self.distrange[1],self.avrange[1]]
        p0 = self.dist,self.av
        if full:
            
            for i in range(len(self.r)):
                if fitstar[i]>0:
                    p0 = (*p0, self.r[i])
                    boundlower.append(self.rrange[i][0])
                    boundupper.append(self.rrange[i][1])
                    if fitteff:
                        p0 = (*p0,self.teff[i])
                        boundlower.append(self.teffrange[i][0])
                        boundupper.append(self.teffrange[i][1])
                    if fitlogg:
                        p0 = (*p0, self.logg[i])
                        boundlower.append(self.loggrange[i][0])
                        boundupper.append(self.loggrange[i][1])
            if fitfeh:
                p0 = (*p0,self.feh)
                boundlower.append(self.fehrange[0])
                boundupper.append(self.fehrange[1])
        else:
            for i in range(len(self.r)):
                if fitstar[i]>0:
                    p0 = (*p0, self.r[i])
                    boundlower.append(self.rrange[i][0])
                    boundupper.append(self.rrange[i][1])
                    
        
        pars=curve_fit(self.wrap, self.f, flux, p0,
                       sigma=eflux,absolute_sigma=True,
                       bounds=(boundlower,boundupper))
        r,teff,logg=[],[],[]
        j=2
        if full:
            for i in range(self.nstar):
                if self.fitstar[i]>0:
                    r.append(pars[0][j])
                    j=j+1
                    if self.fitteff:
                        teff.append(pars[0][j])
                        j=j+1
                    else:
                        teff.append(self.teff[i])
                    if self.fitlogg:
                        logg.append(pars[0][j])
                        j=j+1
                    else:
                        logg.append(self.logg[i])
                else:
                    r.append(self.r[i])
                    teff.append(self.teff[i])
                    logg.append(self.logg[i])
                    
            if self.fitfeh:
                feh=pars[0][-1]
            else:
                feh=self.feh
            
            self.addguesses(dist=pars[0][0],av=pars[0][1],
                            r=r,teff=teff,logg=logg,feh=feh)
        else:
            for i in range(self.nstar):
                if self.fitstar[i]>0:
                    r.append(pars[0][j])
                    j=j+1
                else:
                    r.append(self.r[i])
            
            self.addguesses(dist=pars[0][0],av=pars[0][1],
                            r=r)
        return pars
        
    def wrap(self,f,dist,av,*argv):
        r,teff,logg=[],[],[]
        j=0
        if self.full:
            if self.fitfeh:
                feh=argv[-1]
            else:
                feh=self.feh
            for i in range(self.nstar):
                if self.fitstar[i]>0:
                    r.append(argv[j])
                    j=j+1
                    if self.fitteff:
                        teff.append(argv[j])
                        j=j+1
                    else:
                        teff.append(self.teff[i])
                    if self.fitlogg:
                        logg.append(argv[j])
                        j=j+1
                    else:
                        logg.append(self.logg[i])
                else:
                    r.append(self.r[i])
                    teff.append(self.teff[i])
                    logg.append(self.logg[i])
            f,fx,m,s=self.getfluxsystem(dist,av,r,teff=teff,logg=logg,feh=feh,fullfit=True)
        else:
            for i in range(self.nstar):
                if self.fitstar[i]>0:
                    r.append(argv[j])
                    j=j+1
                else:
                    r.append(self.r[i])
            f,fx,m,s=self.getfluxsystem(dist,av,r)
            
        m=m[self.use_mag]
        if self.use_gaia: m=np.append(m,s)
            
        if (self.fluxratiolambda is not None) & (self.nstar>1):
            x=np.argmin(np.abs(self.la-self.fluxratiolambda))
            a,b=10**fx[0][x],10**fx[1][x]
            m=np.append(m,b/a)
        if (self.teffratio is not None) & (self.nstar>1) & (self.full):
            m=np.append(m,teff[1]/teff[0])
        if (self.teffratio1 is not None) & (self.nstar>2) & (self.full):
            m=np.append(m,teff[2]/teff[0])
        if (self.teffratio2 is not None) & (self.nstar>3) & (self.full):
            m=np.append(m,teff[3]/teff[1])
        if (self.radiusratio is not None) & (self.nstar>1):
            m=np.append(m,r[1]/r[0])
        return m
