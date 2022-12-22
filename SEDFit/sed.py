import astropy.units as u
import astropy.coordinates as coord
from astropy.io import fits
from astropy.table import Table,vstack
from astropy.modeling import models
from specutils import Spectrum1D
from specutils.manipulation import SplineInterpolatedResampler,FluxConservingResampler
from dust_extinction.averages import GCC09_MWAvg,G21_MWAvg
from dust_extinction.parameter_averages import F99
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
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

class SEDFit:
    def __init__(self,ra,dec,radius,frame='icrs',flux_filename='sed.fits',gaia_filename='gaiaxp.fits',
                 download_flux=False,download_gaia=False,use_gaia=True,**kwargs):
        
        if (' ' in ra) or (':' in dec): raunit=u.hourangle
        else: raunit=u.deg
        c=coord.SkyCoord(ra=ra, dec=dec,unit=(raunit, u.deg),frame='icrs')
        self.c=c
        self.ra=c.ra.value
        self.dec=c.dec.value
        self.radius=radius
        
        self.getmaxreddening(c)
        
        if (not download_flux) & (len(glob.glob(flux_filename))>0):
            self.sed=Table.read(flux_filename)
        else:
            self.downloadflux(**kwargs)
            self.sed.write(flux_filename,overwrite=True)
        
        if use_gaia:
            if (not download_gaia) & (len(glob.glob(gaia_filename))>0):
                self.gaia=Table.read(gaia_filename)
            else:
                try:
                    self.gaia=self.getgaia(c)
                    self.gaia.write(gaia_filename,overwrite=True)
                    print('Gaia XP spectra available')
                except:
                    self.gaia=[]
        else:
            self.gaia=[]
            
        
        self.add_new_grid(**kwargs)
        
        self.extuv=GCC09_MWAvg()
        self.extir=G21_MWAvg()
        
        self.refuv=np.where(self.la<1.1*u.micron)
        self.refir=np.where(self.la>=1.1*u.micron)
        
        self.nstar=1
        self.addrange(dist=[0,1e5],av=[0,20],r=[0,2000])
        self.addguesses(dist=1000.,av=0.,r=[1.],teff=[5000.],logg=[4.],feh=0.,alpha=0.)
        
        
    def getmaxreddening(self,coords):
        a=dustmaps.std_paths.data_dir()
        if len(glob.glob(a+'/sfd'))==0:
            dustmaps.sfd.fetch()
        sfd = SFDQuery()
        self.maxav = sfd(coords)*3.1
        return
        
        
    def downloadflux(self,**kwargs):

        target=str(self.ra)+'%20'+str(self.dec)
        self.sed=Table.read(f"https://vizier.cds.unistra.fr/viz-bin/sed?-c={target}&-c.rs={self.radius}")
        self.sed['la']=self.sed['sed_freq'].to(u.AA, equivalencies=u.spectral())
        a=np.where((self.sed['la']<15*u.micron) & (self.sed['la']>1000*u.AA))[0]
        self.sed=self.sed[a]
        
        self.sed["sed_flux"]=self.sed["sed_flux"].to((u.erg/u.s/(u.cm**2)/u.AA),equivalencies=u.spectral_density(self.sed['la'].data*u.AA))
        self.sed["sed_eflux"]=self.sed["sed_eflux"].to((u.erg/u.s/(u.cm**2)/u.AA),equivalencies=u.spectral_density(self.sed['la'].data*u.AA))
        a=np.where((np.isnan(self.sed["sed_eflux"])==True) | (self.sed["sed_eflux"]/self.sed["sed_flux"]<0.02))[0]
        self.sed["sed_eflux"][a]=self.sed["sed_flux"][a]*0.02
        
        self.sed['eflux']=self.sed["sed_eflux"]/self.sed["sed_flux"]/np.log(10)
        self.sed['flux'] =np.log10(self.sed["sed_flux"]*self.sed['la'])
        
        self.definefilter(**kwargs)
        a=np.argsort(self.sed['la'])
        self.sed=self.sed[a]
        self.sed.keep_columns(['sed_filter','la','width','flux','eflux'])
        
        
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
            idx.extend(self.selectflux(filters,width))
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
        if checkpos:
            b=np.min(dist)
            a=np.where(dist==b)[0]
            idx=idx[a]
        return idx
    def getgaia(self,c):
        Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
        j = Gaia.cone_search_async(c, self.radius*u.arcsec,verbose=False)
        sid=j.get_results()['source_id']
        
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

    def load_phoenix_sed(self,name,laname):
        spline = SplineInterpolatedResampler()
        hdul = fits.open(name)
        spec = hdul[0].data*u.erg/u.s/(u.cm**3)
        spec=spec.to(u.erg/u.s/(u.cm**2)/u.AA)
        header=hdul[0].header
        hdul.close()
        hdul = fits.open(laname)
        la=hdul[0].data*u.AA
        lsf=self.lsf_rotate(1,100)
        spec=np.convolve(spec,lsf)/100
        a=np.where((la>1000*u.AA))[0]
        hdul.close()
        sp = Spectrum1D(spectral_axis=la[a], flux=spec[a])
        la=la[a[np.arange(0,len(a),200)]]
        sp=spline(sp,la)
        spec=sp.flux
        
        bb = models.BlackBody(temperature=header['PHXTEFF']*u.K)
        l=np.arange(55000,320000,5000)*u.AA
        b=((bb(l)*u.sr).to(u.erg/u.s/(u.cm**2)/u.AA,equivalencies=u.spectral_density(l)))
        b=b/b[0]*spec[-1]
        la=np.append(la,l[1:])
        spec=np.append(spec,b[1:])
        
        return [np.round(np.log10((spec*la).value),3)],la,[header['PHXTEFF']],[header['PHXLOGG']],[header['PHXM_H']],[header['PHXALPHA']]
    
    def load_atlas_sed(self,name):
        n=name.split('/')[-1].replace('m','-').replace('p','+')
        t=Table.read(name,format='ascii')
        la=t['col1']*u.AA
        spec=t['col2']*u.erg/u.s/(u.cm**2)/u.AA
        
        if n[4]=='t':
            teff=np.float(n.split('t')[1].split('g')[0])
            logg=np.float(n.split('g')[1].split('k')[0])/10
            feh=np.float(n[1:4])/10
            alpha=0
            
        else:
            teff=np.float(n.split('t')[1].split('g')[0])
            logg=np.float(n.split('g')[1].split('k')[0])/10
            feh=np.float(n[1:4])/10
            alpha=0.4
        
        a=np.where((la<10*u.micron) & (la>912*u.AA))[0]
        spec,la=spec[a],la[a]
        l=np.arange(100000,320000,5000)*u.AA
        
        bb = models.BlackBody(temperature=teff*u.K)
        b=((bb(l)*u.sr).to(u.erg/u.s/(u.cm**2)/u.AA,equivalencies=u.spectral_density(l)))
        b=b/b[0]*spec[-1]
        la=np.append(la,l[1:])
        spec=np.append(spec,b[1:])
        
        return [np.round(np.log10((spec*la).value),3)],la,[teff],[logg],[feh],[alpha]


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

    def add_new_grid(self,grid_path='',grid_pickle='',grid_type='coelho',loadgrid=True,laname=None,**kwargs):
        
        if grid_pickle=='':
            grid_pickle='sed_'+grid_type+'.p'
            fullpath=pkg_resources.resource_filename('SEDFit', grid_pickle)
        else:
            fullpath=grid_pickle
        if loadgrid & len(glob.glob(fullpath))>0:
            flux,t,la=pickle.load(open(fullpath,'rb'))
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
            for i in tqdm(path):
                f,la,teff,logg,feh,alpha=self.load_sed(grid_type,i,laname)
                for j in range(len(f)):
                    flux.append(f[j])
                    t.add_row((teff[j],logg[j],feh[j],alpha[j]))
            pickle.dump( [flux,t,la], open(grid_pickle, "wb" ) )
        self.la=la
        self.t=t
        self.flux=flux
        return
    def load_sed(self,grid_type,i,laname):
        if grid_type=='coelho':
            return self.load_coelho_sed(i)
        elif grid_type=='kurucz':
            return self.load_kurucz_sed(i)
        elif grid_type=='atlas9':
            return self.load_atlas_sed(i)
        elif grid_type=='phoenix':
            if laname is None:
                raise Exception("Please provide wavelength array path (laname)")
            return self.load_phoenix_sed(i,laname)
        else:
            raise Exception("Unsupported grid type, available options are 'coelho', 'kurucz', and 'phoenix'")
        
        
    def getclosestpars(self,teff,logg,feh,alpha):
        q=np.abs(self.t['alpha']-alpha)
        a=np.where(q==np.min(q))[0]
        z=np.abs(self.t['feh'][a]-feh)
        a=a[np.where(z==np.min(z))[0]]
        x=np.abs(self.t['teff'][a]-teff)
        a=a[np.where(x==np.min(x))[0]]
        y=np.abs(self.t['logg'][a]-logg)
        a=a[np.argmin(y)]
        return a
        
    
    def addguesses(self,dist=None,av=None,r=None,teff=None,logg=None,feh=None,alpha=None,area=False):
        if dist is None: dist=self.dist
        if av is None: av=self.av
        if r is None: r=self.r
        if teff is None: teff=self.teff
        if logg is None: logg=self.logg
        if feh is None: feh=self.feh
        if alpha is None: alpha=self.alpha
        
        if not type(r) in [list,tuple,np.ndarray]: r=[r]
        if not type(teff) in [list,tuple,np.ndarray]: teff=[teff]
        if not type(logg) in [list,tuple,np.ndarray]: logg=[logg]
            

        self.nstar=np.max([len(r),len(teff),len(logg)])
        
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
        
        
        self.index=[]
        for i in range(len(self.r)):
            self.index.append(self.getclosestpars(self.teff[i],self.logg[i],self.feh,self.alpha))
        self.f,self.fx,self.mags,self.spec=self.getfluxsystem(self.dist,self.av,self.r,area=area)
        self.addrange()
        return
    

    def addrange(self,dist=None,av=None,r=None):
        if dist is None: dist=[self.boundlower[0],self.boundupper[0]]
        if av is None: av=[self.boundlower[1],self.boundupper[1]]
        if (r is None) & (self.nstar==1): r=[self.boundlower[2],self.boundupper[2]]
        if (r is None) & (self.nstar==2): 
            if len(self.boundlower)==3:
                self.boundlower.append(self.boundlower[2])
                self.boundupper.append(self.boundupper[2])
            r=[[self.boundlower[2],self.boundupper[2]],[self.boundlower[3],self.boundupper[3]]]
        if (r is None) & (self.nstar==3): 
            
            if len(self.boundlower)==3:
                self.boundlower.extend([self.boundlower[2],self.boundlower[2]])
                self.boundupper.extend([self.boundupper[2],self.boundupper[2]])
            r=[[self.boundlower[2],self.boundupper[2]],[self.boundlower[3],self.boundupper[3]],
               [self.boundlower[4],self.boundupper[4]]]
            
        if av[1]>self.maxav:
            print('Maximum Av along the line of sight is '+str(np.round(self.maxav,3)))
            av[1]=self.maxav
        
        if self.nstar==1:
            boundlower=[dist[0],av[0],r[0]]
            boundupper=[dist[1],av[1],r[1]]
        if self.nstar==2:
            if np.array(r).size==2: r=[r]*2
            boundlower=[dist[0],av[0],r[0][0],r[1][0]]
            boundupper=[dist[1],av[1],r[0][1],r[1][1]]
        if self.nstar==3:
            if np.array(r).size==2: r=[r]*3
            if np.array(r).size==4: raise Exception('Cannot parse the constraint for 3 stars')
            boundlower=[dist[0],av[0],r[0][0],r[1][0],r[2][0]]
            boundupper=[dist[1],av[1],r[0][1],r[1][1],r[2][1]]
        if self.nstar>3:
            boundupper=[0,0,0,0,0]
            boundlower=[0,0,0,0,0]
            #raise Exception("Systems with 3+ stars are not supported yet")
            
        self.boundlower=boundlower
        self.boundupper=boundupper
        return
        
    def getfluxsystem(self,dist,av,r,area=False):
        redden=self.la.value*0.
        redden[self.refuv] = self.extuv.extinguish(self.la[self.refuv], Av=av)
        redden[self.refir] = self.extir.extinguish(self.la[self.refir], Av=av)
        
        fx=[]
        for x in range(self.nstar):
            if area:
                fx.append(((10**(self.flux[self.index[x]])*(r[x]*(u.Rsun.to(u.cm))**2))/(4*np.pi*((dist*u.pc).to(u.cm))**2)*redden).value)
            else:
                fx.append(((10**(self.flux[self.index[x]])*((r[x]*u.Rsun).to(u.cm))**2)/(((dist*u.pc).to(u.cm))**2)*redden).value)
        f=np.sum(fx,axis=0)
        mags=self.fluxtomag(f)
        if len(self.gaia)>0:
            spec=np.log10(self.gaiatomag(f))
        else:
            spec=[]
        return np.log10(f),np.log10(fx),np.log10(mags),spec
    
    def lsf_rotate(self,deltav,vsini,epsilon=None):
        # based on the IDL routine LSF_ROTATE.PRO
        if epsilon == None:
            epsilon = 0.6
        e1 = 2.0*(1.0 - epsilon)
        e2 = np.pi*epsilon/2.0
        e3 = np.pi*(1.0 - epsilon/3.0)
        npts = np.ceil(2*vsini/deltav)
        if npts % 2 == 0: npts += 1
        nwid = np.floor(npts/2)
        x = np.arange(npts) - nwid
        x = x*deltav/vsini  
        x1 = np.abs(1.0 - x**2)
        return (e1*np.sqrt(x1) + e2*x1)/e3
    
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
    
    def funcsingle(self,f,dist,av,r1):
        f,fx,m,s=self.getfluxsystem(dist,av,[r1])
        m=m[self.use_mag]
        if self.use_gaia: return np.append(m,s)
        else: return m
    def fitsingle(self):
        p0 = self.dist,self.av,self.r[0]
        pars=curve_fit(self.funcsingle, self.f, flux, p0,
                       sigma=eflux,absolute_sigma=True,
                       bounds=(self.boundlower,self.boundupper))
        self.addguesses(dist=pars[0][0],av=pars[0][1],
                        r=[pars[0][2]])
        #print(pars[0])
        return pars

    def funcbinary_r(self,f,dist,av,r1,r2):
        f,fx,m,s=self.getfluxsystem(dist,av,[r1,r1-r2])
        m=m[self.use_mag]
        if self.use_gaia: return np.append(m,s)
        else: return m
    def fitbinary_r(self,flux,eflux):
        p0 = self.dist,self.av,self.r[0],self.r[0]-self.r[1]
        
        lower=self.boundlower
        upper=self.boundupper
        lower[3]=np.max([0,self.r[0]-self.boundupper[3]])
        upper[3]=self.r[0]-self.boundlower[3]
        
        pars=curve_fit(self.funcbinary_r, self.f, flux, p0,
                       sigma=eflux,absolute_sigma=True,
                       bounds=(lower,upper))
        self.addguesses(dist=pars[0][0],av=pars[0][1],
                        r=[pars[0][2],pars[0][2]-pars[0][3]])
        return pars
    
    def funcbinary(self,f,dist,av,r1,r2):
        f,fx,m,s=self.getfluxsystem(dist,av,[r1,r2])
        m=m[self.use_mag]
        if self.use_gaia: return np.append(m,s)
        else: return m
    def fitbinary(self,flux,eflux):
        p0 = self.dist,self.av,self.r[0],self.r[1]
        
        pars=curve_fit(self.funcbinary, self.f, flux, p0,
                       sigma=eflux,absolute_sigma=True,
                       bounds=(self.boundlower,self.boundupper))
        self.addguesses(dist=pars[0][0],av=pars[0][1],
                        r=[pars[0][2],pars[0][3]])
        return pars
    
    def functriple(self,f,dist,av,r1,r2,r3):
        f,fx,m,s=self.getfluxsystem(dist,av,[r1,r2,r3])
        m=m[self.use_mag]
        if self.use_gaia: return np.append(m,s)
        else: return m
    def fittriple(self,flux,eflux):
        p0 = self.dist,self.av,self.r[0],self.r[1],self.r[2]
        
        pars=curve_fit(self.functriple, self.f, flux, p0,
                       sigma=eflux,absolute_sigma=True,
                       bounds=(self.boundlower,self.boundupper))
        self.addguesses(dist=pars[0][0],av=pars[0][1],
                        r=[pars[0][2],pars[0][3],pars[0][4]])
        return pars
    
    def fit(self,order_radii=False,use_gaia=True,use_mag=[]):
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
        
        if self.nstar==1: return self.fitsingle(flux,eflux)
        if (self.nstar==2) & (order_radii): return self.fitbinary_r(flux,eflux)
        if (self.nstar==2) & (not order_radii): return self.fitbinary(flux,eflux)
        if self.nstar==3: return self.fittriple(flux,eflux)
        if self.nstar>3:
            raise Exception('Systems with 4+ stars are not supported yet')
    
    
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
        
    def getchisq(self):
        
        return np.sum(((self.mags-self.sed['flux']) / self.sed["eflux"]) ** 2)/(len(self.sed)+len(self.gaia)-1)
    
    def getr(self):
        return self.r
    
    def getdist(self):
        return self.dist
    
    def getav(self):
        return self.av