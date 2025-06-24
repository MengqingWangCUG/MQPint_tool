import os
import sys
import copy
import matplotlib
from past.utils import old_div
import numpy as np
import pandas as pd
import scipy as sp
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt
import ultraplot as upt
from pmagpy import mqplot as mpt
from matplotlib.pyplot import MultipleLocator
from pmagpy import pmag
from pmagpy import ipmag
from pmagpy import find_pmag_dir
from pmagpy import contribution_builder as cb
from pmagpy import convert_2_magic as convert
import math
import cartopy
from cartopy import crs as ccrs
from numpy import linalg
import matplotlib.ticker as ticker
#import textalloc as ta

import inspect
#import textalloc as ta
from adjustText import adjust_text
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams.update({'font.size': 8})
# matplotlib.rcParams['lines.markeredgewidth'] = 0.5  
# matplotlib.rcParams['lines.linewidth'] = 0.5 
# matplotlib.rcParams['lines.markersize'] = 4
from scipy import interpolate,spatial
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import norm
colors = ['#DE3F2A',
    '#F78404',
    '#F9D83F',
    '#6A961F',
    '#15A99E',
    '#0487E2',
    '#804595',
    '#89678B'
]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors)

def c2np(data):
    if isinstance(data, np.ndarray):
        return data
    elif isinstance(data, pd.DataFrame):
        return data.to_numpy()
    elif isinstance(data, (list, tuple)):
        return np.array(data)
    else:
        raise TypeError(f"Unsupported data type: {type(data)}. Input must be tuple, list, pandas DataFrame, or numpy array.")
    
 

class Specimen:
    def __init__(self,generaldf,momentunit='$Am^2$',redodf=None,Azimuth=None,Dip=None):
        generaldf['moment']=generaldf['moment'].astype(float)/1000 # convert to Am^2
        self.generaldf=generaldf
        self.specimen=generaldf['specimen'].unique()
        self.dec_s=generaldf['dec_s']
        self.inc_s=generaldf['inc_s']
        self.moment=generaldf['moment'] 
        self.NRMmoment=self.moment[0]
        self.momentunit=momentunit
        self.treatment=generaldf['treatment']
        self.treatment_type=generaldf['treatment_type']
        self.normmoment=pd.concat([self.treatment,self.moment/self.NRMmoment],axis=1)
        self.treatment_unit=pd.DataFrame({'unite': ['mT'] * len(self.treatment_type)})
        self.treatment_unit.loc[self.treatment_type=='T','unite']='°C'
        generaldf.loc[:,['dec_g','inc_g']]=erm(generaldf['dec_s'],generaldf['inc_s'],Azimuth,90-Dip)
        self.treatmenttext = pd.DataFrame({'treatmtext': self.treatment.astype(str) + self.treatment_unit['unite'].astype(str)})
        if self.treatment[0]==0:
            self.treatmenttext.loc[0,'treatmtext']='NRM'
        else:
            pass
        self.dec_g=generaldf['dec_g']
        self.inc_g=generaldf['inc_g']

        if redodf is None:
            pass
        else:
            self.change_temp(redodf=redodf)
    def change_temp(self,redodf):
        self.selectbool,self.selectlow,self.selecthigh,self.pca_g,self.fisher_g,self.pca_t,self.fisher_t,self.pca_s,self.fisher_s=redo_pca(redodf,self.generaldf,self.treatment_type,self.specimen)



class Sample:
    def __init__(self, name=None, Azimuth=None, Dip=None, data=None):
        self.sample_name = name
        self.Azimuth = Azimuth
        self.Dip = Dip
        self.data = data
        self.specimens = {}
    def __getitem__(self, keys):
        if isinstance(keys, (tuple, list)):
            if len(keys) == 1:
                specimen_name = keys[0]
                return self.specimens[specimen_name]
            else:
                new_sample = Sample(self.sample_name, self.data)
                for key in keys:
                    if key in self.specimens:
                        new_sample.specimens[key] = self.specimens[key]
                return new_sample
        elif isinstance(keys, str):
            return self.specimens[keys]
        else:
            raise TypeError("Invalid key type")

    def __setitem__(self, keys, value):
        if isinstance(keys, (tuple, list)):
            if len(keys) == 1:
                specimen_name = keys[0]
                self.specimens[specimen_name] = value
            else:
                for key in keys:
                    if key in self.specimens:
                        self.specimens[key] = value
        elif isinstance(keys, str):
            self.specimens[keys] = value
        else:
            raise TypeError("Invalid key type")
            
    def __delitem__(self, keys):
        if isinstance(keys, (tuple, list)):
            if len(keys) == 1:
                specimen_name = keys[0]
                if specimen_name in self.specimens:
                    del self.specimens[specimen_name]
                else:
                    raise KeyError(f"Specimen '{specimen_name}' not found")
            else:
                for key in keys:
                    if key in self.specimens:
                        del self.specimens[key]
        elif isinstance(keys, str):
            if keys in self.specimens:
                del self.specimens[keys]
            else:
                raise KeyError(f"Specimen '{keys}' not found")
        else:
            raise TypeError("Invalid key type")
            
    def sampleslist(self):
        return [self.sample_name]
        
    def specimenslist(self):
        return list(self.specimens.keys())


class Site:
    def __init__(self, name=None, lat=None, lon=None, data=None):
        self.site_name = name
        self.lat = lat.values[0]
        self.lon = lon.values[0]
        self.data = data
        self.samples = {}
        self.site_spe_s_PCA = pd.DataFrame([])
        self.site_spe_g_PCA = pd.DataFrame([])
        self.site_spe_t_PCA = pd.DataFrame([])
        self.site_spe_s_Fisher = pd.DataFrame([])
        self.site_spe_g_Fisher = pd.DataFrame([])
        self.site_spe_t_Fisher = pd.DataFrame([])
        self.sitePCA_s = pd.DataFrame([])
        self.sitePCA_g = pd.DataFrame([])
        self.sitePCA_t = pd.DataFrame([])
        self.siteFisher_s = pd.DataFrame([])
        self.siteFisher_g = pd.DataFrame([])
        self.siteFisher_t = pd.DataFrame([])
    def __getitem__(self, keys):
        if isinstance(keys, (tuple, list)):
            if len(keys) == 1:
                sample_name = keys[0]
                return self.samples[sample_name]
            else:
                new_site = Site(self.site_name, self.data)
                for key in keys:
                    if key in self.samples:
                        new_site.samples[key] = self.samples[key]
                return new_site
        elif isinstance(keys, str):
            return self.samples[keys]
        else:
            raise TypeError("Invalid key type")

    def __setitem__(self, keys, value):
        if isinstance(keys, (tuple, list)):
            if len(keys) == 1:
                sample_name = keys[0]
                self.samples[sample_name] = value
            else:
                for key in keys:
                    if key in self.samples:
                        self.samples[key] = value
        elif isinstance(keys, str):
            self.samples[keys] = value
        else:
            raise TypeError("Invalid key type")
            
    def __delitem__(self, keys):
        if isinstance(keys, (tuple, list)):
            if len(keys) == 1:
                sample_name = keys[0]
                if sample_name in self.samples:
                    del self.samples[sample_name]
                else:
                    raise KeyError(f"Sample '{sample_name}' not found")
            else:
                for key in keys:
                    if key in self.samples:
                        del self.samples[key]
        elif isinstance(keys, str):
            if keys in self.samples:
                del self.samples[keys]
            else:
                raise KeyError(f"Sample '{keys}' not found")
        else:
            raise TypeError("Invalid key type")
            
    def siteslist(self):
        return [self.site_name]
        
    def sampleslist(self):
        return list(self.samples.keys())
        
    def specimenslist(self):
        all_specimens = []
        for sample in self.samples.values():
            all_specimens.extend(sample.specimenslist())
        return all_specimens


class mqdirdata:
    def __init__(self, generalpath, inforpath, redopath=None, Sitelen=3, Samplelen=4,momentunit='$Am^2$',specimenselect=True,treatmentNj=4,specimenMADj=10,specimenDangJ=5):
        self.site = {}
        self.momentchange=pd.DataFrame([])
        self.allsite_PCA_fisher_s=pd.DataFrame([])
        self.allsite_PCA_fisher_g=pd.DataFrame([])
        self.allsite_PCA_fisher_t=pd.DataFrame([])
        self.allsite_fisher_fisher_s=pd.DataFrame([])
        self.allsite_fisher_fisher_g=pd.DataFrame([])
        self.allsite_fisher_fisher_t=pd.DataFrame([])
        self.allsiteinfor=pd.DataFrame([])
        generaldf=pd.read_excel(generalpath)
        allinfordf=pd.read_excel(inforpath)
        specimenlist=pd.Series(generaldf['specimen'].unique())
        samplelist=pd.Series(pd.Series(specimenlist).str[:Samplelen].unique())
        sitelist=pd.Series(pd.Series(samplelist).str[:Sitelen].unique())
        
        redodf=pd.read_csv(redopath,delim_whitespace=True,header=None)
        redodf[0]=redodf[0].str.replace('current_', '')
        try:
            redodf=pd.read_csv(redopath,delim_whitespace=True,header=None)
            redodf[0]=redodf[0].str.replace('current_', '')

        except:
            redodf=None
        for site in sitelist:
            siteinfordf=allinfordf[allinfordf['Site']==site].drop_duplicates(subset='Site')
            self.site[site]=Site(name=site, lat=siteinfordf.Lat, lon=siteinfordf.Lon, data=siteinfordf)
            Single_samplelist=samplelist[samplelist.str[:Sitelen]==site]
            for sample in Single_samplelist:
                sampleinfor=allinfordf[allinfordf['Sample']==sample].drop_duplicates(subset='Sample')
                Single_specimenlist=specimenlist[specimenlist.str[:Samplelen]==sample]

                self.site[site].samples[sample]=Sample(name=sample,Azimuth=sampleinfor['Azimuth'],Dip=sampleinfor['Dip'])
                Azimuth=self.site[site].samples[sample].Azimuth
                Dip=self.site[site].samples[sample].Dip
                for specimen in Single_specimenlist:
                   
                    single_generaldf=generaldf[generaldf['specimen']==specimen].reset_index(drop=True)
                    if redodf is not None:
                        redos=redodf.loc[redodf[0]==specimen,:].reset_index(drop=True)
                    else:
                        redos=None

                    self.site[site].samples[sample].specimens[specimen]=Specimen(single_generaldf,momentunit=momentunit,redodf=redos,Azimuth=Azimuth,Dip=Dip)
                    Smomentchange=self.site[site].samples[sample].specimens[specimen].normmoment.copy()
                    Smomentchange['specimen']=specimen
                    Smomentchange['sample']= sample
                    Smomentchange['site']=site
                    self.momentchange=pd.concat([self.momentchange,Smomentchange],axis=0)
                    try:
                        self.site[site].site_spe_s_PCA=pd.concat([self.site[site].site_spe_s_PCA,self.site[site].samples[sample].specimens[specimen].pca_s],axis=0).reset_index(drop=True)
                        self.site[site].site_spe_s_Fisher=pd.concat([self.site[site].site_spe_s_Fisher,self.site[site].samples[sample].specimens[specimen].fisher_s],axis=0).reset_index(drop=True)
                    except:
                        pass
                    
                    try:
                        self.site[site].site_spe_g_PCA=pd.concat([self.site[site].site_spe_g_PCA,self.site[site].samples[sample].specimens[specimen].pca_g],axis=0).reset_index(drop=True)
                        self.site[site].site_spe_g_Fisher=pd.concat([self.site[site].site_spe_g_Fisher,self.site[site].samples[sample].specimens[specimen].fisher_g],axis=0).reset_index(drop=True)
                    except:
                        pass
                    try:
                        self.site[site].site_spe_t_PCA=pd.concat([self.site[site].site_spe_t_PCA,self.site[site].samples[sample].specimens[specimen].pca_t],axis=0).reset_index(drop=True)
                        self.site[site].site_spe_t_Fisher=pd.concat([self.site[site].site_spe_t_Fisher,self.site[site].samples[sample].specimens[specimen].fisher_t],axis=0).reset_index(drop=True)
                    except:
                        pass

            if specimenselect:
                self.site[site].site_spe_s_PCA=self.site[site].site_spe_s_PCA.loc[(self.site[site].site_spe_s_PCA['specimen_n']>=treatmentNj)&(self.site[site].site_spe_s_PCA['specimen_mad']<=specimenMADj),:].reset_index(drop=True)
                self.site[site].site_spe_g_PCA=self.site[site].site_spe_g_PCA.loc[(self.site[site].site_spe_g_PCA['specimen_n']>=treatmentNj)&(self.site[site].site_spe_g_PCA['specimen_mad']<=specimenMADj),:].reset_index(drop=True)
                try:
                    self.site[site].site_spe_t_PCA=self.site[site].site_spe_t_PCA.loc[(self.site[site].site_spe_t_PCA['specimen_n']>=treatmentNj)&(self.site[site].site_spe_t_PCA['specimen_mad']<=specimenMADj),:].reset_index(drop=True)
                except:
                    pass
                try:
                    self.site[site].site_spe_s_PCA=self.site[site].site_spe_s_PCA.loc[self.site[site].site_spe_s_PCA['specimen_dang']<=specimenDangJ,:].reset_index(drop=True)
                    self.site[site].site_spe_g_PCA=self.site[site].site_spe_g_PCA.loc[self.site[site].site_spe_g_PCA['specimen_dang']<=specimenDangJ,:].reset_index(drop=True)
                except:
                    pass               
                try:
                    self.site[site].site_spe_t_PCA=self.site[site].site_spe_t_PCA.loc[self.site[site].site_spe_t_PCA['specimen_dang']<=specimenDangJ,:].reset_index(drop=True)
                except:
                    pass
            else:
                pass
            try:
                self.site[site].site_spe_s_PCA=self.site[site].site_spe_s_PCA.reset_index(drop=True)
                self.site[site].site_spe_s_Fisher=self.site[site].site_spe_s_Fisher.reset_index(drop=True)
                allDIlistPCA=self.site[site].site_spe_s_PCA[['specimen_dec','specimen_inc']].values.tolist()
                allDIlistFisher=self.site[site].site_spe_s_Fisher[['specimen_dec','specimen_inc']].values.tolist()
                self.site[site].sitePCA_s=pd.DataFrame(fisher_mean(di_block=allDIlistPCA),index=[0]) 
                self.site[site].siteFisher_s=pd.DataFrame(fisher_mean(di_block=allDIlistFisher),index=[0])
                self.allsite_PCA_fisher_s=pd.concat([self.allsite_PCA_fisher_s,self.site[site].sitePCA_s],axis=0).reset_index(drop=True)
                self.allsite_fisher_fisher_s=pd.concat([self.allsite_fisher_fisher_s,self.site[site].siteFisher_s],axis=0).reset_index(drop=True)
            except:
                pass
            try:
                self.site[site].site_spe_g_PCA=self.site[site].site_spe_g_PCA.reset_index(drop=True)
                self.site[site].site_spe_g_Fisher=self.site[site].site_spe_g_Fisher.reset_index(drop=True)
                allDIlistPCA=self.site[site].site_spe_g_PCA[['specimen_dec','specimen_inc']].values.tolist()
                allDIlistFisher=self.site[site].site_spe_g_Fisher[['specimen_dec','specimen_inc']].values.tolist()
                self.site[site].sitePCA_g=pd.DataFrame(fisher_mean(di_block=allDIlistPCA),index=[0])
                self.site[site].siteFisher_g=pd.DataFrame(fisher_mean(di_block=allDIlistFisher),index=[0])
                self.allsite_PCA_fisher_g=pd.concat([self.allsite_PCA_fisher_g,self.site[site].sitePCA_g],axis=0).reset_index(drop=True)
                self.allsite_fisher_fisher_g=pd.concat([self.allsite_fisher_fisher_g,self.site[site].siteFisher_g],axis=0).reset_index(drop=True)
            except:
                pass
            try:
                self.site[site].site_spe_t_PCA=self.site[site].site_spe_t_PCA.reset_index(drop=True)
                self.site[site].site_spe_t_Fisher=self.site[site].site_spe_t_Fisher.reset_index(drop=True)
                allDIlistPCA=self.site[site].site_spe_t_PCA[['specimen_dec','specimen_inc']].values.tolist()
                allDIlistFisher=self.site[site].site_spe_t_Fisher[['specimen_dec','specimen_inc']].values.tolist()
                self.site[site].sitePCA_t=pd.DataFrame(fisher_mean(di_block=allDIlistPCA),index=[0]) 
                self.site[site].siteFisher_t=pd.DataFrame(fisher_mean(di_block=allDIlistFisher),index=[0]) 
                self.allsite_PCA_fisher_t=pd.concat([self.allsite_PCA_fisher_t,self.site[site].sitePCA_t],axis=0).reset_index(drop=True)
                self.allsite_fisher_fisher_t=pd.concat([self.allsite_fisher_fisher_t,self.site[site].siteFisher_t],axis=0).reset_index(drop=True)
            except:
                pass  
            try:    
                allDIsite_PCA_s=self.allsite_PCA_fisher_s[['dec','inc']].values.tolist()
                allDIsite_Fisher_s=self.allsite_fisher_fisher_s[['dec','inc']].values.tolist()
                self.locationfisher_s=pd.DataFrame(fisher_mean(di_block=allDIsite_Fisher_s),index=[0])
                self.locationPCA_s =pd.DataFrame(fisher_mean(di_block=allDIsite_PCA_s),index=[0]) 
            except:
                pass
            try:
                allDIsite_PCA_g=self.allsite_PCA_fisher_g[['dec','inc']].values.tolist()
                allDIsite_Fisher_g=self.allsite_fisher_fisher_g[['dec','inc']].values.tolist()
                self.locationfisher_g=pd.DataFrame(fisher_mean(di_block=allDIsite_Fisher_g),index=[0])
                self.locationPCA_g =pd.DataFrame(fisher_mean(di_block=allDIsite_PCA_g),index=[0])
            except:
                pass
            try:
                allDIsite_PCA_t=self.allsite_PCA_fisher_t[['dec','inc']].values.tolist()
                allDIsite_Fisher_t=self.allsite_fisher_fisher_t[['dec','inc']].values.tolist()
                self.locationfisher_t=pd.DataFrame(fisher_mean(di_block=allDIsite_Fisher_t),index=[0])
                self.locationPCA_t =pd.DataFrame(fisher_mean(di_block=allDIsite_PCA_t),index=[0])
            except:
                pass
        self.calcuVGP()
        self.LocationVGPmean()
#-----------------------------------------------------------------------
        self._is_selection = False
        self._selected_sites = None
    def calcuVGP(self,
        Coordinates='g',
        site_nj=5,
        site_a95j=8,
        site_kj=50, 

        ):
        siteslist=self.siteslist()
        allsiteinfor=pd.DataFrame([])
        for i, site in enumerate(siteslist):
            try:
                if Coordinates == 'g':
                    #print(self.site[site].sitePCA_g)
                    Sitedec=self.site[site].sitePCA_g.dec.iloc[0]
                    Siteinc=self.site[site].sitePCA_g.inc.iloc[0]
                    Sitealpha95=self.site[site].sitePCA_g.alpha95.iloc[0]
                    Sitek=self.site[site].sitePCA_g.k.iloc[0]
                    Siten=self.site[site].sitePCA_g.n.iloc[0]
                elif Coordinates == 't':
                    Sitedec=self.site[site].sitePCA_t.dec.iloc[0]
                    Siteinc=self.site[site].sitePCA_t.inc.iloc[0]
                    Sitealpha95=self.site[site].sitePCA_t.alpha95.iloc[0]
                    Sitek=self.site[site].sitePCA_t.k.iloc[0]
                    Siten=self.site[site].sitePCA_t.n.iloc[0]
                lat=self.site[site].lat
                
                lon=self.site[site].lon

                DIALL_list=[[Sitedec,Siteinc,Sitealpha95,lat,lon]]
                
                lonlat_list=dia_vgp(DIALL_list)
                SiteVGP=pd.DataFrame([])
                SiteVGP.loc[i,'site']=site
                SiteVGP.loc[i,'Lat']=lat
                SiteVGP.loc[i,'k']=Sitek
                SiteVGP.loc[i,'n']=Siten
                

                SiteVGP.loc[i,'Plat']=lonlat_list.loc[0,'Plat']
                self.site[site].Plat=lonlat_list.loc[0,'Plat']
                SiteVGP.loc[i,'Plon']=lonlat_list.loc[0,'Plon']
                self.site[site].Plon=lonlat_list.loc[0,'Plon']
                SiteVGP.loc[i,'dp']=lonlat_list.loc[0,'dp']
                self.site[site].dp=lonlat_list.loc[0,'dp']
                SiteVGP.loc[i,'dm']=lonlat_list.loc[0,'dm']
                self.site[site].dm=lonlat_list.loc[0,'dm']
                SiteVGP.loc[i,'SPA95']=np.sqrt(lonlat_list.loc[0,'dp']*lonlat_list.loc[0,'dm'])
                self.site[site].SPA95=SiteVGP.loc[i,'SPA95']
                if SiteVGP.loc[i,'Plat']>0:
                    SiteVGP.loc[i,'polarity']='N'
                    self.site[site].polarity='N'
                    self.site[site].flipPlat=lonlat_list.loc[0,'Plat']
                    self.site[site].flipPlon=lonlat_list.loc[0,'Plon']
                    SiteVGP.loc[:,['fPlon','fPlat']]=SiteVGP.loc[:,['Plon','Plat']].values
                else:
                    SiteVGP.loc[i,'polarity']='R'
                    self.site[site].polarity='R'
                    flippole=flip(SiteVGP.loc[:,['Plon','Plat']].values.tolist(),combine=True)
                    self.site[site].flipPlat=flippole[0][1]
                    self.site[site].flipPlon=flippole[0][0]
                    SiteVGP.loc[:,['fPlon','fPlat']]=flippole
                goodsitebool=(Sitealpha95<=site_a95j)&(Sitek>=site_kj)&(Siten>=site_nj)
                SiteVGP.loc[:,'goodbool']=goodsitebool
                self.allsiteinfor=pd.concat([self.allsiteinfor,SiteVGP],axis=0).reset_index(drop=True)
            except:
                print('Site',site,'success S data number',len(self.site[site].site_spe_s_PCA))
                print('Site',site,'success G data number',len(self.site[site].site_spe_g_PCA))
                print('Site',site,'success T data number',len(self.site[site].site_spe_t_PCA))

            
        #return self.allsiteinfor
    def LocationVGPmean(self, 
                        kappa=0, n=0, spin=False,boot=True, mm97=True, nb=1000,select=True
                            ):
        """
        coordinates: 'g' for geographic, 't' for tilt corrected form calcuVGP
        vgp_df : Pandas Dataframe with columns
            REQUIRED:
            vgp_lat :  VGP latitude
            ONLY REQUIRED for MM97 correction:
            dir_k : Fisher kappa estimate
            dir_n_samples : number of samples per site
            lat : latitude of the site
            mm97 : if True, will do the correction for within site scatter
            OPTIONAL:
            boot : if True. do bootstrap
            nb : number of bootstraps, default is 1000
        anti : Boolean
            if True, take antipodes of reverse poles
        spin : Boolean
            if True, transform data to spin axis
        rev : Boolean
            if True, take only reverse poles
        v : Boolean
            if True, filter data with Vandamme (1994) cutoff
        boot : Boolean
            if True, use bootstrap for confidence 95% interval
        mm97 : Boolean
            if True, use McFadden McElhinny 1997 correction for S
        nb : int
            number of bootstrapped pseudosamples for confidence estimate
        verbose : Boolean
            if True, print messages
        """
        if select:
            allsiteVGPdf=self.allsiteinfor[self.allsiteinfor['goodbool']].copy()
        else:
            allsiteVGPdf=self.allsiteinfor.copy()
        vgpdf_mean=pd.DataFrame()
        NVGPlist=allsiteVGPdf.loc[allsiteVGPdf['polarity']=='N',['Plon','Plat']].values.tolist()
        if NVGPlist: 
            Nvgpdf=pd.DataFrame(fisher_mean(di_block=NVGPlist))
            Nvgpdf['polarity']='N'
            vgpdf_mean=pd.concat([vgpdf_mean, Nvgpdf], axis=0)
            #vgplist.append(fisher_mean(di_block=NVGPlist))
        else:
            print("Warning: No normal polarity (N) data found.")
            Nvgpdf=pd.DataFrame({'dec': np.nan, 'inc': np.nan, 'n': 0, 'k': np.nan, 'alpha95': np.nan}, index=[0])

        RVGPlist = allsiteVGPdf.loc[allsiteVGPdf['polarity']=='R',['Plon','Plat']].values.tolist()
        if RVGPlist: 
            Rvgpdf=pd.DataFrame(fisher_mean(di_block=RVGPlist))
            Rvgpdf['polarity']='R'
            vgpdf_mean=pd.concat([vgpdf_mean, Rvgpdf], axis=0)
        else:
            print("Warning: No reverse polarity (R) data found.")
            Rvgpdf=pd.DataFrame({'dec': np.nan, 'inc': np.nan, 'n': 0, 'k': np.nan, 'alpha95': np.nan}, index=[0])
        FVGPlist = allsiteVGPdf.loc[:,['fPlon','fPlat']].values.tolist()
        if FVGPlist:
            Fvgpdf=pd.DataFrame(fisher_mean(di_block=FVGPlist))
            Fvgpdf['polarity']='F'
            vgpdf_mean=pd.concat([vgpdf_mean, Fvgpdf], axis=0)
        else:
            print("Warning: No filtered polarity (F) data found.")
            Fvgpdf=pd.DataFrame({'dec': np.nan, 'inc': np.nan, 'n': 0, 'k': np.nan, 'alpha95': np.nan}, index=[0])
        vgpdf_mean.rename(columns={'dec': 'Plon', 'inc': 'Plat', 'alpha95': 'A95'}, inplace=True)
        SBdf=pd.DataFrame([])
        if len(allsiteVGPdf.loc[allsiteVGPdf['polarity']=='N',:]) > 0:
        # only N polarity
            #No cut
            self.N_SBdf = scalc_vgp_df(allsiteVGPdf.loc[allsiteVGPdf['polarity']=='N',:], anti=False, rev=False, cutoff=180., 
                        kappa=kappa, n=n, spin=spin, v=False, boot=boot, 
                        mm97=mm97, nb=nb)
            self.N_SBdf['polarity'] = 'N'
            #45°Cut
            self.N_SBdf_45 = scalc_vgp_df(allsiteVGPdf.loc[allsiteVGPdf['polarity']=='N',:], anti=False, rev=False, cutoff=45., 
                        kappa=kappa, n=n, spin=spin, v=False, boot=boot, 
                        mm97=mm97, nb=nb)
            self.N_SBdf_45['polarity'] = 'N'
            #Vcut
            self.N_SBdf_V = scalc_vgp_df(allsiteVGPdf.loc[allsiteVGPdf['polarity']=='N',:], anti=False, rev=False, cutoff=180., 
                        kappa=kappa, n=n, spin=spin, v=True, boot=boot, 
                        mm97=mm97, nb=nb)
            self.N_SBdf_V['polarity'] = 'N'
            self.N_SBdf = pd.concat([self.N_SBdf, self.N_SBdf_45, self.N_SBdf_V], axis=0).reset_index(drop=True)
        else:
            print("Warning: No normal polarity (N) data found for S calculations.")
            self.N_SBdf = pd.DataFrame({'nVGP': [0], 'VGP_S': [np.nan], 'VGP_Sb': [np.nan], 
                                'Sb_low': [np.nan], 'Sb_high': [np.nan], 'Cutoff': [np.nan], 
                                'polarity': ['N']})
        if len(allsiteVGPdf.loc[allsiteVGPdf['polarity']=='R',:]) > 0:
            #No cut
            self.R_SBdf = scalc_vgp_df(allsiteVGPdf.loc[allsiteVGPdf['polarity']=='R',:], anti=True, rev=False, cutoff=180, 
                        kappa=kappa, n=n, spin=spin, v=False, boot=boot, 
                        mm97=mm97, nb=nb)
            self.R_SBdf['polarity'] = 'R'
            #45°Cut
            self.R_SBdf_45 = scalc_vgp_df(allsiteVGPdf.loc[allsiteVGPdf['polarity']=='R',:], anti=True, rev=False, cutoff=45., 
                        kappa=kappa, n=n, spin=spin, v=False, boot=boot, 
                        mm97=mm97, nb=nb)
            self.R_SBdf_45['polarity'] = 'R'
            #Vcut
            self.R_SBdf_V = scalc_vgp_df(allsiteVGPdf.loc[allsiteVGPdf['polarity']=='R',:], anti=True, rev=False, cutoff=180., 
                        kappa=kappa, n=n, spin=spin, v=True, boot=boot, 
                        mm97=mm97, nb=nb)
            self.R_SBdf_V['polarity'] = 'R'
            self.R_SBdf = pd.concat([self.R_SBdf, self.R_SBdf_45, self.R_SBdf_V], axis=0).reset_index(drop=True)
        else:
            print("Warning: No reversed polarity (R) data found for S calculations.")
            self.R_SBdf = pd.DataFrame({'nVGP': [0], 'VGP_S': [np.nan], 'VGP_Sb': [np.nan], 
                                'Sb_low': [np.nan], 'Sb_high': [np.nan], 'Cutoff': [np.nan], 
                                'polarity': ['R']})

        if len(allsiteVGPdf.loc[:,['fPlon','fPlat']]) > 0:
            #No cut
            flipVGPdf=allsiteVGPdf.copy()
            flipVGPdf.loc[:,['Plon','Plat']]= allsiteVGPdf.loc[:,['fPlon','fPlat']].values

            self.F_SBdf = scalc_vgp_df(flipVGPdf, anti=False, rev=False, cutoff=180.0, 
                        kappa=kappa, n=n, spin=spin, v=False, boot=boot, 
                        mm97=mm97, nb=nb)
            self.F_SBdf['polarity'] = 'F'
            #45°Cut
            self.F_SBdf_45 = scalc_vgp_df(flipVGPdf, anti=False, rev=False, cutoff=45., 
                        kappa=kappa, n=n, spin=spin, v=False, boot=boot, 
                        mm97=mm97, nb=nb)
            self.F_SBdf_45['polarity'] = 'F'
            #Vcut
            self.F_SBdf_V = scalc_vgp_df(flipVGPdf, anti=False, rev=False, cutoff=180.,
                        kappa=kappa, n=n, spin=spin, v=True, boot=boot, 
                        mm97=mm97, nb=nb)
            self.F_SBdf_V['polarity'] = 'F'
            self.F_SBdf = pd.concat([self.F_SBdf, self.F_SBdf_45, self.F_SBdf_V], axis=0).reset_index(drop=True)
        else:

            print("Warning: No filtered data found for S calculations.")
            self.F_SBdf = pd.DataFrame({'nVGP': [0], 'VGP_S': [np.nan], 'VGP_Sb': [np.nan], 
                                'Sb_low': [np.nan], 'Sb_high': [np.nan], 'Cutoff': [np.nan], 
                                'polarity': ['F']})

        self.SBdf = pd.concat([self.N_SBdf, self.R_SBdf, self.F_SBdf], axis=0).reset_index(drop=True)
        self.VGPdf_mean = vgpdf_mean.reset_index(drop=True)


    def siteslist(self):
        return list(self.site.keys())
    
    def sampleslist(self):
        all_samples = []
        for site in self.site.values():
            all_samples.extend(site.sampleslist())
        return all_samples
    
    def specimenslist(self):
        all_specimens = []
        for site in self.site.values():
            all_specimens.extend(site.specimenslist())
        return all_specimens
    
    def __getitem__(self, keys):
        if self._is_selection and self._selected_sites:
            specimen_keys = keys
            result = mqdirdata()   

            for site_name in self._selected_sites:
                if site_name in self.site:
                    result.site[site_name] = self.site[site_name]
            
            if isinstance(specimen_keys, (tuple, list)) and len(specimen_keys) > 0:
                for site_name in list(result.site.keys()):
                    site_obj = result.site[site_name]
                    for sample_name in list(site_obj.samples.keys()):
                        sample_obj = site_obj.samples[sample_name]
                        all_specimens = sample_obj.specimenslist()
                        to_remove = [spec for spec in all_specimens if spec not in specimen_keys]
                        
                        for spec_name in to_remove:
                            del sample_obj.specimens[spec_name]
                            
                        if not sample_obj.specimens:
                            del site_obj.samples[sample_name]
                            
                    if not site_obj.samples:
                        del result.site[site_name]
            
            elif isinstance(specimen_keys, str):
                for site_name in list(result.site.keys()):
                    site_obj = result.site[site_name]
                    for sample_name in list(site_obj.samples.keys()):
                        sample_obj = site_obj.samples[sample_name]
                        all_specimens = sample_obj.specimenslist()
                        to_remove = [spec for spec in all_specimens if spec != specimen_keys]
                        
                        for spec_name in to_remove:
                            del sample_obj.specimens[spec_name]
                            
                        if not sample_obj.specimens:
                            del site_obj.samples[sample_name]
                            
                    if not site_obj.samples:
                        del result.site[site_name]

            if hasattr(self, 'params_df') and result.site:
                result.create_params_df()
                if hasattr(self, 'PassName') and hasattr(self, 'PassBresult'):
                    result.pass_reset()

            self._is_selection = False
            self._selected_sites = None
            
            return result

        if isinstance(keys, (tuple, list)):
            if len(keys) == 3:  # site, sample, specimen
                site_name, sample_name, specimen_name = keys
                return self.site[site_name].samples[sample_name].specimens[specimen_name]
            elif len(keys) == 2:  # site, sample
                site_name, sample_name = keys
                return self.site[site_name].samples[sample_name]
            elif len(keys) == 1:  # site
                site_name = keys[0]
                return self.site[site_name]
            else:
                self._is_selection = True
                self._selected_sites = [k for k in keys if k in self.site]
                return self
        elif isinstance(keys, str):
            if keys in self.site:
                return self.site[keys]
            raise KeyError(f"Site '{keys}' not found")
        else:
            raise TypeError("Invalid key type")
            
    def __setitem__(self, keys, value):
        if isinstance(keys, (tuple, list)):
            if len(keys) == 3:  # site, sample, specimen
                site_name, sample_name, specimen_name = keys
                if site_name not in self.site:
                    self.site[site_name] = Site(site_name)
                if sample_name not in self.site[site_name].samples:
                    self.site[site_name].samples[sample_name] = Sample(sample_name)
                self.site[site_name].samples[sample_name].specimens[specimen_name] = value
            elif len(keys) == 2:  # site, sample
                site_name, sample_name = keys
                if site_name not in self.site:
                    self.site[site_name] = Site(site_name)
                self.site[site_name].samples[sample_name] = value
            elif len(keys) == 1:  # site
                site_name = keys[0]
                self.site[site_name] = value
            else:
                for key in keys:
                    if key in self.site:
                        self.site[key] = value
        elif isinstance(keys, str):
            self.site[keys] = value
        else:
            raise TypeError("Invalid key type")
            
    def __delitem__(self, keys):
        if isinstance(keys, (tuple, list)):
            if len(keys) == 3:  # site, sample, specimen
                site_name, sample_name, specimen_name = keys
                if site_name in self.site:
                    if sample_name in self.site[site_name].samples:
                        self.site[site_name].samples[sample_name].__delitem__(specimen_name)
                        if not self.site[site_name].samples[sample_name].specimens:
                            del self.site[site_name].samples[sample_name]
                        if not self.site[site_name].samples:
                            del self.site[site_name]
                    else:
                        raise KeyError(f"Sample '{sample_name}' not found in site '{site_name}'")
                else:
                    raise KeyError(f"Site '{site_name}' not found")
            elif len(keys) == 2:  # site, sample
                site_name, sample_name = keys
                if site_name in self.site:
                    self.site[site_name].__delitem__(sample_name)
                    if not self.site[site_name].samples:
                        del self.site[site_name]
                else:
                    raise KeyError(f"Site '{site_name}' not found")
            elif len(keys) == 1:  # site
                site_name = keys[0]
                if site_name in self.site:
                    del self.site[site_name]
                else:
                    raise KeyError(f"Site '{site_name}' not found")
            else:
                for site_name in keys:
                    if site_name in self.site:
                        del self.site[site_name]
        elif isinstance(keys, str):
            if keys in self.site:
                del self.site[keys]
            else:
                raise KeyError(f"Site '{keys}' not found")
        else:
            raise TypeError("Invalid key type")
            
        if hasattr(self, 'params_df'):
            self.create_params_df()
        if hasattr(self, 'PassName') and hasattr(self, 'PassBresult'):
            self.pass_reset()


def Makedirobj(generalpath, inforpath, redopath=None, Sitelen=3, Samplelen=4,momentunit='$Am^2$'):
    mqdirdataobj=mqdirdata(generalpath, inforpath, redopath=redopath, Sitelen=Sitelen, Samplelen=Samplelen,momentunit=momentunit)
    return mqdirdataobj


def redo_pca(redodf,generaldf,treatment_type,specimen_name):

    redodf[0]=redodf[0].str.replace('current_', '')
    #Single_redo=redodf
    vals2 = redodf.iloc[:, 2].to_numpy()
    vals3 = redodf.iloc[:, 3].to_numpy()
    redodf.iloc[:, 2] = np.where(vals2 < 1, vals2 * 1000, vals2 - 273)
    redodf.iloc[:, 3] = np.where(vals3 < 1, vals3 * 1000, vals3 - 273)
    #Single_redo = redodf.copy()
    fit_number=len(redodf)
    #PCAmethod=pd.DataFrame(['DE-BFL', 'DE-BFL-A', 'DE-BFL-O','DE-BFP', 'DE-FM'])
    for u in range(fit_number):
        fitname=redodf.iloc[u,5]
        Single_redo=redodf.loc[u,:]
        starti=generaldf.loc[generaldf['treatment']==Single_redo[2],:].index[0]
        endi=generaldf.loc[generaldf['treatment']==Single_redo[3],:].index[0]
        selectlow=Single_redo[2]
        selecthigh=Single_redo[3]
        selectbool=(generaldf['treatment']>=selectlow)&(generaldf['treatment']<=selecthigh)
        PCAmethod=Single_redo[1]
        try:
            testq=generaldf['quality']
        except:
            generaldf['quality']='g'
        try:
            data_g=generaldf[['treatment','dec_g','inc_g','moment','quality']].values.tolist()
            pca_g=pd.DataFrame(pmag.domean(data_g,starti,endi,PCAmethod))
            fisher_g=pd.DataFrame(pmag.domean(data_g,starti,endi,'DE-FM'))
            pca_g['specimen']=specimen_name[0]
            pca_g['fit_name']=fitname
            fisher_g['specimen']=specimen_name[0]
            fisher_g['fit_name']=fitname
            fisher_g=fisher_g.drop_duplicates(keep='first')
            center_of_mass_col = pca_g.pop('center_of_mass')
            pca_g=pca_g.drop_duplicates(keep='first')
            transposed_center_of_mass = center_of_mass_col.values.reshape(1, -1)
            new_columns = ['center_of_mass1', 'center_of_mass2', 'center_of_mass3']
            df_transposed = pd.DataFrame(transposed_center_of_mass, columns=new_columns)
            pca_g = pd.concat([pca_g, df_transposed], axis=1)
        except:
            pca_g=pd.DataFrame([])
            fisher_g=pd.DataFrame([])
        try:
            data_t=generaldf[['treatment','dec_t','inc_t','moment','quality']].values.tolist()
            pca_t=pd.DataFrame(pmag.domean(data_t,starti,endi,PCAmethod))
            fisher_t=pd.DataFrame(pmag.domean(data_t,starti,endi,'DE-FM'))
            pca_t['specimen']=specimen_name[0]
            pca_t['fit_name']=fitname
            fisher_t['specimen']=specimen_name[0]
            fisher_t['fit_name']=fitname
            fisher_t=fisher_t.drop_duplicates(keep='first')
            center_of_mass_col = pca_t.pop('center_of_mass')
            pca_t=pca_t.drop_duplicates(keep='first')
            transposed_center_of_mass = center_of_mass_col.values.reshape(1, -1)
            new_columns = ['center_of_mass1', 'center_of_mass2', 'center_of_mass3']
            df_transposed = pd.DataFrame(transposed_center_of_mass, columns=new_columns)
            pca_t = pd.concat([pca_t, df_transposed], axis=1)
        except:
            pca_t=pd.DataFrame([])
            fisher_t=pd.DataFrame([])
        try:
            data_s=generaldf[['treatment','dec_s','inc_s','moment','quality']].values.tolist()
            pca_s=pd.DataFrame(pmag.domean(data_s,starti,endi,PCAmethod))
            fisher_s=pd.DataFrame(pmag.domean(data_s,starti,endi,'DE-FM'))
            pca_s['specimen']=specimen_name[0]
            pca_s['fit_name']=fitname
            fisher_s['specimen']=specimen_name[0]
            fisher_s['fit_name']=fitname
            fisher_s=fisher_s.drop_duplicates(keep='first')
            center_of_mass_col = pca_s.pop('center_of_mass')
            pca_s=pca_s.drop_duplicates(keep='first')
            transposed_center_of_mass = center_of_mass_col.values.reshape(1, -1)
            new_columns = ['center_of_mass1', 'center_of_mass2', 'center_of_mass3']
            df_transposed = pd.DataFrame(transposed_center_of_mass, columns=new_columns)
            pca_s = pd.concat([pca_s, df_transposed], axis=1)
        except:
            pca_s=pd.DataFrame([])
            fisher_s=pd.DataFrame([])
    return selectbool,selectlow,selecthigh,pca_g,fisher_g,pca_t,fisher_t,pca_s,fisher_s



def Tmatrix(X):
    T = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]
    for row in X:
        for k in range(3):
            for l in range(3):
                T[k][l] += row[k] * row[l]
    return T
def tauV(T):

    t, V, tr = [], [], 0.
    ind1, ind2, ind3 = 0, 1, 2
    evalues, evectmps = linalg.eig(T)
    # to make compatible with Numeric convention
    evectors = np.transpose(evectmps)
    for tau in evalues:
        tr += tau
    if tr != 0:
        for i in range(3):
            evalues[i] = evalues[i] / tr
    else:
        return t, V
# sort evalues,evectors
    t1, t2, t3 = 0., 0., 1.
    for k in range(3):
        if evalues[k] > t1:
            t1, ind1 = evalues[k], k
        if evalues[k] < t3:
            t3, ind3 = evalues[k], k
    for k in range(3):
        if evalues[k] != t1 and evalues[k] != t3:
            t2, ind2 = evalues[k], k
    V.append(evectors[ind1])
    V.append(evectors[ind2])
    V.append(evectors[ind3])
    t.append(t1)
    t.append(t2)
    t.append(t3)
    return t, V
def vector_mean(data):
    Xbar = np.zeros((3))
    X = dir2cart(data).transpose()
    for i in range(3):
        Xbar[i] = X[i].sum()
    R = np.sqrt(Xbar[0]**2+Xbar[1]**2+Xbar[2]**2)
    Xbar = Xbar/R
    dir = cart2dir(Xbar)
    return dir, R
def domean(data, start, end, calculation_type):

    mpars = {}
    datablock = []
    start0, end0 = start, end
    indata = []
    for rec in data:
        if len(rec) < 6:
            rec.append('g')
        indata.append(rec)
    if indata[start0][5] == 'b':
        print("Can't select 'bad' point as start for PCA")
    flags = [x[5] for x in indata]
    bad_before_start = flags[:start0].count('b')
    bad_in_mean = flags[start0:end0 + 1].count('b')
    start = start0 - bad_before_start
    end = end0 - bad_before_start - bad_in_mean
    datablock = [x for x in indata if x[5] == 'g']
    if indata[start0] != datablock[start]:
        print('problem removing bad data in pmag.domean start of datablock shifted:\noriginal: %d\nafter removal: %d' % (
            start0, indata.index(datablock[start])))
    if indata[end0] != datablock[end]:
        print('problem removing bad data in pmag.domean end of datablock shifted:\noriginal: %d\nafter removal: %d' % (
            end0, indata.index(datablock[end])))
    mpars["calculation_type"] = calculation_type
    rad = np.pi/180.
    if end > len(datablock) - 1 or end < start:
        end = len(datablock) - 1
    control, data, X, Nrec = [], [], [], float(end - start + 1)
    cm = [0., 0., 0.]
    fdata = []
    for k in range(start, end + 1):
        if calculation_type == 'DE-BFL' or calculation_type == 'DE-BFL-A' or calculation_type == 'DE-BFL-O':  # best-fit line
            data = [datablock[k][1], datablock[k][2], datablock[k][3]]
        else:
            data = [datablock[k][1], datablock[k][2], 1.0]  # unit weight
        fdata.append(data)
        cart = dir2cart(data)
        X.append(cart)
    if calculation_type == 'DE-BFL-O':  # include origin as point
        X.append([0., 0., 0.])
        # pass
    if calculation_type == 'DE-FM':  # for fisher means
        fpars = fisher_mean(fdata)
        mpars["specimen_direction_type"] = 'l'
        mpars["specimen_dec"] = fpars["dec"]
        mpars["specimen_inc"] = fpars["inc"]
        mpars["specimen_alpha95"] = fpars["alpha95"]
        mpars["specimen_n"] = fpars["n"]
        mpars["specimen_r"] = fpars["r"]
        mpars["measurement_step_min"] = indata[start0][0]
        mpars["measurement_step_max"] = indata[end0][0]
        mpars["center_of_mass"] = cm
        mpars["specimen_dang"] = -1
        return mpars
    for cart in X:
        for l in range(3):
            cm[l] += cart[l] / Nrec
    mpars["center_of_mass"] = cm
    if calculation_type != 'DE-BFP':
        mpars["specimen_direction_type"] = 'l'
    if calculation_type == 'DE-BFL' or calculation_type == 'DE-BFL-O':  # not for planes or anchored lines
        for k in range(len(X)):
            for l in range(3):
                X[k][l] = X[k][l] - cm[l]
    else:
        mpars["specimen_direction_type"] = 'p'
    T = np.array(Tmatrix(X))
    t, V = tauV(T)
    if t == []:
        mpars["specimen_direction_type"] = "Error"
        print("Error in calculation")
        return mpars
    v1, v3 = V[0], V[2]
    if t[2] < 0:
        t[2] = 0  # make positive
    if calculation_type == 'DE-BFL-A':
        Dir, R = vector_mean(fdata)
        mpars["specimen_direction_type"] = 'l'
        mpars["specimen_dec"] = Dir[0]
        mpars["specimen_inc"] = Dir[1]
        mpars["specimen_n"] = len(fdata)
        mpars["measurement_step_min"] = indata[start0][0]
        mpars["measurement_step_max"] = indata[end0][0]
        mpars["center_of_mass"] = cm
        s1 = np.sqrt(t[0])
        MAD = np.arctan(np.sqrt(t[1] + t[2]) / s1) / rad
        if np.iscomplexobj(MAD):
            MAD = MAD.real
        # I think this is how it is done - i never anchor the "PCA" - check
        mpars["specimen_mad"] = MAD
        return mpars
    if calculation_type != 'DE-BFP':
        #
        #   get control vector for principal component direction
        #
        rec = [datablock[start][1], datablock[start][2], datablock[start][3]]
        P1 = dir2cart(rec)
        rec = [datablock[end][1], datablock[end][2], datablock[end][3]]
        P2 = dir2cart(rec)
#
#   get right direction along principal component
##
        for k in range(3):
            control.append(P1[k] - P2[k])
        # changed by rshaar
        # control is taken as the center of mass
        # control=cm

        dot = 0
        for k in range(3):
            dot += v1[k] * control[k]
        if dot < -1:
            dot = -1
        if dot > 1:
            dot = 1
        if np.arccos(dot) > old_div(np.pi, 2.):
            for k in range(3):
                v1[k] = -v1[k]
#   get right direction along principal component
#
        s1 = np.sqrt(t[0])
        Dir = cart2dir(v1)
        MAD = old_div(np.arctan(old_div(np.sqrt(t[1] + t[2]), s1)), rad)
        if np.iscomplexobj(MAD):
            MAD = MAD.real
    if calculation_type == "DE-BFP":
        Dir = cart2dir(v3)
        MAD = old_div(
            np.arctan(np.sqrt(old_div(t[2], t[1]) + old_div(t[2], t[0]))), rad)
        if np.iscomplexobj(MAD):
            MAD = MAD.real
#
#   get angle with  center of mass
#
    CMdir = cart2dir(cm)
    Dirp = [Dir[0], Dir[1], 1.]
    dang = angle(CMdir, Dirp)
    mpars["specimen_dec"] = Dir[0]
    mpars["specimen_inc"] = Dir[1]
    mpars["specimen_mad"] = MAD
    # mpars["specimen_n"]=int(Nrec)
    mpars["specimen_n"] = len(X)
    mpars["specimen_dang"] = dang[0]
    mpars["measurement_step_min"] = indata[start0][0]
    mpars["measurement_step_max"] = indata[end0][0]
    return mpars







def equalareaplot(dec=None,inc=None,a95=None, legendbool=True,
                  cirlabel=None,circolor='k',
                  fisher=False,fisherdf=None,fishertextloc=None,
                  ax=None,
                  type=0,line=True,showticks=False,markersize=4,linewidth=1,steptext=None,
                  starlabel=None,starcolor='k',):
    """Create an equal-area plot.
       inc>=0:down k
        else:up w
    """
    if dec is not None:
        if isinstance(dec, (int, float)):
            dec = [dec]
        elif not isinstance(dec, list):
            dec = list(dec)  # 处理numpy数组等其他类型
    
    # 检查并转换 inc
    if inc is not None:
        if isinstance(inc, (int, float)):
            inc = [inc]
        elif not isinstance(inc, list):
            inc = list(inc)
    
    # 检查并转换 a95
    if a95 is not None:
        if isinstance(a95, (int, float)):
            a95 = [a95]
        elif not isinstance(a95, list):
            a95 = list(a95)
    if ax is None:
        fig, ax = upt.subplots(figsize=(18/2.54, 18/2.54),
        top='1.3cm',
        bottom='0.8cm',
        right='0.5cm',
        left='0.5cm',
         ncols=1, nrows=1,share=False,proj="polar", )
    else:
        pass
    if dec is None or inc is None:
        pass
    else:
        dec_c=np.radians(dec)
        inc_c=EqualArea(inc)
        downbool = np.array(inc) >= 0
        upbool = np.array(inc) < 0
        downc= None
        upc=None
        try:
            downc=ax.plot(dec_c[downbool],inc_c[downbool],'o',facecolor=circolor,edgecolor=circolor, markersize=markersize, label='Down',zorder=100,clip_on=False)
        except:
            pass
        try:
            upc=ax.plot(dec_c[upbool],inc_c[upbool],'o',facecolor='w', edgecolor=circolor,markersize=markersize, label='Up',zorder=100,clip_on=False)
        except:
            pass
        if line:
            ax.plot(dec_c,inc_c,'k--',linewidth=linewidth) 
        if steptext is not None:
            pass
        # EQsteptext(ax, dec_c, inc_c, steptext, step=3)
    if type == 0:
        angles = [0, 90, 180, 270]
        angle_ranges = [(10, 100, 10)] + [(10, 90, 10)] * 3 
        symbols = []
        for i, (angle, (start, end, step)) in enumerate(zip(angles, angle_ranges)):
            Xsym, Ysym = [], []
            for Pl in range(start, end, step):
                ix = EqualArea(Pl)
                # if angle == 0:  
                #     print(ix)
                Xsym.append(np.radians(angle))
                Ysym.append(ix)
            symbol = ax.plot(Xsym, Ysym, marker='+', color='gray', markersize=markersize+2, linestyle='None')
            symbols.append(symbol)
        list1=None
        list2=None
        grid=False
        ax.set_rticks([])
    elif type == 1:
        list1=[EqualArea(x) for x in range(10,90,20)]
        if showticks:
            list2=[str(x) for x in range(10,90,20)]
        else:
            list2="null"#[str(x) for x in range(10,90,20)]
        grid=True

    if a95 is not None:
        for i, a95_value in enumerate(a95):
            if a95_value is not None:
                dec_i = dec[i]
                inc_i = inc[i]
                
                xpos, ypos, xneg, yneg = alphacirc(dec_i, inc_i, alpha=a95_value)
                
                try:
                    ax.plot(xpos, ypos, 'k-', linewidth=0.5,alpha=0.7)
                except:
                    pass
                try:

                    ax.plot(xneg, yneg, 'k-', linewidth=0.5,alpha=0.7)
                except:
                    pass
    if fisher:
        # Plot Fisher mean direction
        if fisherdf is not None:
            fisherresualt=fisherdf
            try:
                dec_f, inc_f,alpha95 = fisherresualt.specimen_dec[0], fisherresualt.specimen_inc[0], fisherresualt.specimen_alpha95[0]
            except:
                dec_f, inc_f,alpha95 = fisherresualt.dec[0], fisherresualt.inc[0], fisherresualt.alpha95[0]
            
        else:
            try:
                fisherresualt =fisher_mean(dec, inc)
                dec_f, inc_f,alpha95 = fisherresualt.dec[0], fisherresualt.inc[0], fisherresualt.alpha95[0]
            except:
                pass
        downS = None
        upperS = None
        if inc_f < 0:
            #Up direction
            # upperS=ax.plot(np.radians(dec_f), EqualArea(inc_f), '*',facecolor='w', markersize=markersize+3,markeredgewidth=0.2, edgecolor=starcolor,clip_on=False,zorder=1001,label='Fisher Mean Up',)
            upperS = ax.scatter(np.radians(dec_f), EqualArea(inc_f), marker='*', 
                facecolors='w', edgecolors=starcolor, s=(markersize+3)**2,
                linewidths=0.2, clip_on=False, zorder=10001, 
                label='Fisher Mean Up')
        else:
            #Down direction
            #downS=ax.plot(np.radians(dec_f), EqualArea(inc_f), '*',facecolor=starcolor, markersize=markersize+3,markeredgewidth=0.2,edgecolor=starcolor, clip_on=False,zorder=1001,label='Fisher Mean Down',)
            downS = ax.scatter(np.radians(dec_f), EqualArea(inc_f), marker='*', 
                facecolors=starcolor, edgecolors=starcolor, 
                s=(markersize+3)**2, linewidths=0.2, clip_on=False, 
                zorder=10001, label='Fisher Mean Down')
        
        fxpos, fypos, fxneg, fyneg = alphacirc(dec_f, inc_f, alpha=alpha95)
        ax.fill(fxpos, fypos, color='gray', alpha=0.5,zorder=10000,)
        ax.fill(fxneg, fyneg, color='gray', alpha=0.5,zorder=10000)
        if starlabel is not None and (downS is not None or upperS is not None):
            handles = []
            if downS is not None:

                handles.append(downS)
            if upperS is not None:

                handles.append(upperS)

            if handles:
                Slegend=ax.legend(handles=handles, label=starlabel,loc="t", ncols=2, frame=False)
        else:
            pass

        if fishertextloc is None:
            textf=False
        elif fishertextloc == 't':
            textf=True
            textx=0.5
            texty=0.8
        elif fishertextloc == 'b':
            textf=True
            textx=0.5
            texty=0.2
        elif fishertextloc == 'l':
            textf=True
            textx=0.2
            texty=0.5
        elif fishertextloc == 'r':
            textf=True
            textx=0.8
            texty=0.5
        else:
            textf=False
        if textf:
            nt='n='+str(int(fisherresualt.n.values))
            Dt='D='+str(round(float(fisherresualt.dec.values),2))
            It='I='+str(round(float(fisherresualt.inc.values),2))
            a95t='α_{95}='+str(round(float(fisherresualt.alpha95.values),2))
            kt='κ='+str(round(float(fisherresualt.k.values),2))
            alltext= r'${}${}${}${}${}${}${}${}${}$'.format(nt,'\n',Dt,'\n',It,'\n',a95t,'\n',kt)
            alltextP=ax.text(textx,texty,alltext,ha='center', va='center',color='k',fontsize=12,transform=ax.transAxes,clip_on=False)
            alltextP.set_bbox(dict(facecolor='w',linewidth=0, alpha=0.5, pad=0)) 
        else:
            pass
    else:
        pass
    downl=ax.scatter(0, 100, marker='o',markersize=markersize+10,color='k', edgecolor='k',label='Down')
    upperl=ax.scatter(0, 100, marker='o',markersize=markersize+10,color='w', edgecolor='k',label='Up') 

    if cirlabel is not None:
        llegend=ax.legend(handles=[downc, upc], label=cirlabel,loc="t", ncols=2, frame=False)
    else:
        if legendbool:
            plegend=ax.legend(handles=[downl, upperl],loc="t", ncols=2, frame=False)
        else:
            pass
    ax.format(
        thetadir=-1,
        thetalines=90,
        theta0="N",
        rlim=(90.,0.),
        grid=grid,
        rlines=list1,
        rformatter=list2,
    )

    return ax
def EQsteptext(ax, dec, inc, steptext, step=3):
    steps = []
    treatmtext =  steptext

    
    # Create array of valid indices (not NaN in x or y)
    x_array = np.array(dec)
    y_array = np.array(inc)
    valid_mask = ~np.isnan(x_array) & ~np.isnan(y_array)
    valid_indices = np.where(valid_mask)[0]
    xtext = x_array[::2]
    ytext = y_array[::2]
    steptextplot=np.array(steptext.values).flatten().tolist()
    steptextplot=steptextplot[::2]
    ax.text(xtext, ytext, steptextplot,
            ha='left', va='bottom', color='gray', fontsize=6, clip_on=False)
    # ta.allocate(ax,xtext,ytext,
    #         steptextplot,
    #         x_scatter=x_array, y_scatter=y_array,
    #         textsize=10)
    # if len(valid_indices) > step:
    #     # Select indices that are approximately evenly spaced
    #     step_size = max(1, len(valid_indices) // step)
    #     INDIC = valid_indices[::step_size][:step]
    # else:
    #     INDIC = valid_indices
    
    
    
    
    # for t in range(len(INDIC)):
    #     idx = INDIC[t]
    #     if idx < len(treatmtext):
    #         steps.append(ax.text(x_array[idx], y_array[idx], treatmtext.loc[idx, 'treatmtext'], 
    #                             ha='center', va='center', color='gray', fontsize=6, clip_on=False))

    # if steps:
    #     adjust_text(steps,ax=ax,)

def EqualArea(Pl):
    Pl_array = np.atleast_1d(Pl)
    Pl_adjusted = np.abs(Pl_array)
    result =90- np.sqrt(2.) * 90. * np.sin(np.radians(90. - Pl_adjusted) / 2.)

    
    if np.isscalar(Pl):
        return float(result[0])
    else:
        return result




def alphacirc(c_dec,c_inc,alpha):
    
    alpha=np.radians(alpha) 
    t=np.zeros((3,3)) 
    t[2]=dir2cart([c_dec,c_inc])
    plane1=[c_dec,c_inc-90.]
    plane2=[c_dec+90.,0]
    
    t[0]=dir2cart(plane1)
    t[1]=dir2cart(plane2)
    t=t.transpose()
    npts=201
    xnum=float(npts-1.)/2.
    v=[0,0,0]
    PTS=[]
    for i in range(npts): 
            psi=float(i)*np.pi/xnum
            v[0]=np.sin(alpha)*np.cos(psi)
            v[1]=np.sin(alpha)*np.sin(psi)
            if alpha==np.pi/2.:
                v[2]=0.
            else:
                v[2]=np.sqrt(1.-v[0]**2 - v[1]**2)
            elli=[0,0,0]
            for j in range(3):
                for k in range(3):
                    elli[j]=elli[j] + t[j][k]*v[k] 
            PTS.append(cart2dir(elli))
    pts_array = np.array(PTS)
    upper_mask = pts_array[:, 1] > 0
    lower_mask = ~upper_mask

    # Upper hemisphere (negative inclination)
    xpos = np.radians(pts_array[upper_mask, 0]).tolist()
    ypos = EqualArea(-pts_array[upper_mask, 1]).tolist()

    # Lower hemisphere (positive inclination)
    xneg = np.radians(pts_array[lower_mask, 0]).tolist()
    yneg = EqualArea(pts_array[lower_mask, 1]).tolist()
    if len(xpos) > 0 and len(xneg) > 0:
        rxpos=xpos[::-1]
        
        half = len(xneg) // 2
        xneg = xneg[half:] + xneg[:half]
        yneg = yneg[half:] + yneg[:half]
        rxneg=xneg[::-1]
        rypos=[0.0] * len(xpos)
        ryneg=[0.0] * len(xneg)
        xpos.extend(rxpos)
        ypos.extend(rypos)
        xneg.extend(rxneg)
        yneg.extend(ryneg)
    else:
        pass

    return xpos, ypos, xneg, yneg


def AllSiteeqplot(MQdirobj, ax=None,Coordinates='g',
                  marksize=2.5,savepath=None, meantype='PCA',
                  abc=True, abcloc='l',eqtype=0):
    #matplotlib.rcParams.update({'font.size':10})
    siteslist=MQdirobj.siteslist()
    fig,axes=mplot(pn=len(siteslist),cn=2,proj="polar",topc='1.5cm',bottomc='1cm',rightc='0.5cm',
        leftc='0.5cm',
        wspacec='1.5cm',
        hspacec='1.5cm',  )
    for i,site in enumerate(siteslist): 
        ax=axes[i]  
        singlesite=MQdirobj[site]
        if meantype == 'PCA':
            if Coordinates == 's':
                fisherdf=singlesite.sitePCA_s
                dec=singlesite.site_spe_s_PCA.specimen_dec.values
                inc=singlesite.site_spe_s_PCA.specimen_inc.values
                
                a95=singlesite.site_spe_s_PCA.specimen_mad.values
            elif Coordinates == 'g':
                fisherdf=singlesite.sitePCA_g
                dec=singlesite.site_spe_g_PCA.specimen_dec.values
                inc=singlesite.site_spe_g_PCA.specimen_inc.values
                a95=singlesite.site_spe_g_PCA.specimen_mad.values
            elif Coordinates == 't':
                fisherdf=singlesite.sitePCA_t
                dec=singlesite.site_spe_t_PCA.specimen_dec.values
                inc=singlesite.site_spe_t_PCA.specimen_inc.values
                a95=singlesite.site_spe_t_PCA.specimen_mad.values
        else:
            if Coordinates == 's':
                fisherdf=singlesite.siteFisher_s
                dec=singlesite.site_spe_s_Fisher.specimen_dec.values
                inc=singlesite.site_spe_s_Fisher.specimen_inc.values
                a95=singlesite.site_spe_s_Fisher.specimen_alpha95.values
            elif Coordinates == 'g':
                fisherdf=singlesite.siteFisher_g
                dec=singlesite.site_spe_g_Fisher.specimen_dec.values
                inc=singlesite.site_spe_g_Fisher.specimen_inc.values
                a95=singlesite.site_spe_g_Fisher.specimen_alpha95.values
            elif Coordinates == 't':
                fisherdf=singlesite.siteFisher_t
                dec=singlesite.site_spe_t_Fisher.specimen_dec.values
                inc=singlesite.site_spe_t_Fisher.specimen_inc.values
                a95=singlesite.site_spe_t_Fisher.specimen_alpha95.values

        ax=equalareaplot(dec,inc,a95=a95,
                              fisher=True,fisherdf=fisherdf,ax=ax,
                              type=eqtype,line=False,showticks=False,
                              markersize=marksize)
        nt='n='+str(int(fisherdf.n.values))
        Dt='D='+str(round(float(fisherdf.dec.values),2))
        It='I='+str(round(float(fisherdf.inc.values),2))
        a95t='α_{95}='+str(round(float(fisherdf.alpha95.values),2))
        kt='κ='+str(round(float(fisherdf.k.values),2))
        alltext= r'${}${}${}${}${}${}${}${}${}$'.format(nt,'\n',Dt,'\n',It,'\n',a95t,'\n',kt)
        if 225<float(fisherdf.inc.values)<315:
            textx=90
        else:
            textx=270
        alltextP=ax.text(np.radians(textx),30,alltext,ha='center', va='center',color='k',fontsize=12)
        alltextP.set_bbox(dict(facecolor='w',linewidth=0, alpha=0.5, pad=0))
        ax.format(title=site,titleloc='l',)
    axes.format(    
    abc=abc,
    abcloc=abcloc,
    )   
    if isinstance(axes.number, int):
        subplotpars=1
    else:
        subplotpars  = len(axes)
    pliax=subplotpars-i-1
    if pliax == 0:
        pass
    else:
        for of in range(0, pliax):
            axes[i+1+of].remove()
            #axes[i+1+of].axis('off')
    if savepath is not None:
        savepdf=savepath+'/'+'AllSiteresult.pdf'
        savepng=savepath+'/'+'AllSiteresult.png'
        fig.savefig(savepdf)
        fig.savefig(savepng)
    else:
        pass
    return fig,axes










def NSzplot(singlespecimen_obj,
            Coordinates='g',
            ax=None,
            PCA=False,
            markersize=4,
            linewidth=1,
            ):
    if ax is None:
        fig, ax = upt.subplots(figsize=(18/2.54, 18/2.54),ncols=1, nrows=1,share=False)
    else:
        juge=1
    specimen=singlespecimen_obj.specimen
    if Coordinates=='s':
        dec=singlespecimen_obj.dec_s
        inc=singlespecimen_obj.inc_s
        if PCA:
            PCAdf=singlespecimen_obj.pca_s
    elif Coordinates=='g':
        dec=singlespecimen_obj.dec_g
        inc=singlespecimen_obj.inc_g
        if PCA:
            PCAdf=singlespecimen_obj.pca_g

    elif Coordinates=='t':
        dec=singlespecimen_obj.dec_t
        inc=singlespecimen_obj.inc_t
        if PCA:
            PCAdf=singlespecimen_obj.pca_t
    moment=singlespecimen_obj.moment
    # treatment_typelist=singlespecimen_obj.treatment_type
    # treatment=singlespecimen_obj.treatment
    # treatment_unit=singlespecimen_obj.treatment_unit
    # treatmtext = pd.DataFrame({'treatmtext': treatment.astype(str) + treatment_unit['unite'].astype(str)})
    treatmtext = singlespecimen_obj.treatmenttext
    
    Iunit= singlespecimen_obj.momentunit

    VH = pmag.dir2cart(pd.concat([90-dec, -inc, moment], axis=1))
    ax.plot(VH[:,0], VH[:,1],
        color='black', linestyle='-', marker='o',markersize=markersize,linewidth=linewidth,markerfacecolor='black', markeredgecolor='black',label='N')
    ax.plot(VH[:,0], VH[:,2],
        color='black', linestyle='-', marker='o',markersize=markersize,linewidth=linewidth,markerfacecolor='white', markeredgecolor='black',label='UP')
    ax.plot(0, 0, alpha=0)
    #redo_pcaplot(ax,singlespecimen_obj.generaldf,VH,PCAdf,singlespecimen_obj.selectlow,singlespecimen_obj.selecthigh,singlespecimen_obj.selectbool,markersize=4)
    if PCA:
        try:

            redo_pcaplot(ax,singlespecimen_obj.generaldf,VH,PCAdf,singlespecimen_obj.selectlow,singlespecimen_obj.selecthigh,singlespecimen_obj.selectbool,markersize=markersize,markeredgewidth=None)
        except Exception as e:
            print("Error occurred while performing PCA:", e)
    else:
        pass

    maxX = np.max(np.maximum(0, VH[:,0]))
    maxY = np.max(np.maximum(0, np.maximum(VH[:,1], VH[:,2])))
    maxall = max(maxX, maxY)
    minX = np.min(np.minimum(0, VH[:,0]))
    minY = np.min(np.minimum(0, np.minimum(VH[:,1], VH[:,2])))
    minall = min(minX, minY)
    tick_step = calculate_nice_tick_step(minall, maxall)
    realrange=[minall-tick_step,maxall+1.5*tick_step]


    xvisrange=[minX-tick_step, maxX+tick_step]
    yvisrange=[minY-tick_step, maxY+tick_step]

    setup_partial_axis(ax, x_full_range=realrange, y_full_range=realrange, 
                        x_visible_range=xvisrange, y_visible_range=yvisrange)
    ax.xaxis.set_major_locator(MultipleLocator(tick_step))
    ax.yaxis.set_major_locator(MultipleLocator(tick_step))
    max_exp=int(np.floor(np.log10(tick_step)))
    
    ax.format(title=singlespecimen_obj.specimen[0], 
            grid=False,xloc='zero', yloc='zero',tickminor=False,ticklabelsize=6)
    if abs(max_exp) >= 2:
        ax.xaxis.set_major_formatter(precision_formatter(ax, 'x', max_exp))
        ax.yaxis.set_major_formatter(precision_formatter(ax, 'y', max_exp))

        xlabeltext = r'$\mathrm{E\;(\times 10^{' + str(max_exp) + r'} ' + Iunit.strip('$') + r')}$'
        ylabeltext = r'N ● UP ○ ($\times 10^{' + str(max_exp) + r'} ' + Iunit.strip('$') + r'$)'
    else:
        xlabeltext = r'$\mathrm{E\;(' + Iunit.strip('$') + r')}$'
        ylabeltext = r'N ● UP ○ ($' + Iunit.strip('$') + r'$)'        

    ax.text(0,yvisrange[1]+tick_step*0.15, ylabeltext, fontname='DejaVu Sans',ha='center', va='bottom',fontsize=6, clip_on=False)
    ax.text(xvisrange[1]+tick_step*0.3, -tick_step*0.15,xlabeltext,fontname='DejaVu Sans',va='center',fontsize=6, rotation=270, clip_on=False)

    steps = []
    treatmtext=np.array(treatmtext.values).flatten().tolist()
    plottreatmtext=treatmtext[::3]
    textVH=VH[::3]
    # if len(treatmtext) > 3:
    #     INDIC = np.arange(1, len(treatmtext), int(len(treatmtext)/3))
    # else:
    #     INDIC = np.arange(1, len(treatmtext),1)
    # NRMpositon=ax.text(np.array(VH)[0,0],np.array(VH)[0,1],'NRM',ha='center', va='center',color='gray',fontsize=6)
    # steps.append(NRMpositon)
    for t in range(len(textVH)):
        x = float(np.array(textVH)[t, 0])
        y = float(np.array(textVH)[t, 1])
        steps.append(ax.text(x, y,plottreatmtext[t],ha='center', va='center',color='gray',fontsize=6))
    adjust_text(steps, ax=ax)
    return ax
def calculate_nice_tick_step(data_min, data_max, density=0.5):
    data_range = max(data_max - data_min, 1e-10) 
    magnitude = 10 ** np.floor(np.log10(data_range))
    normalized_range = data_range / magnitude
    if normalized_range < 1.5:
        step = 0.25 * magnitude / density
    elif normalized_range < 3:
        step = 0.5 * magnitude / density 
    elif normalized_range < 7:
        step = 1.0 * magnitude / density
    else:
        step = 2.5 * magnitude / density
    min_ticks = 4
    if (data_max - data_min) / step < min_ticks:
        step /= 2
        
    return step


def setup_partial_axis(ax, x_full_range=None, y_full_range=None, x_visible_range=None, y_visible_range=None):

    if x_full_range is not None:
        ax.set_xlim(x_full_range)
        if x_visible_range is not None:
            class XPartialLocator(ticker.MaxNLocator):
                def __call__(self):
                    ticks = super().__call__()
                    return [t for t in ticks if x_visible_range[0] <= t <= x_visible_range[1]]
            ax.xaxis.set_major_locator(XPartialLocator())
            ax.spines['bottom'].set_bounds(x_visible_range[0], x_visible_range[1])
    
    if y_full_range is not None:

        ax.set_ylim(y_full_range)
        if y_visible_range is not None:
            class YPartialLocator(ticker.MaxNLocator):
                def __call__(self):
                    ticks = super().__call__()
                    return [t for t in ticks if y_visible_range[0] <= t <= y_visible_range[1]]
            ax.yaxis.set_major_locator(YPartialLocator())
            ax.spines['left'].set_bounds(y_visible_range[0], y_visible_range[1])
    return ax
# def precision_formatter(ax, axis='x', exp=0):
#     if axis == 'x':
#         ticks = ax.get_xticks()
#     else:
#         ticks = ax.get_yticks()
#     scaled_ticks = [t / (10**exp) for t in ticks if t != 0]
#     if len(scaled_ticks) > 1:
#         intervals = [abs(scaled_ticks[i] - scaled_ticks[i-1]) 
#                     for i in range(1, len(scaled_ticks))]
#         min_interval = min(intervals)
#         decimal_places = max(2, min(4, int(-np.log10(min_interval)) + 1))
#     else:
#         decimal_places = 2
#     def formatter(x, pos):
#         if x == 0:
#             return '0'
#         scaled_value = x/(10**exp)
#         if abs(scaled_value - round(scaled_value)) < 1e-10:
#             return f'{int(round(scaled_value))}'
#         formatted = f'{scaled_value:.{decimal_places}f}'
#         parts = formatted.split('.')
#         if len(parts) > 1:
#             decimal_part = parts[1]
#             first_zero_pos = len(decimal_part)
#             for i, digit in enumerate(decimal_part):
#                 if digit == '0':
#                     first_zero_pos = i
#                     break
#             if first_zero_pos > 0:
#                 parts[1] = decimal_part[:first_zero_pos]
#                 return f'{parts[0]}.{parts[1]}'
#             else:
#                 return parts[0]
#         return formatted
#     return ticker.FuncFormatter(formatter)
def precision_formatter(ax, axis='x', exp=0):
    if axis == 'x':
        ticks = ax.get_xticks()
    else:
        ticks = ax.get_yticks()
    scaled_ticks = [t / (10**exp) for t in ticks if abs(t) > 1e-15]
    
    if len(scaled_ticks) > 1:
        intervals = [abs(scaled_ticks[i] - scaled_ticks[i-1]) 
                    for i in range(1, len(scaled_ticks))]
        min_interval = min(intervals)
        decimal_places = max(2, min(4, int(-np.log10(min_interval)) + 1))
    else:
        decimal_places = 2
    
    def formatter(x, pos):
        if abs(x) < 1e-15:
            return '0'
        
        scaled_value = x / (10**exp)
        if abs(scaled_value - round(scaled_value)) < 1e-10:
            return f'{int(round(scaled_value))}'
        formatted = f'{scaled_value:.{decimal_places}f}'

        parts = formatted.split('.')
        if len(parts) > 1:
            decimal_part = parts[1].rstrip('0')  
            if decimal_part:  
                return f'{parts[0]}.{decimal_part}'
            else:  
                return parts[0]
        return formatted
    
    return ticker.FuncFormatter(formatter)
def redo_pcaplot(ax,SDATAI,VH,PCAdf,startstep,endstep,Select_bool,markersize=4,markeredgewidth=0.3):
    # Perform PCA and plot the results
    for u in range(len(PCAdf)):
        
        Spcadf=PCAdf.loc[[u],:]
        PREVH=VH[Select_bool,:]
        starti=SDATAI.loc[SDATAI['treatment']==startstep,:].index[0]
        endi=SDATAI.loc[SDATAI['treatment']==endstep,:].index[0]
        startint=SDATAI.loc[starti,'moment']
        endint=SDATAI.loc[endi,'moment']
        if Spcadf['calculation_type'][u]=='DE-BFL-A':
            HPCADIR=pmag.dir2cart([90-Spcadf.loc[u,'specimen_dec'],-Spcadf.loc[u,'specimen_inc'],startint])#HARROWSCALE
            Varrowx=HPCADIR[0]
            Varrowy=HPCADIR[1]
            Harrowx=HPCADIR[0]
            Harrowy=HPCADIR[2]
            Vtextx=0
            Vtexty=0
            Htextx=0
            Htexty=0
        elif Spcadf['calculation_type'][u]=='DE-BFL-O':

            HPCADIR=pmag.dir2cart([90-Spcadf.loc[u,'specimen_dec'],-Spcadf.loc[u,'specimen_inc'],startint])#HARROWSCALE
            VPCADIR=pmag.dir2cart([90-Spcadf.loc[u,'specimen_dec'],-Spcadf.loc[u,'specimen_inc'],endint])
            halfalllength=(HPCADIR-VPCADIR)/2*1.2
            zero_row = np.zeros((1, PREVH.shape[1]))
            meanpoint=np.mean(np.vstack((PREVH, zero_row)), axis=0)
            arrowpoint=meanpoint+halfalllength
            textpoint=meanpoint-halfalllength
            Varrowx=arrowpoint[0]
            Varrowy=arrowpoint[1]
            Harrowx=arrowpoint[0]
            Harrowy=arrowpoint[2]
            Vtextx=textpoint[0]
            Vtexty=textpoint[1]
            Htextx=textpoint[0]
            Htexty=textpoint[2]
        else:
            HPCADIR=pmag.dir2cart([90-Spcadf.loc[u,'specimen_dec'],-Spcadf.loc[u,'specimen_inc'],startint])#HARROWSCALE
            VPCADIR=pmag.dir2cart([90-Spcadf.loc[u,'specimen_dec'],-Spcadf.loc[u,'specimen_inc'],endint])
            halfalllength=(HPCADIR-VPCADIR)/2*1.2
            meanpoint=np.mean(PREVH, axis=0)
            arrowpoint=meanpoint+halfalllength
            textpoint=meanpoint-halfalllength
            Varrowx=arrowpoint[0]
            Varrowy=arrowpoint[1]
            Harrowx=arrowpoint[0]
            Harrowy=arrowpoint[2]
            Vtextx=textpoint[0]
            Vtexty=textpoint[1]
            Htextx=textpoint[0]
            Htexty=textpoint[2]
        arrow_style = "simple,head_length=0.6,head_width=0.6,tail_width=0.2"
        ax.plot(PREVH[:,0],PREVH[:,1],'o',markersize=markersize,markerfacecolor='none',markeredgecolor='#E45F2B',markeredgewidth=0.5,)   
        ax.plot(PREVH[:,0],PREVH[:,2],'o',markersize=markersize,markerfacecolor='none',markeredgecolor='#2A88FA',markeredgewidth=0.5,)   

        ax.annotate("",xy=(Harrowx*1.3,Harrowy*1.3),\
                    xytext=(Htextx,Htexty),\
                    arrowprops=dict(arrowstyle=arrow_style,mutation_scale=10, linewidth=0.1,color='#2A88FA',alpha=0.7),
                    annotation_clip=False)
        ax.annotate("",xy=(Varrowx*1.3,Varrowy*1.3),\
                    xytext=(Vtextx,Vtexty),\
                    arrowprops=dict(arrowstyle=arrow_style,mutation_scale=10,linewidth=0.1,color='#E45F2B',alpha=0.7),
                    annotation_clip=False)

def Smagplot(SingSpecimenObj,select=None,markersize=4,linewidth=1,
             ax=None,norm=True,):
    """
    treatment unit AF=mT, TD=°C
    """
    specimenname= SingSpecimenObj.specimen
    treatment= SingSpecimenObj.treatment
    treatment_type= SingSpecimenObj.treatment_type
    moment= SingSpecimenObj.moment
    momentunit= SingSpecimenObj.momentunit
    NRM=SingSpecimenObj.moment[0]
    c1 = "red8"
    c2 = "blue8"
    if ax is None:
        fig, ax = upt.subplots(figsize=(18/2.54, 18/2.54),ncols=1, nrows=1,share=False)
    else:
        juge=1
    
    AFbool=treatment_type=='A'
    TDbool=treatment_type=='T'
    AFmoment=moment[AFbool]
    TDmoment=moment[TDbool]
    AFtreatment=treatment[AFbool]
    TDtreatment=treatment[TDbool]
    if norm:
        ymax=NRM
        ylabeltext="$M/M_{NRM}$"
        try:
            AFmoment=AFmoment/ymax
        except:
            pass
        try:
            TDmoment=TDmoment/ymax
        except:
            pass
    else:
        pass
    AFl=None
    TDl=None
    if select is None:
        if len(AFtreatment) == 0:
            #pur TD
            TDl=ax.plot(TDtreatment,TDmoment , 'o-', color=c1, markeredgewidth= 0.5 , linewidth=linewidth, markersize= markersize , alpha=0.7,label="TD")
            maxY = max(TDmoment)
            ax.format(xlabel="Step (°C)",xcolor=c1, xmin=0, ymin=0) 
        elif len(TDtreatment) == 0:
            #pur AF
            AFl=ax.plot(AFtreatment,AFmoment , 'o-', color=c2, markeredgewidth= 0.5 , linewidth=linewidth, markersize= markersize , alpha=0.7,label="AFD")
            maxY = max(AFmoment)
            ax.format(xlabel="Step (mT)",xcolor=c1, xmin=0, ymin=0)
        else:

            #both AF and TD
            ox = ax.altx(color=c2, label="Step (mT)")#,xmin=0
            TDl=ax.plot(TDtreatment,TDmoment , 'o-', color=c1, markeredgewidth= 0.5 , linewidth= linewidth, markersize=markersize , alpha=0.7,label="TD")
            AFl=ox.plot(AFtreatment,AFmoment, 'o-', color=c2, markeredgewidth= 0.5 , linewidth=linewidth, markersize=markersize , alpha=0.7,label="AFD")
            maxY = max(np.max(TDmoment), np.max(AFmoment))
            ax.format(xlabel="Step (°C)",xcolor=c1, xmin=0, ymin=0)
    elif select=='TD':
        #pur TD
        TDl=ax.plot(TDtreatment,TDmoment , 'o-', color=c1, markeredgewidth= 0.5 , linewidth=linewidth, markersize=markersize , alpha=0.7,label="TD")
        maxY = max(TDmoment)
        ax.format(xlabel="Step (°C)",xcolor=c1, xmin=0, ymin=0)
    elif select=='AF':
        #pur AF
        AFl=ax.plot(AFtreatment,AFmoment , 'o-', color=c2, markeredgewidth= 0.5 , linewidth=linewidth, markersize=markersize , alpha=0.7,label="AFD")
        maxY = max(AFmoment)
        ax.format(xlabel="Step (mT)",xcolor=c1, xmin=0, ymin=0)
    else:
        pass
    ax.legend(handles=[TDl, AFl], loc='upper right', fontsize=8, frameon=False)
    if norm:
        pass
    else:
        y_exp = int(np.floor(np.log10(maxY))) 
        max_exp = y_exp 
        if abs(max_exp) >= 2:
            ax.yaxis.set_major_formatter(precision_formatter(ax, 'y', max_exp))
            ylabeltext = r'Moment ($\times 10^{' + str(max_exp) + r'} ' + momentunit.strip('$') + r'$)'
        else:
        
            ylabeltext = r'Moment ($' + momentunit.strip('$') + r'$)'
    ax.format(ylabel=ylabeltext)     
    return ax



def Multizemplot(MQdirobj,    
    PCA=True,
    Coordinates='g',
    axes=None,
    fishertype=0,
    fisher=False,
    savepath=None,
    joinpdf=None,
    marksize=2,
    linewidth=0.5,
    magnorm=True,
    samplelength=4,
    # textsigma=True,
    abc=True,
    abcloc="l",Iunit='$Am^2$'):
    figures=[]
    siteslist=MQdirobj.siteslist()
    for site in siteslist:
        if isinstance(MQdirobj, Site):
            specimenslist=MQdirobj.specimenslist()
        else:
            specimenslist=MQdirobj[site].specimenslist()
        
        chunk_size = 5

        split_specimen = [specimenslist[i:i + chunk_size] for i in range(0, len(specimenslist), chunk_size)]
        for figpart in range(len(split_specimen)):
            specimenslist=split_specimen[figpart]
            fig,axes=zemplot(rn=len(specimenslist),cn=3)
            for i,specimen in enumerate(specimenslist):
                axZ=axes[i*3]     
                axEQ=axes[i*3+1]  
                axM=axes[i*3+2] 
                sample=specimen[:samplelength]  
                if isinstance(MQdirobj, Site):
                    singlespecimen=MQdirobj[sample][specimen]
                else:
                    
                    singlespecimen=MQdirobj[site][sample][specimen]
                
                szax=NSzplot(singlespecimen,
                            Coordinates=Coordinates,
                            ax=axZ,
                            PCA=PCA,
                            markersize=marksize,
                            linewidth=linewidth,
                            )
                if Coordinates=='s':
                    eqdec=singlespecimen.dec_s
                    eqinc=singlespecimen.inc_s
                    fisherdf=singlespecimen.fisher_s
                elif Coordinates=='g':
                    eqdec=singlespecimen.dec_g
                    eqinc=singlespecimen.inc_g
                    fisherdf=singlespecimen.fisher_g
                elif Coordinates=='t':
                    eqdec=singlespecimen.dec_t
                    eqinc=singlespecimen.inc_t
                    fisherdf=singlespecimen.fisher_t
                seqax=equalareaplot(eqdec,eqinc,a95=None,
                                    fisher=fisher,fisherdf=fisherdf,
                                    ax=axEQ,type=fishertype,
                                    line=True,showticks=False,markersize=marksize,linewidth=linewidth,steptext=singlespecimen.treatmenttext,)
                smax=Smagplot(singlespecimen,select=None,markersize=marksize,
                                    ax=axM,norm=magnorm,)
            fig.format(    
            abc=abc,
            abcloc=abcloc,
            suptitle='Site: '+site,
            )    
            if savepath is not None:
                savepdf=savepath+'/'+site+'_Part'+str(figpart)+'Resetplot.pdf'
                savepng=savepath+'/'+site+'_Part'+str(figpart)+'Resetplot.png'
                fig.savefig(savepdf)
                fig.savefig(savepng)
            else:
                pass
            figures.append(fig)
    if joinpdf is None:
        pass   
    else:
        with PdfPages(joinpdf) as pdf: 
            for fig in figures:
                pdf.savefig(fig)
    return fig,axes
def Multizplot(MQdirobj,    
    PCA=True,
    Coordinates='g',
    axes=None,
    savepath=None,
    joinpdf=None,
    marksize=4,
    linewidth=0.5,
    samplelength=4,
    # textsigma=True,
    abc=True,
    abcloc="l"):
    figures=[]
    siteslist=MQdirobj.siteslist()
    for site in siteslist:
        if isinstance(MQdirobj, Site):
            specimenslist=MQdirobj.specimenslist()
        else:
            specimenslist=MQdirobj[site].specimenslist()
        
        # chunk_size = 5

        # split_specimen = [specimenslist[i:i + chunk_size] for i in range(0, len(specimenslist), chunk_size)]
        # for figpart in range(len(specimenslist)):
        #     specimenslist=specimenslist[figpart]
        fig,axes=mplot(pn=len(specimenslist),cn=2)
        for i,specimen in enumerate(specimenslist):
            ax=axes[i]  
            sample=specimen[:samplelength]  
            if isinstance(MQdirobj, Site):
                singlespecimen=MQdirobj[sample][specimen]
            else:
                
                singlespecimen=MQdirobj[site][sample][specimen]
            
            szax=NSzplot(singlespecimen,
                        Coordinates=Coordinates,
                        ax=ax,
                        PCA=PCA,
                        markersize=marksize,
                        linewidth=linewidth,
                        )
        fig.format(    
        abc=abc,
        abcloc=abcloc,
        suptitle='Site: '+site,
        )   
        if isinstance(axes.number, int):
            subplotpars=1
        else:
            subplotpars  = len(axes)
        pliax=subplotpars-i-1
        if pliax == 0:
            pass
        else:
            for of in range(0, pliax):
                axes[i+1+of].remove()
                #axes[i+1+of].axis('off')     

        if savepath is not None:
            savepdf=savepath+'/'+site+'Zplot.pdf'
            savepng=savepath+'/'+site+'Zplot.png'
            fig.savefig(savepdf)
            fig.savefig(savepng)
        else:
            pass
        figures.append(fig)
    if joinpdf is None:
        pass   
    else:
        with PdfPages(joinpdf) as pdf: 
            for fig in figures:
                pdf.savefig(fig)
    return fig,axes
def mplot(pn=None,rn=None,cn=3,topc='1.3cm',bottomc='0.5cm',rightc='0.5cm',
        leftc='0.5cm',
        wspacec='0.5cm',
        hspacec='0.5cm', 
        spanc=False, 
        sharec=False,
        proj=None):
    if pn==None and rn!=None:
        rown=rn
    elif pn!=None and rn==None:
        rown=math.ceil(pn/cn)
    axesh=30/cn
    fig, axes = upt.subplots(figsize=(30 / 2.54, (rown *axesh+1) / 2.54),ncols=cn, nrows=rown,
        top=topc,
        bottom=bottomc,
        right=rightc,
        left=leftc,
        wspace=wspacec,
        hspace=hspacec, 
        span=spanc, 
        share=sharec,
        proj=proj
        )
    return fig , axes
def zemplot(pn=None,rn=None,cn=3,topc='1cm',bottomc='1cm',rightc='0.7cm',
        leftc='0.5cm',
        wspacec='2cm',
        hspacec='0.5cm', 
        spanc=False, 
        sharec=False):
    if pn==None and rn!=None:
        rown=rn
    elif pn!=None and rn==None:
        rown=math.ceil(pn/cn)
    axesh=30/cn
    fig= upt.figure(figsize=(30 / 2.54, (rown *axesh+2.6) / 2.54),#ncols=cn, nrows=rown,    
        top=topc,
        bottom=bottomc,
        right=rightc,
        left=leftc,
        wspace=wspacec,
        hspace=hspacec, 
        span=spanc, 
        share=sharec
        )
    gs = upt.GridSpec(ncols=cn, nrows=rown,)
    for i in range(1,rn*cn,3):
        ax = fig.subplot(gs[i-1])
        ax = fig.subplot(gs[i], proj='polar')
        ax = fig.subplot(gs[i+1])
        axes=fig.axes
    return fig ,axes
#Basic function



#File convert 
def SQ2GUI(DataPath=None,SavePath=None):
    """
    convert 超导.DAT file to pmag gui stand file as dataframe

    Parameters
    ----------
    PATH 超导文件路径 例如C:\\Users\\wangm\\Desktop\\huizong\\TAISU\\SQDATA treat='A'AFDEMAG treat='T'THERMAL DEMAG

    Returns
    ---------
    dataframe as  pmag gui stand
    """
    if DataPath is None:
        DATpath = os.getcwd()
    else:
        DATpath =DataPath
    datafiles = [file for file in os.listdir(DATpath) if file.endswith('.DAT')]
    SQDAT = pd.DataFrame()
    stdf = pd.DataFrame([])
    for file in datafiles:
        file_path = os.path.join(DATpath, file)
        df = pd.read_csv(file_path, delimiter='\t', on_bad_lines='skip')
        SQDAT = pd.concat([SQDAT, df], ignore_index=True)
    stdf['specimen']=SQDAT.iloc[:,0]
    TYPE1 = (SQDAT.iloc[:, 11] == 'Degauss X, Y, & Z') | (SQDAT.iloc[:, 11].astype(str) == 'None')
    TYPE2=(SQDAT.iloc[:,11]=='Thermal Demag')
    stdf.loc[TYPE1,'treatment'] = SQDAT.loc[TYPE1,'AF X']/10
    stdf.loc[TYPE2,'treatment'] = SQDAT.loc[TYPE2,'Temp C']
    stdf.loc[TYPE1,'treatment_type']='A'
    stdf.loc[TYPE2,'treatment_type']='T'
    
    stdf['moment']=SQDAT.iloc[:,3].astype(float)
    stdf['dec_s']=SQDAT.iloc[:,1].astype(float)
    stdf['inc_s']=SQDAT.iloc[:,2].astype(float)
    stdf=stdf.assign(dec_g=pd.Series(dtype='float'), inc_g=pd.Series(dtype='float'), dec_t=pd.Series(dtype='float'), inc_t=pd.Series(dtype='float'))
    stdf=stdf.reset_index(drop=True)
    if SavePath is None:
        pass
    else:
        stdf.to_excel(SavePath, index=False)
    return stdf
def PSM2GUI(filename,PATH=None,treattype='A'):
    """
    convert PSM.DAT file to pmag gui stand file as dataframe

    Parameters
    ----------
    PATH PSM FILE PATH AS C:\\Users\\wangm\\Desktop\\huizong\\TAISU\\SQDATA treat='A'AFDEMAG treat='T'THERMAL DEMAG
    treattype='A' OR  'T' OR 'B' B means both of thermal and af excist, origin record is true
    Returns
    ---------
    dataframe as  pmag gui stand
    """
    if PATH is None:
        DATpath = os.getcwd()
    else:
        DATpath =PATH
    data=pd.read_csv(DATpath+'\\'+filename, usecols=range(11))
    data['Name'].fillna(method='ffill', inplace=True)
    Nnan_Data= data[data['Measurement'].notna()].iloc[1:,:]
    Nnan_Data.loc[Nnan_Data['Value'].isna(), 'Value']=0 
    stdf = pd.DataFrame([])
    stdf['specimen']=Nnan_Data.loc[:,'Name']
    TYPE1 = (Nnan_Data.iloc[:, 1] == 'AFD') | (Nnan_Data.iloc[:, 1] == 'NRM')
    TYPE2=(Nnan_Data.iloc[:,1]=='TRM')
    stdf.loc[TYPE1,'treatment'] = Nnan_Data.loc[TYPE1,'Value']
    stdf.loc[TYPE2,'treatment'] = Nnan_Data.loc[TYPE2,'Value']
    if treattype=='A':
        stdf.loc[:,'treatment_type']='A'
    elif treattype=='T':
        stdf.loc[:,'treatment_type']='T'
    elif treattype=='B':
        stdf.loc[TYPE1,'treatment_type']='A'
        stdf.loc[TYPE2,'treatment_type']='T'
    stdf['moment']=Nnan_Data.loc[:,'Mag.'].astype(float)
    stdf['dec_s']=Nnan_Data.loc[:,'Dec'].astype(float)
    stdf['inc_s']=Nnan_Data.loc[:,'Inc'].astype(float)
    stdf=stdf.assign(dec_g=pd.Series(dtype='float'), inc_g=pd.Series(dtype='float'), dec_t=pd.Series(dtype='float'), inc_t=pd.Series(dtype='float'))
    stdf=stdf.reset_index(drop=True)
    stdf=stdf.groupby('specimen', group_keys=False).apply(lambda x: x.sort_values(by='treatment'))
    return stdf
def jr62GUI(file,SavePath=None):
    """
    convert jr6 ASCII TEXT file to pmag gui stand file as dataframe

    Parameters
    ----------
    jr6 ASCII TEXT file name treat='A'AFDEMAG treat='T'THERMAL DEMAG
    AFstep smallest heat step 
    
    Returns
    ---------
    dataframe as  pmag gui stand
    """
    stdf = pd.DataFrame([])
    predata = pd.read_csv(file,delim_whitespace=' ')
    stdf['specimen']=predata.iloc[:,0]
    stdf['treatment'] = predata.iloc[:,1].str.replace(r'[a-zA-Z]+','', regex=True)
    stdf['treatment'] = stdf['treatment'].replace('',0)
    stdf['treatment'] = pd.to_numeric(stdf['treatment'])
    #stdf['treatment_type'] = predata.iloc[:,1].str.extract(r'([A-Za-z]+)')
    stdf['treatment_type']=predata.iloc[:,1].str.replace(r'AF|AD', 'A',regex=True).str.replace('TD', 'T').str.replace(r'\d+', '', regex=True).str.replace(r'.', '').str.replace(r'NRM', 'T')
    # if AFstep is None:
    #     pass
    # else:
    #     stdf.loc[stdf['treatment']<AFstep,'treatment_type']='A'
    stdf['moment']=predata.iloc[:,2].astype(float)
    stdf['dec_s']=predata.iloc[:,3].astype(float)
    stdf['inc_s']=predata.iloc[:,4].astype(float)
    stdf=stdf.assign(dec_g=pd.Series(dtype='float'), inc_g=pd.Series(dtype='float'), dec_t=pd.Series(dtype='float'), inc_t=pd.Series(dtype='float'))
    stdf=stdf.reset_index(drop=True)
    stdf = stdf.loc[~((stdf.iloc[:, 0] == 'CY_STD') | (stdf.iloc[:, 0] == 'HOLDER'))]
    if SavePath is None:
        pass
    else:
        stdf.to_excel(SavePath, index=False)
    return stdf





def pca_df(DATAI,Coor='s'):
    scode=DATAI.iloc[:, 0].unique()
    dpca=pd.DataFrame([])
    fisherpcapd=pd.DataFrame([])
    pcapreREC=pd.DataFrame([])
    for spe in scode:
        DATAIs=DATAI.loc[DATAI['specimen']==spe,:].reset_index(drop=True)
        value = DATAIs['PCAP']
        if value.isna().all():
                pass
        else:
                DATAIs['quality']='g'
                PCAmethod=pd.DataFrame(['DE-BFL', 'DE-BFL-A', 'DE-BFL-O','DE-BFP', 'DE-FM'])
                PCA_array= pd.DataFrame(DATAIs['PCAP'].dropna().index)
                for u in range(0,len(PCA_array)-1):
                        starti = PCA_array.iloc[u,0]
                        endi = PCA_array.iloc[u+1,0]
                        if Coor=='g':
                            dataP=DATAIs[['treatment','dec_g','inc_g','moment','quality']].values.tolist()
                        elif Coor=='t':
                            dataP=DATAIs[['treatment','dec_t','inc_t','moment','quality']].values.tolist()
                        else:
                            dataP=DATAIs[['treatment','dec_s','inc_s','moment','quality']].values.tolist()
                        PCAmethodC=PCAmethod.iloc[int(DATAIs['PCAP'].iloc[endi]-1)]
                        pcapre=pd.DataFrame(pmag.domean(dataP,starti,endi,PCAmethodC.to_string(index=False)))
                        fisherpcap=pd.DataFrame(pmag.domean(dataP,starti,endi,'DE-FM'))
                        pcapre['specimen']=spe
                        fisherpcap['specimen']=spe
                        dpca=pd.concat([dpca,pcapre.iloc[0,3:9]]).reset_index(drop=True)
                        fisherpcapd=pd.concat([fisherpcapd,fisherpcap.drop_duplicates(keep='first')]).reset_index(drop=True)
                        pcapreREC=pd.concat([pcapreREC,pcapre]).reset_index(drop=True)
    return pcapreREC,fisherpcapd





def dovandamme(vgp_df):
    """
    Determine the S_b value for VGPs using the Vandamme (1994) method
    for determining cutoff value for "outliers".
    
    Parameters
    ----------
    vgp_df : pandas DataFrame with required column "vgp_lat"
             This should be in the desired coordinate system and assumes one polarity

    Returns
    -------
    vgp_df : after applying cutoff
    cutoff : colatitude cutoff
    S_b : S_b of vgp_df  after applying cutoff
    """
    vgp_df['delta'] = 90.-vgp_df['Plat'].values
    ASD = np.sqrt(np.sum(vgp_df.delta**2)/(vgp_df.shape[0]-1))
    A = 1.8 * ASD + 5.
    delta_max = vgp_df.delta.max()
    if (delta_max<A):
        return vgp_df, A, ASD
    while delta_max > A:
        delta_max = vgp_df.delta.max()
        if delta_max < A:
            return vgp_df, A, ASD
        vgp_df = vgp_df[vgp_df.delta < delta_max]
        ASD = np.sqrt(np.sum(vgp_df.delta**2)/(vgp_df.shape[0]-1))
        A = 1.8 * ASD + 5.
        
def get_sb_df(df, mm97=False):
    """
    Calculates Sf for a dataframe with VGP Lat., and optional Fisher's k, site latitude and N information can be used to correct for within site scatter (McElhinny & McFadden, 1997)

    Parameters
    ----------
    df : Pandas Dataframe with columns
        REQUIRED:
            vgp_lat :  VGP latitude
        ONLY REQUIRED for MM97 correction:
            dir_k : Fisher kappa estimate
            dir_n : number of specimens (samples) per site
            lat : latitude of the site
    mm97 : if True, will do the correction for within site scatter

    Returns
    -------
    Sf : float value for the Sf
    """
    df['delta'] = 90.-df.Plat
    Sp2 = np.sum(df.delta**2)/(df.shape[0]-1)
    if 'k' in df.columns and mm97:        
        ks = df.k
        Ns = df.n
        Ls = np.radians(df.Lat)
        A95s = 140./np.sqrt(ks*Ns)
        Sw2_n = 0.335*(A95s**2)*(2.*(1.+3.*np.sin(Ls)**2) /
                                 (5.-3.*np.sin(Ls)**2))
        return [np.sqrt(Sp2-Sw2_n.mean()),np.sqrt(Sp2)]
    else:
        return [np.sqrt(Sp2),np.sqrt(Sp2)]
    
def scalc_vgp_df(vgp_df, anti=True, rev=0, cutoff=180., kappa=0, n=0, spin=False, v=True, boot=True, mm97=True, nb=1000):
    """
    Calculates Sf for a dataframe with VGP Lat., and optional Fisher's k,
    site latitude and N information can be used to correct for within site
    scatter (McElhinny & McFadden, 1997)

    Parameters
    ----------
    vgp_df : Pandas Dataframe with columns
        REQUIRED:
        vgp_lat :  VGP latitude
        ONLY REQUIRED for MM97 correction:
        dir_k : Fisher kappa estimate
        dir_n_samples : number of samples per site
        lat : latitude of the site
        mm97 : if True, will do the correction for within site scatter
        OPTIONAL:
        boot : if True. do bootstrap
        nb : number of bootstraps, default is 1000
    anti : Boolean
        if True, take antipodes of reverse poles
    spin : Boolean
        if True, transform data to spin axis
    rev : Boolean
        if True, take only reverse poles
    v : Boolean
        if True, filter data with Vandamme (1994) cutoff
    boot : Boolean
        if True, use bootstrap for confidence 95% interval
    mm97 : Boolean
        if True, use McFadden McElhinny 1997 correction for S
    nb : int
        number of bootstrapped pseudosamples for confidence estimate
    verbose : Boolean
        if True, print messages

    Returns
    -------
    N : number of VGPs used in calculation
    S_B : S value
    low : 95% confidence lower bound [0 if boot=0]
    high : 95% confidence upper bound [0 if boot=0]
    cutoff : cutoff used in calculation of  S
    """
    vgp_df['delta'] = 90.-vgp_df.Plat.values
    # filter by cutoff, kappa, and n if desired
    if v: vgp_df = vgp_df[vgp_df.delta <= cutoff]
    if mm97:vgp_df = vgp_df[vgp_df.k >= kappa]
    if n: vgp_df = vgp_df[vgp_df.n >= n]
    if spin:  # do transformation to pole
        Pvgps = vgp_df[['Plon', 'Plat']].values
        ppars = pmag.doprinc(Pvgps)
        Bdirs = np.full((Pvgps.shape[0]), ppars['dec']-180.)
        Bdips = np.full((Pvgps.shape[0]), 90.-ppars['inc'])
        Pvgps = np.column_stack((Pvgps, Bdirs, Bdips))
        lons, lats = pmag.dotilt_V(Pvgps)
        vgp_df['Plon'] = lons
        vgp_df['Plat'] = lats
        vgp_df['delta'] = 90.-vgp_df.Plat
    if anti:
        vgp_rev = vgp_df[vgp_df.Plat < 0]
        vgp_norm = vgp_df[vgp_df.Plat >= 0]
        vgp_anti = vgp_rev
        vgp_anti['Plat'] = -vgp_anti['Plat']
        vgp_anti['Plon'] = (vgp_anti['Plon']-180) % 360
        vgp_df = pd.concat([vgp_norm, vgp_anti], sort=True)
    if rev:
        vgp_df = vgp_df[vgp_df.Plat < 0]  # use only reverse data
    if v:
        vgp_df, cutoff, S_v = dovandamme(vgp_df)  # do vandamme cutoff
    S_B = get_sb_df(vgp_df, mm97=mm97)  # get
    N = vgp_df.shape[0]
    SBs, low, high = [], 0, 0
    if boot:
        for i in range(nb):  # now do bootstrap
            bs_df = vgp_df.sample(n=N, replace=True)
            Sb_bs = get_sb_df(bs_df)
            SBs.append(Sb_bs[0])
        SBs.sort()
        low = SBs[int(.025 * nb)]
        high = SBs[int(.975 * nb)]
    SBdf=pd.DataFrame([])
    SBdf['nVGP']=[N]
    SBdf['VGP_S']=[S_B[1]]
    SBdf['VGP_Sb']=[S_B[0]]
    SBdf['Sb_low']=[low]
    SBdf['Sb_high']=[high]
    SBdf['Cutoff']=[cutoff]
    return SBdf
def dia_vgp(*args):  # new function interface by J.Holmes, SIO, 6/1/2011
    """
    Converts directional data (declination, inclination, alpha95) at a given
    location (Site latitude, Site longitude) to pole position (pole longitude,
    pole latitude, dp, dm).

    Parameters
    ----------
    Takes input as (Dec, Inc, a95, Site latitude, Site longitude)
    Input can be as individual values (5 parameters)
    or
    as a list of lists: [[Dec, Inc, a95, lat, lon],[Dec, Inc, a95, lat, lon]]

    Returns
    -------
    if input is individual values for one pole the return is:
    pole longitude, pole latitude, dp, dm

    if input is list of lists the return is:
    list of pole longitudes, list of pole latitudes, list of dp, list of dm
    
    Examples
    --------
    >>> pmag.dia_vgp(4, 41, 0, 33, -117)
    (41.68629415047637, 79.86259998889103, 0.0, 0.0)
    """
    # test whether arguments are one 2-D list or 5 floats
    if len(args) == 1:  # args comes in as a tuple of multi-dim lists.
        largs = list(args).pop()  # scrap the tuple.
        # reorganize the lists so that we get columns of data in each var.
        (decs, dips, a95s, slats, slongs) = list(zip(*largs))
    else:
        # When args > 1, we are receiving five floats. This usually happens when the invoking script is
        # executed in interactive mode.
        (decs, dips, a95s, slats, slongs) = (args)

    # We send all incoming data to numpy in an array form. Even if it means a
    # 1x1 matrix. That's OKAY. Really.
    (dec, dip, a95, slat, slong) = (np.array(decs), np.array(dips), np.array(a95s),
                                    np.array(slats), np.array(slongs))  # package columns into arrays
    rad = np.pi / 180.  # convert to radians
    dec, dip, a95, slat, slong = dec * rad, dip * \
        rad, a95 * rad, slat * rad, slong * rad
    p = np.arctan2(2.0, np.tan(dip))
    plat = np.arcsin(np.sin(slat) * np.cos(p) +
                     np.cos(slat) * np.sin(p) * np.cos(dec))
    beta = (np.sin(p) * np.sin(dec)) / np.cos(plat)

    # -------------------------------------------------------------------------
    # The deal with "boolmask":
    # We needed a quick way to assign matrix values based on a logic decision, in this case setting boundaries
    # on out-of-bounds conditions. Creating a matrix of boolean values the size of the original matrix and using
    # it to "mask" the assignment solves this problem nicely. The downside to this is that Numpy complains if you
    # attempt to mask a non-matrix, so we have to check for array type and do a normal assignment if the type is
    # scalar. These checks are made before calculating for the rest of the function.
    # -------------------------------------------------------------------------

    boolmask = beta > 1.  # create a mask of boolean values
    if isinstance(beta, np.ndarray):
        beta[boolmask] = 1.  # assigns 1 only to elements that mask TRUE.
    # Numpy gets upset if you try our masking trick with a scalar or a 0-D
    # matrix.
    else:
        if boolmask:
            beta = 1.
    boolmask = beta < -1.
    if isinstance(beta, np.ndarray):
        beta[boolmask] = -1.  # assigns -1 only to elements that mask TRUE.
    else:
        if boolmask:
            beta = -1.

    beta = np.arcsin(beta)
    plong = slong + np.pi - beta
    if (np.cos(p) > np.sin(slat) * np.sin(plat)).any():
        boolmask = (np.cos(p) > (np.sin(slat) * np.sin(plat)))
        if isinstance(plong, np.ndarray):
            plong[boolmask] = (slong + beta)[boolmask]
        else:
            if boolmask:
                plong = slong + beta

    boolmask = (plong < 0)
    if isinstance(plong, np.ndarray):
        plong[boolmask] = plong[boolmask] + 2 * np.pi
    else:
        if boolmask:
            plong = plong + 2 * np.pi

    boolmask = (plong > 2 * np.pi)
    if isinstance(plong, np.ndarray):
        plong[boolmask] = plong[boolmask] - 2 * np.pi
    else:
        if boolmask:
            plong = plong - 2 * np.pi

    dm = np.rad2deg(a95 * (np.sin(p) / np.cos(dip)))
    dp = np.rad2deg(a95 * (pmag.old_div((1 + 3 * (np.cos(p)**2)), 2)))
    plat = np.rad2deg(plat)
    plong = np.rad2deg(plong)

    try:
        VGPdf = pd.DataFrame({'Plon': plong[0], 'Plat': plat[0], 'dp': dp, 'dm': dm})
    except:
        VGPdf = pd.DataFrame({'Plon': plong[0], 'Plat': plat[0], 'dp': dp[0], 'dm': dm[0]})
    return VGPdf#plong.tolist(), plat.tolist(), dp.tolist(), dm.tolist()


def IGRF(input_list, mod='shawqIA', ghfile=""):
    """
    Determine declination, inclination and intensity from a geomagnetic
    field model. The default model used is the IGRF model 
    (http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html) with other
    models available for selection with the available options detailed
    in the mod parameter below.
    
    Parameters:
        input_list : list with format [Date, Altitude, Latitude, Longitude]
            date must be in decimal year format XXXX.XXXX (Common Era)
            altitude is in kilometers
        mod :  desired model
            "" : Use the IGRF13 model by default
            'custom' : use values supplied in ghfile
            or choose from this list
            ['arch3k','cals3k','pfm9k','hfm10k','cals10k.2','cals10k.1b','shadif14','shawq2k','shawqIA','ggf100k']
            where:
                - arch3k (Korte et al., 2009)
                - cals3k (Korte and Constable, 2011)
                - cals10k.1b (Korte et al., 2011)
                - pfm9k  (Nilsson et al., 2014)
                - hfm10k is the hfm.OL1.A1 of Constable et al. (2016)
                - cals10k.2 (Constable et al., 2016)
                - shadif14 (Pavon-Carrasco et al., 2014)
                - shawq2k (Campuzano et al., 2019)
                - shawqIA (Osete et al., 2020)
                - ggk100k (Panovska et al., 2018)[only models from -99950 in 200 year increments allowed)
                - the first four of these models, are constrained to agree
                - with gufm1 (Jackson et al., 2000) for the past four centuries
        gh : path to file with l m g h data

    Returns:
        igrf_array 
            array of magnetic field values (0: dec; 1: inc; 2: intensity (in nT))
    """
    if ghfile != "":
        lmgh = np.loadtxt(ghfile)
        gh = []
        lmgh = np.loadtxt(ghfile).transpose()
        gh.append(lmgh[2][0])
        for i in range(1, lmgh.shape[1]):
            gh.append(lmgh[2][i])
            gh.append(lmgh[3][i])
        if len(gh) == 0:
            print('no valid gh file')
            return
        mod = 'custom'
    if mod == "":
        x, y, z, f = pmag.doigrf(
            input_list[3] % 360., input_list[2], input_list[1], input_list[0])
    elif mod != 'custom':
        x, y, z, f = pmag.doigrf(
            input_list[3] % 360., input_list[2], input_list[1], input_list[0], mod=mod)
    else:
        x, y, z, f = pmag.docustom(
            input_list[3] % 360., input_list[2], input_list[1], gh)

    igrf_array = pmag.cart2dir((x, y, z))
    return igrf_array
def erm(dec_s,inc_s,Azimuth,Dip,Geo2samp=False):
    """
    自学
    Azimuth:北方向顺时针到走向角度
    Dip:向下为负向上为正,与水平方向的夹角
    Ms顺时针旋转,Azimuth,如果dip为正则为顺时针旋转其绝对值,如果dip为负则为逆时针旋转其绝对值，
    参考butl p67 4.3
    -Dip=lisa pmagpy
    90-Dip=palemag.org
    R表示是否是逆推,输入的仍是原始Azimuth和dip
    """
    DIG=pd.DataFrame([])
    if not isinstance(dec_s, pd.DataFrame):
        dec_s = pd.DataFrame(np.transpose([dec_s]))
    if not isinstance(inc_s, pd.DataFrame):
        inc_s = pd.DataFrame(np.transpose([inc_s]))
    DIS=pd.concat([dec_s,inc_s],axis=1)


    
    if not isinstance(Azimuth, pd.DataFrame):
        Azimuth = pd.DataFrame(np.transpose([Azimuth]))
    if not isinstance(Dip, pd.DataFrame):
        Dip = pd.DataFrame(np.transpose([Dip]))
    if len(DIS)!=len(Azimuth):
        Azimuth = pd.DataFrame(np.transpose(np.full(len(DIS), Azimuth.iloc[0])))
        Dip = pd.DataFrame(np.transpose(np.full(len(DIS), Dip.iloc[0])))

    else:
        pass

    step=len(DIS)
    SK_ite = True
    for i in range(step):
        RAzimuth = np.deg2rad(np.array(Azimuth.iloc[i]))
        RDip = np.deg2rad(np.array(Dip.iloc[i]))
        ORIm = pd.DataFrame(pmag.dir2cart(DIS.iloc[i,:]))
        #顺时针从屁股看
        Rz = np.array([[float(np.cos(RAzimuth)),float(-np.sin(RAzimuth)),0],
                       [float(np.sin(RAzimuth)),float(np.cos(RAzimuth)),0],
                       [0,0,1]])
        #逆时针从屁股看
        Ry = np.array([[float(np.cos(RDip)),0,float(np.sin(RDip))],
                   [0,1,0],
                   [float(-np.sin(RDip)),0,float(np.cos(RDip))]])
        if Geo2samp is False:
            RZR=Ry@np.array(np.transpose(ORIm.iloc[0,:]))
            R=np.transpose(pd.DataFrame(pmag.cart2dir(Rz@RZR)))
        elif Geo2samp is True:
            #你旋转倒转矩阵相乘顺序
            RZR=np.array(np.transpose(ORIm.iloc[0,:]))@Rz
            R=np.transpose(pd.DataFrame(pmag.cart2dir(RZR@Ry)))
        
        if SK_ite:
            DIG=R
            SK_ite = False
        else:
            DIG=pd.concat([DIG, R])
    DIG.columns = ['dec_g', 'inc_g','int']
    return DIG.iloc[:,0:2].reset_index(drop=True)
def Titlcor(dec_s,inc_s,Azimuth,Dip):
    """
    利用from scipy.spatial.transform import Rotation as R进行旋转
    Azimuth:北方向顺时针到走向角度
    Dip:向下为负向上为正,与水平方向的夹角
    Ms顺时针旋转,Azimuth,如果dip为正则为顺时针旋转其绝对值,如果dip为负则为逆时针旋转其绝对值，
    参考butl p67 4.3
    Azimuth +180=lisa
    90-Dip=palemag.org
    """
    r = R.from_euler('zyx', [90, 45, 30], degrees=True)
    v = [1, 2, 3]
    r.apply(v)
    DIG=pd.DataFrame([])
    if not isinstance(dec_s, pd.DataFrame):
        dec_s = pd.DataFrame(np.transpose([dec_s]))
    if not isinstance(inc_s, pd.DataFrame):
        inc_s = pd.DataFrame(np.transpose([inc_s]))
    if not isinstance(Azimuth, pd.DataFrame):
        Azimuth = pd.DataFrame(np.transpose([Azimuth]))
    if not isinstance(Dip, pd.DataFrame):
        Dip = pd.DataFrame(np.transpose([Dip]))
    DIS=pd.concat([dec_s,inc_s],axis=1)
    step=len(DIS)
    for i in range(step):
        RAzimuth = np.deg2rad(np.array(Azimuth.iloc[i]))
        RDip = np.deg2rad(np.array(Dip.iloc[i]))
        ORIm = pd.DataFrame(pmag.dir2cart(DIS.iloc[i,:]))
        #顺时针从屁股看
        Rz = np.array([[float(np.cos(RAzimuth)),float(-np.sin(RAzimuth)),0],
                       [float(np.sin(RAzimuth)),float(np.cos(RAzimuth)),0],
                       [0,0,1]])
        #逆时针从屁股看
        Ry = np.array([[float(np.cos(RDip)),0,float(np.sin(RDip))],
                   [0,1,0],
                   [float(-np.sin(RDip)),0,float(np.cos(RDip))]])
        RZR=Ry@np.array(np.transpose(ORIm.iloc[0,:]))
        R=np.transpose(pd.DataFrame(pmag.cart2dir(Rz@RZR)))
        DIG = DIG.append(R)
    return DIG.iloc[:,0:2].reset_index(drop=True)
def circ(dec, inc, alpha,npts=201):
    """
    dec : float
        declination of vector
    dip : float
        dip of vector
    alpha : float
        angle of small circle - 90 if vector is pole to great circle
    npts : int
        number of points on the circle, default 201
    
    Returns
    -------
    D_out, I_out : list
    """
    D_out, I_out = [], []
    dec, inc, alpha = np.radians(dec), np.radians(inc), np.radians(alpha)
    dec1 = dec + np.pi/2.
    isign = 1
    if inc != 0:
        isign = abs(inc)/ inc
    dip1 = inc - isign * (np.pi/ 2.)
    t = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    v = [0, 0, 0]
    t[0][2] = np.cos(dec) * np.cos(inc)
    t[1][2] = np.sin(dec) * np.cos(inc)
    t[2][2] = np.sin(inc)
    t[0][1] = np.cos(dec) * np.cos(dip1)
    t[1][1] = np.sin(dec) * np.cos(dip1)
    t[2][1] = np.sin(dip1)
    t[0][0] = np.cos(dec1)
    t[1][0] = np.sin(dec1)
    t[2][0] = 0
    for i in range(npts):
        psi = float(i) * np.pi / ((npts-1)/2)
        v[0] = np.sin(alpha) * np.cos(psi)
        v[1] = np.sin(alpha) * np.sin(psi)
        v[2] = np.sqrt(abs(1. - v[0]**2 - v[1]**2))
        elli = [0, 0, 0]
        for j in range(3):
            for k in range(3):
                elli[j] = elli[j] + t[j][k] * v[k]
        ellir=np.array(elli).reshape(-1)
        Dir = pmag.cart2dir(ellir)
        D_out.append(Dir[0])
        I_out.append(Dir[1])
    return D_out, I_out
def in_ellipse(center,x_co,y_co,points,scatter='no'):
    """
    Determines if points is located inside ellipse (or circle).
    Parameters
    -----------
    center : list of x and y coordinate of center of ellipse (or circle)
    x_co   : list of x coordinates defining ellipse (or circle)
    y_co   : list of y coordinates defining ellipse (or circle)
    points  : nested list of x and y coordinates of points that will be tested'
    scatter : if 'yes', plots the input points red (if inside ellipse) or blue (outside ellipse)
    
    Return
    -----------
    location: list of 'in' or 'out' statements to indicate if point is inside ellipse

    """
    d=[]
    x_center,y_center=center[0],center[1]   # store x and y coordinate of center
    xs_out,ys_out,xs_in,ys_in=[],[],[],[]
        
    for k in range(len(points)):
        x_point,y_point=points[k][0],points[k][1]
        for i in range(len(x_co)):  # loop through number of ellipse points
         d_to_ellipse_point=np.sqrt((x_co[i]-x_point)**2+(y_co[i]-y_point)**2)    # calculate distance from point to ellipse point
         d.append(d_to_ellipse_point)   # store distances
        d_min=min(d)                    # distance of point to closest ellipse point
        indx=d.index(d_min)                 # get index of closest ellipse point
        x_ell,y_ell=x_co[indx],y_co[indx]   # get coordinates of closest ellipse point

        R1=np.sqrt((x_point-x_center)**2+(y_point-y_center)**2)     # calculate distance from point to center
        R2=np.sqrt((x_ell-x_center)**2+(y_ell-y_center)**2)         # calculate distance from ellipse point to center
        if R1>R2:       # check if point is outside ellipse
            xs_out.append(x_point)
            ys_out.append(y_point)
            location=0
        else:
            xs_in.append(x_point)
            ys_in.append(y_point)
            location=1
        d=[]    # empty list d
    
    if scatter=='yes':
        plt.scatter(xs_out,ys_out,c='blue',s=5)
        plt.scatter(xs_in,ys_in,c='red',s=5)

    return location
def docutoff(vgp_block,angle=45,di_block=0):
    """
    Applies fixed cutoff to set of Virtual Geomagnetic Poles. Default is 45 degrees. 

    Parameters
    -----------
    vgp_block : a nested list of VGPs
    
    Optional Parameters (defaults are used if not specified)
    -----------
    angle : fixed cutoff angle from mean of VGPs
    di_block : a nested list of directions 
    plot : plots VGPs per iteration of cutoff procedure (red if rejected)

    Return
    -----------
    vgp_block : a nested list of VGPs after applying fixed angle cutoff
    vgp_block_rejected : a nested list of rejected VGPs
    di_block : a nested list of accepted directions (if di_block is specified)
    di_block_rejected : a nested list of rejected directions (if di_block is specified)
    """

    vgp_block_rejected=[]
    di_block_rejected=[]

    while True:
        vgp_mean = pmag.fisher_mean(vgp_block)  # calculate Fisher (1953) mean of VGPs
        ref_mean=[vgp_mean['dec'],vgp_mean['inc']]  # store mean longitude and latitude

        cutoffangle=angle   # store fixed cutoff angle, default is 45 degrees
        for i in range(len(vgp_block)): # loop through all VGPs
            angletomean=pmag.angle(vgp_block[i],ref_mean)   # calculate angle between VGP and mean
            if angletomean > cutoffangle:   # check if angle exceeds cutoff angle
                cutoffangle=angletomean     # store maximum angle 
                index=i                     # store index of VGP with maximum angle

        if cutoffangle==angle:  # if no angle exceeds cutoff angle, end algorithm
            break

        vgp_block_rejected.append(vgp_block[index]) # add rejected VGP to list
        del vgp_block[index] # remove rejected VGP from original list of VGPs

        if di_block:
            di_block_rejected.append(di_block[index]) # add rejected direction to list
            del di_block[index] # remove rejected VGP from original list of directions


    if di_block:
        return vgp_block,vgp_block_rejected,di_block,di_block_rejected
    else:
        return vgp_block,vgp_block_rejected
def donewvandamme(vgp_block,di_block=0):
    """
    Applies the Vandamme (1994) to set of Virtual Geomagnetic Poles. Follows recursive method of figure 5 of original publication.

    Parameters
    -----------
    vgp_block : a nested list of VGPs
    
    Optional Parameters (defaults are used if not specified)
    -----------
    plot : plots VGPs per iteration of cutoff procedure (red if rejected)
    di_block : a nested list of directions 

    Return
    -----------
    vgp_block : a nested list of VGPs after applying fixed angle cutoff
    vgp_block_rejected : a nested list of rejected VGPs
    ASD : angular standard deviation after cutoff
    A : cutoff angle
    di_block : a nested list of accepted directions (if di_block is specified)
    di_block_rejected : a nested list of rejected directions (if di_block is specified)
    """
    
    vgp_block_rejected=[]
    di_block_rejected=[]
    
    while True:
        vgp_mean = pmag.fisher_mean(vgp_block)  # calculate Fisher (1953) mean of VGPs
        ref_mean=[vgp_mean['dec'],vgp_mean['inc']]  # store mean longitude and latitude

        maxangle=0
        deltasum=0
        N=len(vgp_block)
        for i in range(N): # loop through all VGPs
            angletomean=pmag.angle(vgp_block[i],ref_mean)   # calculate angle between VGP and mean
            if angletomean[0] > maxangle:   # check if angle exceeds largest angle
                maxangle=angletomean[0]     # store maximum angle 
                index=i                        # store index of VGP with maximum angle
            deltasum+=math.pow(angletomean[0],2) # add angle to sum of squared angles

        # Apply method of Vandamme (1994)
        ASD=math.sqrt(deltasum/(N-1)) # calculate angular standard deviation
        A=1.8 * ASD + 5. # calculate optimum cutoff angle (A)

        if maxangle<A:  # if largest angle does not exceed optimum cutoff angle, end algorithm
            break

        vgp_block_rejected.append(vgp_block[index]) # add rejected VGP to list
        del vgp_block[index] # remove rejected VGP from original list of VGPs

        if di_block:
            di_block_rejected.append(di_block[index]) # add rejected direction to list
            del di_block[index] # remove rejected VGP from original list of directions


    if di_block:
        return vgp_block,vgp_block_rejected,ASD,A,di_block,di_block_rejected
    else:
        return vgp_block,vgp_block_rejected,ASD,A
def di2vgp(di_block,slat=0,slon=0,samp_loc=0):
    """
    Converts directional data (declination, inclination) at a given sampling
    location (Site latitude, Site longitude) to pole position (pole longitude,
    pole latitude)

    Required Parameters
    -----------
    di_block : nested list of direction given as [declination,inclination] 
    slat : sampling latitude
    slon : sampling longitude
    samp_loc : sampling location(s) given as nested list of [site LATITUDE, site LONGITUDE]

    Return
    -----------
    plat : latitude of VGP (or pole) associates with paleomagnetic direction assuming a GAD field
    plon : longitude of VGP (or pole) associates with paleomagnetic direction assuming a GAD field
    """
    # check number of directions
    n_dirs = len(di_block)

    # check if there are multiple sampling locations
    if samp_loc==0:

        input_list=[[di_block[k][0],di_block[k][1],0,slat,slon] for k in range(n_dirs)]
        output_list=pmag.dia_vgp(input_list)
        vgp_block=[[output_list[0][k],output_list[1][k]] for k in range(n_dirs)]

    else:
        n_locs = len(samp_loc)

        if n_dirs!=n_locs:
            print('ERROR: NUMBER OF LOCATIONS ARE NOT EQUAL TO NUMBER OF DIRECTIONS')
        else:
            input_list=[[di_block[k][0],di_block[k][1],0,samp_loc[k][0],samp_loc[k][1]] for k in range(n_dirs)]
            output_list=pmag.dia_vgp(input_list)
            vgp_block=[[output_list[0][k],output_list[1][k]] for k in range(n_dirs)]

    return vgp_block
def fish_dev(k,n=1):
    """
    Generate a random draw from a Fisher distribution with mean declination
    of 0 and inclination of 90 with a specified kappa.
    Parameters
    ----------
    k : kappa (precision parameter) of the distribution
    Returns
    ----------
    dec, inc : declination and inclination of random Fisher distribution draw
               if k is an array, dec, inc are returned as arrays, otherwise, single values
    """
    k = np.array(k)
    R1 = random.random(size=n)
    R2 = random.random(size=n)
    L = np.exp(-2 * k)
    a = R1 * (1 - L) + L
    fac = np.sqrt(-np.log(a)/(2 * k))
    inc = 90. - np.degrees(2 * np.arcsin(fac))
    dec = np.degrees(2 * np.pi * R2)
    print(k,R1,R2,L,a,fac,inc,dec)

    return dec, inc
def fish_VGPs(K=20, n=100, lon=0, lat=90):
    """
    Generates Fisher distributed unit vectors from a specified distribution
    using the pmag.py fshdev and dodirot functions.
    Parameters
    ----------
    k : kappa precision parameter (default is 20)
    n : number of vectors to determine (default is 100)
    lon : mean longitude of distribution (default is 0)
    lat : mean latitude of distribution (default is 90)
    di_block : this function returns a nested list of [lon,lat] as the default
    if di_block = False it will return a list of lon and a list of lat
    Returns
    ---------
    di_block : a nested list of [lon,lat] (default)
    lon,lat : a list of lon and a list of lat (if di_block = False)
    """

    k = np.array(K)
    R1 = random.random(size=n)
    R2 = random.random(size=n)
    L = np.exp(-2 * k)
    a = R1 * (1 - L) + L
    fac = np.sqrt(-np.log(a)/(2 * k))
    inc = 90. - np.degrees(2 * np.arcsin(fac))
    dec = np.degrees(2 * np.pi * R2)

    DipDir, Dip = np.ones(n, dtype=np.float).transpose(
    )*(lon-180.), np.ones(n, dtype=np.float).transpose()*(90.-lat)
    data = np.array([dec, inc, DipDir, Dip]).transpose()
    drot, irot = pmag.dotilt_V(data)
    drot = (drot-180.) % 360.  #
    VGPs = np.column_stack((drot, irot))
    #rot_data = np.column_stack((drot, irot))
    #VGPs = rot_data.tolist()

    # for data in range(n):
    #     lo, la = pmag.fshdev(K)
    #     drot, irot = pmag.dodirot(lo, la, lon, lat)
    #     VGPs.append([drot, irot])
    return VGPs
def gen_pseudopoles(ref_poles,kappas,ns,N_ref_poles,N_test,equal_N=True):
    """
    Generates pseudopole from reference dataset of paleopoles with N and K
    """

    # generate nested list of parametrically sampled VGPs
    nested_VGPs = [fish_VGPs(K=kappas[j],n=ns[j],lon=ref_poles[j][0],lat=ref_poles[j][1]) for j in range(N_ref_poles)]
    # create single list with all simulated VGPs
    sim_VGPs = [nested_VGPs[i][j] for i in range(N_ref_poles) for j in range(ns[i])]
    #print(sim_VGPs)

    if equal_N==True:
        # select random VGPs from dataset
        Inds = np.random.randint(len(sim_VGPs), size=N_test)
        #print(Inds)
        D = np.array(sim_VGPs)
        sample = D[Inds]
    else:
        sample = sim_VGPs

    # compute pseudopole
    polepars = pmag.fisher_mean(sample)
    return [polepars['dec'], polepars['inc']]

def gen_pseudopoles_file(ref_VGPs,N_test):
    """
    Generates pseudopole from reference VGPs stored in csv file
    """

    # select random VGPs from dataset
    Inds = np.random.randint(len(ref_VGPs), size=N_test)
    #print(Inds)
    D = np.array(ref_VGPs)
    sample = D[Inds]

    # compute pseudopole
    polepars = pmag.fisher_mean(sample)
    return [polepars['dec'], polepars['inc']]
    
def par_VGPs(ref_poles,kappas,ns,N_ref_poles):
    """
    Generates nb sets of parametrically sampled VGPs from set of reference poles
    """
    # generate nested list of parametrically sampled VGPs
    nested_VGPs = [fish_VGPs(K=kappas[j],n=ns[j],lon=ref_poles[j][0],lat=ref_poles[j][1]) for j in range(N_ref_poles)]
    #print(nested_VGPs)
    # create single list with all simulated VGPs
    sim_VGPs = [nested_VGPs[i][j] for i in range(N_ref_poles) for j in range(ns[i])]
    #print(sim_VGPs)

    return sim_VGPs

def pseudo_N(DIs, random_seed=None, Nref=None):
    """
    Draw a bootstrap sample of directions returning as many bootstrapped samples
    as in the input directions
    Parameters
    ----------
    DIs : nested list of dec, inc lists (known as a di_block)
    random_seed : set random seed for reproducible number generation (default is None)
    Returns
    -------
    Bootstrap_directions : nested list of dec, inc lists that have been
    bootstrapped resampled
    """
    #print(Nref)
    
    if random_seed != None:
        np.random.seed(random_seed)
    if Nref == None:
        sample_size = len(DIs)
    else:
        sample_size = Nref
    #print(Nref,sample_size)
    Inds = np.random.randint(len(DIs), size=sample_size)
    D = np.array(DIs)
    return D[Inds]

def di_boot(DIs, nb=500, Nref=None):
    """
     returns bootstrap means  for Directional data
     Parameters
     _________________
     DIs : nested list of Dec,Inc pairs
     nb : number of bootstrap pseudosamples
     Returns
    -------
     BDIs:   nested list of bootstrapped mean Dec,Inc pairs
    """
#
# now do bootstrap to collect BDIs  bootstrap means
#
    BDIs = []  # number of bootstraps, list of bootstrap directions
#

    for k in range(nb):  # repeat nb times
        #        if k%50==0:print k,' out of ',nb
        if Nref==None:
            pDIs = pseudo_N(DIs)
        else:
            pDIs = pseudo_N(DIs,Nref=Nref)  # get a pseudosample
        # plt.figure(num=1,figsize=(8,8))
        # ipmag.plot_net(1)   # plot equal area projection
        # ipmag.plot_di(di_block=VGPs,color='k',edge='k',marker='o',markersize=80)
        # ipmag.plot_di(di_block=pDIs,color='red',edge='k',marker='D',markersize=50)
        # plt.show()
        bfpars = pmag.fisher_mean(pDIs)  # get bootstrap mean bootstrap sample
        BDIs.append([bfpars['dec'], bfpars['inc']])
    return BDIs

def par_boot(DIs, kappas, nns, nb=500, Nref=None):
    """
     returns bootstrap means for Directional data,
     calculated using a parametric bootstrap (Tauxe, 2010)

     Parameters
     _________________
     DIs : nested list of Dec,Inc pairs
     nb : number of bootstrap pseudosamples
     Returns
    -------
     BDIs:   nested list of bootstrapped mean Dec,Inc pairs
    """
    #
# now do bootstrap to collect BDIs  bootstrap means
#
    BDIs = []  # number of bootstraps, list of bootstrap directions
    inds = []
#
    for k in range(nb):  # repeat nb times
        sample_size=len(DIs) # determine size of pseudosample
        pDIs=[]
        for i in range(sample_size): 
            ind = np.random.randint(sample_size, size=1) # create random index
            j = ind[0]

            pseudo_dirs = ipmag.fishrot(k=kappas[j],n=nns[j],dec=DIs[j][0],inc=DIs[j][1])
            pseudo_site = ipmag.fisher_mean(di_block=pseudo_dirs)
            pDI = [pseudo_site['dec'],pseudo_site['inc']]
            
            #print(j,DIs[j],kappas[j],nns[j])
            #print(pDI,pseudo_site['k'],pseudo_site['n'])

            pDIs.append(pDI)
            inds.append(j)
    
        bfpars = pmag.fisher_mean(pDIs)  # get bootstrap mean bootstrap sample
        BDIs.append([bfpars['dec'], bfpars['inc']])

        if k % 10 == 0:
            print(k)

    #print(len(inds),sample_size)

    # plt.figure(num=1,figsize=(8,8),dpi=100)
    # plt.hist(inds,bins=25,range=(0,2500),rwidth=0.9,density=False,color='steelblue')
    # plt.show()
    return BDIs

def max_GCD_test(ref_data, test_data, nb=500, boot_test=False):
    """
    Conducts a maximum GCD test
    """
    N_ref = len(ref_data) # store number of VGPs in reference dataset
    N_test = len(test_data) # store number of VGPs in test dataset
    ppoles,ppoles_A95=[],[]
    
    # GENERATE nb PSEUDOPOLES FROM REFERENCE DATASET
    for i in range(nb):
        # SELECT RANDOM VGPS FROM DATASET
        Inds = np.random.randint(N_ref, size=N_test)
        #print(Inds)
        D = np.array(ref_data)
        sample = D[Inds]

        # COMPUTE PSEUDOPOLES
        polepars = pmag.fisher_mean(sample)
        ppoles.append([polepars['dec'], polepars['inc']])
        ppoles_A95.append(polepars['alpha95'])
        #print(synth_poles)

    # COMPUTE REFERENCE MEAN
    VGP_mean = pmag.fisher_mean(ref_data) # compute mean of VGPs
    ref_mean = [VGP_mean['dec'],VGP_mean['inc']] # store reference mean
    
    # ppole_mean = pmag.fisher_mean(ppoles) # compute mean of pseudopoles
    # ref_mean = [ppole_mean['dec'],ppole_mean['inc']] # store reference mean

    # COMPUTE ANGULAR DISTANCE OF PSEUDOPOLES TO REFERENCE MEAN
    D_ppoles=[]
    for j in range(nb):
        ang_distance=pmag.angle(ppoles[j],ref_mean)
        D_ppoles.append(ang_distance[0])

    # COMPUTE STATISTICS
    D_ppoles_mean=mean(D_ppoles)
    D_ppoles_median=median(D_ppoles)
    
    D_ppoles.sort()
    ind_95perc=int(0.95*nb)
    B95 = D_ppoles[ind_95perc]
    print('B95=',B95)

    A95_median=median(ppoles_A95)
    print('Median A95=',A95_median)

    # COMPUTE TEST POLE
    test_mean = pmag.fisher_mean(test_data) # compute mean of VGPs
    test_pole = [test_mean['dec'],test_mean['inc']] # store reference mean
    test_A95 = test_mean['alpha95']
    GCD = pmag.angle(test_pole,ref_mean)
    print('GCD=',GCD) 

    # BOOTSTRAP TEST DATA
    if boot_test==True:
        boot_test_data = di_boot(test_data,nb)
        boot_test_mean = pmag.fisher_mean(boot_test_data)

        D_tpoles=[]
        for j in range(nb):
            angle=pmag.angle(boot_test_data[j],[boot_test_mean['dec'],boot_test_mean['inc']])
            D_tpoles.append(angle[0])

        D_tpoles.sort()
        ind_test_95perc=int(0.95*nb)
        B95_test = D_tpoles[ind_test_95perc]



    if boot_test==True:
        return B95, B95_test, A95_median
    else:
        return B95, test_A95, A95_median

def par_GCD_test(ref_poles, ref_kappas, ref_ns, test_data, N_pp=None, nb=500,boot_test=False, filename=False):
    """
    Conducts a maximum GCD test using a reference dataset of parametrically sampled VGPs
    """
    N_poles = len(ref_poles) # store number of VGPs in reference dataset
    if len(test_data)>1:
        N_test = len(test_data) # store number of VGPs in test dataset
    else:
        N_test = N_pp
    
    # GENERATE nb PSEUDOPOLES FROM REFERENCE DATASET
    if filename==False:
        ppoles = [gen_pseudopoles(ref_poles, ref_kappas, ref_ns, N_poles, N_test) for i in range(nb)]
    else:
        sim_VGPs = genfromtxt(filename, delimiter=',')
        ppoles = [gen_pseudopoles_file(sim_VGPs,N_test) for i in range (nb)]

    # COMPUTE REFERENCE MEAN
    VGP_mean = pmag.fisher_mean(ppoles) # compute mean of VGPs
    ref_mean = [VGP_mean['dec'],VGP_mean['inc']] # store reference mean

    # COMPUTE ANGULAR DISTANCE OF PSEUDOPOLES TO REFERENCE MEAN
    D_ppoles=[]
    for j in range(nb):
        ang_distance=pmag.angle(ppoles[j],ref_mean)
        D_ppoles.append(ang_distance[0])

    # COMPUTE STATISTICS
    D_ppoles.sort()
    ind_95perc=int(0.95*nb)
    B95 = D_ppoles[ind_95perc]

    # COMPUTE TEST POLE
    if len(test_data)>1:
        test_mean = pmag.fisher_mean(test_data) # compute mean of VGPs
        test_pole = [test_mean['dec'],test_mean['inc']] # store reference mean
    else:
        test_pole = test_data

    # BOOTSTRAP TEST DATA
    if boot_test==True:
        boot_test_data = di_boot(test_data,nb)
        boot_test_mean = pmag.fisher_mean(boot_test_data)

        D_tpoles=[]
        for j in range(nb):
            angle=pmag.angle(boot_test_data[j],[boot_test_mean['dec'],boot_test_mean['inc']])
            D_tpoles.append(angle[0])

        D_tpoles.sort()
        ind_test_95perc=int(0.95*nb)
        B95_test=D_tpoles[ind_test_95perc]


    if boot_test==True:
        return B95, B95_test
    else:
        return B95




#作图
def multi_direq_point(ax,DIdf=None,dec=None,inc=None,
                      msn=10,line='N',palpha=0.7,
                      symn=0,LABLOC='best',LABLE='',labelj='Y'):
    ''' 
    等面积点图，不可单独使用
    dfDI df of dec inc
    '''
    symbols = {
    0: {'upper': ['o', 'w', 'k', msn], 'lower': ['o', 'k', 'k', msn]},
    1: {'upper': ['o', 'w', '#0072BD', msn], 'lower': ['o', '#0072BD', '#0072BD', msn]},
    2: {'upper': ['o', 'w', '#D95319', msn], 'lower': ['o', '#D95319', '#D95319', msn]},
    3: {'upper': ['o', 'w', '#EDB120', msn], 'lower': ['o', '#EDB120', '#EDB120', msn]},
    4: {'upper': ['o', 'w', '#7E2F8E', msn], 'lower': ['o', '#7E2F8E', '#7E2F8E', msn]},
    5: {'upper': ['o', 'w', '#77AC30', msn], 'lower': ['o', '#77AC30', '#77AC30', msn]},
    6: {'upper': ['o', 'w', '#4DBEEE', msn], 'lower': ['o', '#4DBEEE', '#4DBEEE', msn]},
    7: {'upper': ['o', 'w', '#A2142F', msn], 'lower': ['o', '#A2142F', '#A2142F', msn]},
    8: {'upper': ['s', 'w', 'k', msn], 'lower': ['s', 'k', 'k', msn]},
    9: {'upper': ['s', 'w', '#0072BD', msn], 'lower': ['s', '#0072BD', '#0072BD', msn]},
    10: {'upper': ['s', 'w', '#D95319', msn], 'lower': ['s', '#D95319', '#D95319', msn]},
    11: {'upper': ['s', 'w', '#EDB120', msn], 'lower': ['s', '#EDB120', '#EDB120', msn]},
    12: {'upper': ['s', 'w', '#7E2F8E', msn], 'lower': ['s', '#7E2F8E', '#7E2F8E', msn]},
    13: {'upper': ['s', 'w', '#77AC30', msn], 'lower': ['s', '#77AC30', '#77AC30', msn]},
    14: {'upper': ['s', 'w', '#4DBEEE', msn], 'lower': ['s', '#4DBEEE', '#4DBEEE', msn]},
    15: {'upper': ['s', 'w', '#A2142F', msn], 'lower': ['s', '#A2142F', '#A2142F', msn]},
    16: {'upper': ['D', 'w', 'k', msn], 'lower': ['D', 'k', 'k', msn]},
    17: {'upper': ['D', 'w', '#0072BD', msn], 'lower': ['D', '#0072BD', '#0072BD', msn]},
    18: {'upper': ['D', 'w', '#D95319', msn], 'lower': ['D', '#D95319', '#D95319', msn]},
    19: {'upper': ['D', 'w', '#EDB120', msn], 'lower': ['D', '#EDB120', '#EDB120', msn]},
    20: {'upper': ['D', 'w', '#7E2F8E', msn], 'lower': ['D', '#7E2F8E', '#7E2F8E', msn]},
    21: {'upper': ['D', 'w', '#77AC30', msn], 'lower': ['D', '#77AC30', '#77AC30', msn]},
    22: {'upper': ['D', 'w', '#4DBEEE', msn], 'lower': ['D', '#4DBEEE', '#4DBEEE', msn]},
    23: {'upper': ['D', 'w', '#A2142F', msn], 'lower': ['D', '#A2142F', '#A2142F', msn]},
    24: {'upper': ['*', 'w', 'k', msn], 'lower': ['*', 'k', 'k', msn]},
    25: {'upper': ['*', 'w', '#0072BD', msn], 'lower': ['*', '#0072BD', '#0072BD', msn]},
    26: {'upper': ['*', 'w', '#D95319', msn], 'lower': ['*', '#D95319', '#D95319', msn]},
    27: {'upper': ['*', 'w', '#EDB120', msn], 'lower': ['*', '#EDB120', '#EDB120', msn]},
    28: {'upper': ['*', 'w', '#7E2F8E', msn], 'lower': ['*', '#7E2F8E', '#7E2F8E', msn]},
    29: {'upper': ['*', 'w', '#77AC30', msn], 'lower': ['*', '#77AC30', '#77AC30', msn]},
    30: {'upper': ['*', 'w', '#4DBEEE', msn], 'lower': ['*', '#4DBEEE', '#4DBEEE', msn]},
    31: {'upper': ['*', 'w', '#A2142F', msn], 'lower': ['*', '#A2142F', '#A2142F', msn]},
    }
    sym=symbols[symn]
    X_down, X_up, Y_down, Y_up ,x,y= [], [], [], [], [], []
    if DIdf is None:
        DIdf=pd.DataFrame([])
        DIdf['dec']=dec
        DIdf['inc']=inc
    for i in range(len(DIdf)):
        XY = pmag.dimap(DIdf.iloc[i,0], DIdf.iloc[i,1])
        x.append(XY[0])
        y.append(XY[1])
        if DIdf.iloc[i,1]>=0:
            X_down.append(XY[0])
            Y_down.append(XY[1])
        else:
            X_up.append(XY[0])
            Y_up.append(XY[1])
    if line=='Y':
        ax.plot(x, y, 'k--',linewidth= 0.5,alpha=0.7, )     
    down=None
    upper=None
    downl=None
    upperl=None
    symb=symbols[0]
    downl=ax.scatter(0, 0, marker=symb['lower'][0],
                c=symb['lower'][1], 
                edgecolor=symb['lower'][2],
                s=5,alpha=0.7,label=' Down')
    
    upperl=ax.scatter(0, 0, marker=symb['upper'][0],
                      c=symb['upper'][1], 
                      edgecolor=symb['upper'][2], 
                      s=5,alpha=0.7,label=' Up')    
    if len(X_down) > 0:
        down=ax.scatter(X_down, Y_down, marker=sym['lower'][0],
                    c=sym['lower'][1], edgecolor=sym['lower'][2],s=sym['lower'][3],alpha=palpha,label=' Down')

    if len(X_up) > 0:
        upper=ax.scatter(X_up, Y_up, marker=sym['upper'][0],
                    c=sym['upper'][1], edgecolor=sym['upper'][2], s=sym['upper'][3],alpha=palpha,label=' Up')
    if labelj=='Y':
        plegend=ax.legend(handles=[downl, upperl],loc=LABLOC,bbox_to_anchor=(0.9, 1.22), ncol=2,frameon=False,markerfirst=False,columnspacing=0.1,handletextpad=0.1,title=LABLE)
    downl.remove()
    upperl.remove()
    return x,y
def Single_direq_pc(ax,averageDIal,
                    symn=0,palpha=1,
                    color='k',linestyle='-',alpha=0.4,lw=1,avemsn=6,labelj='Y',LABLOC='best',LABLE='',):
    ''' 
    等面积alpha95图,不可单独使用
    dfDI df of dec inc
    '''
    averageDIal.reset_index(drop=True, inplace=True)
    dec=averageDIal.loc[0,'specimen_dec'].astype(float)
    inc=averageDIal.loc[0,'specimen_inc'].astype(float)
    alpha95=averageDIal.loc[0,'specimen_alpha95'].astype(float)
    Da95, Ia95 = circ(dec, 
                    inc, 
                    alpha95)
    Xcirc, Ycirc = [], []
    for k in range(len(Da95)):
        XY = pmag.dimap(Da95[k], Ia95[k])
        Xcirc.append(XY[0])
        Ycirc.append(XY[1])
    ax.plot(Xcirc, Ycirc,linestyle=linestyle, c=color, alpha=alpha,lw=lw)
    avex,avey=multi_direq_point(ax,averageDIal.iloc[:,1:3],line='N',symn=symn,msn=avemsn,palpha=palpha,labelj=labelj,LABLOC=LABLOC,LABLE='',)

def single_fisher_area(ax,specimen_pca,alpha=0.4, facecolor='#EDB120', edgecolor='#EDB120'):
    allDIlist=specimen_pca[['specimen_dec','specimen_inc']].values.tolist()
    Site_all_fisher=pd.DataFrame(ipmag.fisher_mean(di_block=allDIlist),index=[0])
    
    Da95, Ia95 = circ(Site_all_fisher.dec.values, 
                    Site_all_fisher.inc.values, 
                    Site_all_fisher.alpha95.values)
    Xcirc, Ycirc = [], []
    for k in range(len(Da95)):
        XY = pmag.dimap(Da95[k], Ia95[k])
        Xcirc.append(XY[0])
        Ycirc.append(XY[1])
    ax.fill(Xcirc, Ycirc, alpha=alpha, facecolor=facecolor, edgecolor=edgecolor, zorder=100000)
    XY = pmag.dimap(Site_all_fisher.dec.values, Site_all_fisher.inc.values)
    nt='n='+str(int(Site_all_fisher.n.values))
    Dt='D='+str(round(float(Site_all_fisher.dec.values),2))
    It='I='+str(round(float(Site_all_fisher.inc.values),2))
    a95t='α_{95}='+str(round(float(Site_all_fisher.alpha95.values),2))
    kt='k='+str(round(float(Site_all_fisher.k.values),2))
    alltext= r'${}${}${}${}${}${}${}${}${}$'.format(nt,'\n',Dt,'\n',It,'\n',a95t,'\n',kt)
    tx,ty=-math.ceil(XY[0])/2,-math.ceil(XY[1])/2
    alltextP=ax.text(tx,0,alltext,ha='center', va='center',color='k',fontsize=6)
    alltextP.set_bbox(dict(facecolor='w',linewidth=0, alpha=0.5, pad=0))
    DI=pd.DataFrame([])

    DI['dec_s']=Site_all_fisher.loc[:,'dec'].astype(float)
    DI['inc_s']=Site_all_fisher.loc[:,'inc'].astype(float)
    avex,avey=multi_direq_point(ax,DI,
                      msn=7,line='N',palpha=1,
                      symn=24,LABLOC='best',LABLE='',labelj='Y')
    
    return Site_all_fisher

def Single_direq(SingleDATAI=None,
                PCApd=None,#PD.DF[NAME,DEC,INC,A95] 
                name=None,
                dec=None,
                inc=None,
                avedec=None,
                aveinc=None,
                alpha95=None,
                plotpca='N',
                line='Y',
                steptext='Y',
                Coor='s',#与DATAI联用
                ax=None,
                symn=0,#选择symbols列表中的样式
                msn=7,
                avemsn=10,
                palpha=0.7,
                avpalpha=1,
                net=True,
                alcolor='k',
                allinestyle='-',
                alalpha=0.4,
                allw=1,
                labelj='Y',
                LABLE='',#指定labelname
                LABLOC='best'):
    """
    directions on equal area plot
    DATAI: pmaggui DF
    site dec inc:采点投图使用,与DATAI不同时,DATAI出现时失效
    symn : set matplotlib symbol 预定义0-31
    STEP是否显示 步骤
    PCApd Dataframe with dec , inc ,a95,color
    return pca df for only single specimen
    """
    matplotlib.rcParams.update({'font.size': 6})
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(18/2.54, 18/2.54))
    if net: 
    	EQNET(ax=ax)
    if SingleDATAI is None:
        if name is None:
            scode = PCApd[['specimen']]
        else:
            scode = name
    else:
    	scode = SingleDATAI['specimen'].unique()
    if len (pd.DataFrame([scode]))==1:   
        DI=pd.DataFrame([])
        aveDI=pd.DataFrame([])
        aveDI['specimen']=scode
        if SingleDATAI is None and dec is not None: 
            steptext='N'
            DI['dec_s']=pd.DataFrame(dec)
            DI['inc_s']=pd.DataFrame(inc)
            x,y=multi_direq_point(ax,DI,line=line,msn=msn,symn=symn,labelj=labelj,LABLOC=LABLOC,LABLE=LABLE,)

        elif SingleDATAI is not None and dec is None:

            if Coor=='g':
                DI = pd.concat([SingleDATAI.iloc[:, 6], SingleDATAI.iloc[:, 7]], axis=1)
            elif Coor=='t':
                DI = pd.concat([SingleDATAI.iloc[:, 8], SingleDATAI.iloc[:, 9]], axis=1)
            elif Coor=='s':
                DI = pd.concat([SingleDATAI.iloc[:, 4], SingleDATAI.iloc[:, 5]], axis=1)
            x,y=multi_direq_point(ax,DI,line=line,msn=msn,symn=symn,palpha=palpha,labelj=labelj,LABLOC=LABLOC,LABLE=LABLE,)    
        if PCApd is None:
            if alpha95 is None: 
                plotpca='N'
            else:
                aveDI['specimen_dec']=pd.DataFrame(avedec)
                aveDI['specimen_inc']=pd.DataFrame(aveinc)
                aveDI['specimen_alpha95']=pd.DataFrame(alpha95)
        else:

            aveDI['specimen_dec']=PCApd['specimen_dec']
            aveDI['specimen_inc']=PCApd['specimen_inc']
            aveDI['specimen_alpha95']=PCApd['specimen_alpha95']
 
        if steptext=='Y':
            GUIstepTEXT(ax,x,y,SingleDATAI,step=3)
        if plotpca=='Y': 
            Single_direq_pc(ax,aveDI,symn=symn+1,avemsn=avemsn,color=alcolor,palpha=avpalpha,
                            linestyle=allinestyle,alpha=alalpha,lw=allw)
        #title =ax.set_title('Down:\u25CB UP:\u25CF', fontsize=7,x=0.5,y=1.1)

    else:
        print('Multiple sample presence use Mzplot')
    return aveDI
def multi_equalarea(DIpd=None,
                    PCApd=None,
                    plotpca='Y',
                    line='N',
                    steptext='N',
                    Coor='s',
                    msn=5,
                    avemsn=5,
                    alcolor='k',
                    allinestyle='-',
                    palpha=1,
                    avpalpha=1,
                    alalpha=0.4,
                    allw=1,
                    symn=0,
                    coln=3,
                    ax=None,
                    fisher='Y',
                    Falpha=0.6, 
                    Ffacecolor='#EDB120', 
                    Fedgecolor='#EDB120',
                    labelj='Y',
                    LABLE='',#指定labelname
                    LABLOC='best'):
    """_summary_

    Args:
        DIpd (DataFrame, optional): point df if generate df can line and step text. Defaults to None.
        PCApd (DataFrame, optional): specimen step demag pca result .
        plotpca (str, optional):plot pca circ or not. Defaults to 'Y'.
        line (str, optional): plot step line or not. Defaults to 'N'.
        steptext (str, optional): plot step text or not. Defaults to 'N'.
        Coor (str, optional): 
            Coor select DIpd plot colum if its int plot site area. t
            he int indicate the number of site from specimen. 
            Defaults to 's'.
                s specimen coor,g geographic coor,t bed corr,int type site plot 
                int is the site str number from PCApd
                loc is genrate one ax plot all from PCApd
        msn (int, optional):size of point. Defaults to 10.
        avemsn (int, optional): size of average point. Defaults to 5.
        alcolor (str, optional): pca circ line color. Defaults to 'k'.
        allinestyle (str, optional): pca line style. Defaults to '-'.
        alalpha (float, optional): pca line alpha. Defaults to 0.4.
        allw (int, optional):pca line width. Defaults to 1.
        symn (int, optional): select symbol style. Defaults to 0.
        coln (int, optional): number of colum. Defaults to 3.
        ax (_type_, optional): _description_. Defaults to None.
        fisher (str, optional): calculate fisher and plot face or not. Defaults to 'Y'.
        Falpha (float, optional): fisher face alpha. Defaults to 0.4.
        Ffacecolor (str, optional): fisher face color. Defaults to 'green'.
        Fedgecolor (str, optional): fisher edge color. Defaults to 'green'.
        labelj (str, optional):label or not. Defaults to 'Y'.
        LABLE (str, optional): label name. Defaults to ''.

    Returns:
        allaveDI: DF[specimen specimen_dec specimen_inc specimen_alpha]
        figure:
    """
    #matplotlib.rcParams.update({'font.size': 6})
    allaveDI=pd.DataFrame([])
    siteallaveDI=pd.DataFrame([])

    if isinstance(Coor, int) or  Coor=='loc':#mean site plot only for pca result            
        if isinstance(Coor, int):
            site=PCApd['specimen'].astype(str).str[:Coor].unique()
        else:
            site=np.array(['empty'])
            PCApd['specimen']='empty'
        scn=len(site)
        if ax is None:
            fig, axes=mpt.xmplot(scn,coln=coln,fsw=18,lcm=0.5,bcm=1,wcm=0.6,hcm=2)
        else:
            axes=ax
            fig=axes.get_figure()
        for i in range(0, scn):
            SPCApd=PCApd.loc[PCApd['specimen'].str.contains(str(site[i])),:].reset_index(drop=True)
            if len(np.shape(axes))==0:
                axi=axes
            else:
                axi = axes[math.ceil((i+1)/coln)-1, i%coln]
            sscn=len(SPCApd)
            for ic in range(0, sscn):
                if ic==0:
                    net=True
                else:
                    net=False
                aveDI=Single_direq(SingleDATAI=None,
                    PCApd=SPCApd.iloc[ic,:],#PD.DF[NAME,DEC,INC,A95] 
                    plotpca=plotpca,
                    line='N',
                    steptext='N',
                    ax=axi,
                    msn=msn,
                    avemsn=avemsn,
                    palpha=palpha,
                    avpalpha=avpalpha,
                    alcolor=alcolor,
                    allinestyle=allinestyle,
                    alalpha=alalpha,allw=allw,
                    symn=symn,#选择symbols列表中的样式
                    net=net,
                    labelj=labelj,
                    LABLE=LABLE,#指定labelname
                    LABLOC=LABLOC,
                    )
            if np.any(site[i]=='empty'):
                pass
            else:
                axi.text(0.1, 1.1,site[i],fontsize=6,transform=axi.transAxes,  ha='left', va='center')
            if fisher=='Y' and len(PCApd)>1:
                single_site_fisher=single_fisher_area(axi,SPCApd,alpha=Falpha, facecolor=Ffacecolor, edgecolor=Fedgecolor)
                single_site_fisher['site']=site[i]
                siteallaveDI=pd.concat([siteallaveDI,single_site_fisher], ignore_index=True)
    else:
        #specimen show plot can plot pca resualt
        scode = DIpd.iloc[:, 0].unique()
        scn = len(scode)
        if ax is None:
            fig, axes=mpt.xmplot(scn,coln=coln,fsw=18,lcm=0.5,bcm=1,wcm=0.6,hcm=2)
        else:
            axes=ax
            fig=axes.get_figure()
        for i in range(0, scn):
            DATAI = DIpd[DIpd.iloc[:, 0] == scode[i]].reset_index(drop=True)
            SPCApd=PCApd.loc[PCApd['specimen'] == scode[i],:].reset_index(drop=True)
            if len(np.shape(axes))==0:
                axi=axes
            else:
                axi = axes[math.ceil((i+1)/coln)-1, i%coln]
            aveDI=Single_direq(SingleDATAI=DATAI,
                PCApd=SPCApd,#PD.DF[NAME,DEC,INC,A95] 
                plotpca=plotpca,
                line=line,
                steptext=steptext,
                Coor=Coor,#与DATAI联用
                ax=axi,
                palpha=palpha,
                avpalpha=avpalpha,
                msn=msn,
                avemsn=avemsn,
                alcolor=alcolor,
                allinestyle=allinestyle,
                alalpha=alalpha,allw=allw,
                symn=symn,#选择symbols列表中的样式
                net=True,
                labelj=labelj,
                LABLE=LABLE,#指定labelname
                LABLOC=LABLOC)
            allaveDI=pd.concat([allaveDI,aveDI], ignore_index=True)
            axi.text(0.1,1.1, scode[i], fontsize=6,transform=axi.transAxes,  ha='left', va='center')
    if len(np.shape(axes))==0:
        pass
    else:
        pliax=coln*math.ceil(scn/coln)-scn
        if pliax == 0:
            pass
        else:
            for of in range(0, pliax):
                axes[math.ceil((scn)/ coln)-1,coln-1-of].axis('off')
        mpt.fignumb(axes,scn,numbposX=0,numbposY=1.1)
    return allaveDI,siteallaveDI,fig

def allsite_equalarea(allsitefisher,ax=None,net=True,linestyle='-',color=None,alpha=0.7,lw=1,
                      averagesitelist=None):
    matplotlib.rcParams.update({'font.size': 6})
    # dec=allsitefisher.dec
    # inc=allsitefisher.inc
    # alpha95=allsitefisher.alpha95
    dec=allsitefisher.loc[0,'dec'].astype(float)
    inc=allsitefisher.loc[0,'inc'].astype(float)
    alpha95=allsitefisher.loc[0,'alpha95'].astype(float)  
    
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(9/2.54, 9/2.54))
    else:
        pass
    if net: 
    	EQNET(ax=ax)
    else:
    	pass
    for i in range(len(allsitefisher)):
        deci=allsitefisher.loc[i,'dec'].astype(float)
        inci=allsitefisher.loc[i,'inc'].astype(float)
        alpha95i=allsitefisher.loc[i,'alpha95'].astype(float)
        Da95, Ia95 = circ(deci, 
                        inci, 
                        alpha95i)
        Xcirc, Ycirc = [], []
        for k in range(len(Da95)):
            XY = pmag.dimap(Da95[k], Ia95[k])
            Xcirc.append(XY[0])
            Ycirc.append(XY[1])
        ax.plot(Xcirc, Ycirc,linestyle=linestyle, c=color, alpha=alpha,lw=lw)
    DI=pd.DataFrame([])

    DI['dec_s']=allsitefisher.loc[:,'dec'].astype(float)
    DI['inc_s']=allsitefisher.loc[:,'inc'].astype(float)
    avex,avey=multi_direq_point(ax,DI,
                      msn=7,line='N',palpha=1,
                      symn=0,LABLOC='best',LABLE='',labelj='Y')
    selectsitedf=allsitefisher.loc[allsitefisher['site'].isin(averagesitelist),:]
    transdf=pd.DataFrame([])
    transdf['specimen_dec']=selectsitedf.loc[:,'dec']
    transdf['specimen_inc']=selectsitedf.loc[:,'inc']

    
    Site_all_fisher=single_fisher_area(ax,transdf,alpha=0.4, facecolor='#EDB120', edgecolor='k')
    return fig,Site_all_fisher












def Single_plot_paleopole(
                    VGPdf,
                    VGPfisher,
                    axes=None,
                    cutoff=None,
                    land_color='darkgrey',
                    land_edge_color='black',
                    ocean_color='lightblue',
                    confcircle=['N','R','F'],
                    Falpha=0,
                    Ffacecolor='r',
                    Labloc=(0.5,0.8),
                    Fedgecolor=None,
                    img_proj = ccrs.Orthographic(0, 90),
                    lat_grid=[0.,30.,60., 80.],
                    lon_grid=[-180., -150., -120.,  -90.,  -60.,-30.,  
                            0.,30.,60., 90., 120., 150.,180.],):
    if axes is None:
        fig = plt.figure(figsize=[18 / 2.54, 18  / 2.54])
        ax = fig.add_subplot(1, 1, 1, projection=img_proj)
    else:
        ax=axes
    
    pointlabeldf = pd.DataFrame({
    'N': ['Normal polarity VGPs'],
    'R': ['Reversal polarity VGPs'],
    'F': ['Fliped polarity VGPs']
    })
    circlelabeldf = pd.DataFrame({
    'N': ['Normal polarity α95'],
    'R': ['Reversal polarity α95'],
    'F': ['Fliped polarity α95']
    })
    pole=confcircle[0]
    vgp=VGPfisher.loc[VGPfisher['polarity']==pole,:]
    ALLP=ax.scatter(vgp.dec,vgp.inc,s=40,marker='p',color='y',edgecolor='k',linewidth=0.5,transform=ccrs.PlateCarree(),label='Mean VGP')
    #for i,pole in enumerate(confcircle):
    
    VGPdf=VGPdf.loc[VGPdf['polarity']==pole,:]
    Plabel=ax.scatter(VGPdf.Plon,VGPdf.Plat,s=6,edgecolor='k',linewidth=0.5,color='#E45F2B',transform=ccrs.PlateCarree(),label=pointlabeldf[pole].values[0])
    
    NDa95, NIa95 = circ(vgp.dec.values, 
                    vgp.inc.values, 
                    vgp.A95.values)
    #ax.fill(Da95, Ia95, alpha=Falpha, facecolor=Ffacecolor, edgecolor=Fedgecolor, zorder=100000,transform=ccrs.PlateCarree())
    circleline=ax.plot(NDa95, NIa95, color=Ffacecolor,lw=1, linestyle='--',zorder=100000,transform=ccrs.PlateCarree(),label=circlelabeldf[pole])
    if cutoff is not None:
        cutline=[]
    else:
        cutd, cuti = circ(0, 
                90, 
                cutoff)
        cutline=ax.plot(cutd, cuti, color='b',lw=1, linestyle='--',zorder=100000,transform=ccrs.PlateCarree(),label='Vcutoff')

    nt='N='+str(int(vgp.n.values))
    Dt='Plat='+str(round(float(vgp.inc.values),2))
    It='Plon='+str(round(float(vgp.dec.values),2))
    a95t='A_{95}='+str(round(float(vgp.A95.values),2))
    kt='K='+str(round(float(vgp.k.values),2))
    alltext= r'${}${}${}${}${}${}${}${}${}$'.format(nt,'\n',Dt,'\n',It,'\n',a95t,'\n',kt)
    Poleinfor=ax.text(0.3,0.8,alltext,ha='left', va='bottom',color='k',fontsize=6,transform=ax.transAxes)
    Poleinfor.set_bbox(dict(facecolor='w',linewidth=0, alpha=0.5, pad=0))
    ax.add_feature(cartopy.feature.OCEAN, zorder=0, facecolor=ocean_color, alpha=0.2)
    ax.add_feature(cartopy.feature.LAND, zorder=0,facecolor=land_color, edgecolor=land_edge_color, alpha=0.2)
    ax.add_feature(cartopy.feature.COASTLINE, zorder=1,color=land_edge_color,alpha=0.2)
    ax.gridlines(xlocs=lon_grid, ylocs=lat_grid,draw_labels=True, linewidth=1,alpha=0.3,color='black', linestyle='dotted')
    ax.set_global()
        #
    polelegend=ax.legend(handles=[Plabel,circleline[0],ALLP],loc=Labloc, ncol=1,frameon=True,columnspacing=0.1,handletextpad=0.1,alignment='left')#
    fig=ax.get_figure()
    return fig




def fisher_mean(dec=None, inc=None, di_block=None, unit_vector=True):

    if di_block is None:
        if dec is None or inc is None:
            raise ValueError("Must provide either di_block or both dec and inc")
        di_block = []
        if unit_vector:
            for n in range(len(dec)):
                di_block.append([dec[n], inc[n], 1.0])
        else:
            for n in range(len(dec)):
                di_block.append([dec[n], inc[n]])
    N = len(di_block)
    fpars = {}
    if N < 2: 
        return {'dec': di_block[0], 
                'inc': di_block[1]}
    X = []
    for direction in di_block:
        cart_coords = dir2cart(direction)
        X.append(cart_coords)
    X = np.array(X)

    Xbar = X.sum(axis=0)

    R = np.linalg.norm(Xbar)
    Xbar = Xbar / R
    direction = cart2dir(Xbar)
    # print(direction)
    try:
        fpars["dec"] = direction[0][0]
        fpars["inc"] = direction[0][1]
    except:
        fpars["dec"] = direction[0]
        fpars["inc"] = direction[1]
    fpars["n"] = N
    fpars["r"] = R
    
    if N != R:
        k = (N - 1.) / (N - R)
        fpars["k"] = k
        csd = 81. / np.sqrt(k)
    else:
        fpars['k'] = 'inf'
        csd = 0.
        
    b = 20.**(1. / (N - 1.)) - 1
    a = 1 - b * (N - R) / R
    if a < -1:
        a = -1
    a95 = np.degrees(np.arccos(a))
    fpars["alpha95"] = a95
    fpars["csd"] = csd
    if a < 0:
        fpars["alpha95"] = 180.0

    return pd.DataFrame(fpars, index=[0])

# def cart2dir(X):
#     R=np.sqrt(X[0]**2+X[1]**2+X[2]**2) # calculate resultant vector length
#     Az=np.degrees(np.arctan2(X[1],X[0]))%360. # calculate declination taking care of correct quadrants (arctan2) and making modulo 360.
#     Pl=np.degrees(np.arcsin(X[2]/R)) # calculate inclination (converting to degrees) #
#     return [Az,Pl,R]
# def dir2cart(Dir):
#     if len(Dir)>2:
#         R=Dir[2]
#     else:
#         R=1.
#     Az,Pl=np.radians(Dir[0]),np.radians(Dir[1])
#     return [R*np.cos(Az)*np.cos(Pl),R*np.sin(Az)*np.cos(Pl),R*np.sin(Pl)]

def flip(di_block, combine=False):
    """
    Determines 'normal' direction along the principle eigenvector, then flips 
    the reverse mode to the antipode.

    Parameters
    ----------
    di_block : nested list of directions
    combine : whether to return directions as one di_block (default is False)
        
    Returns
    -------
    D1 : normal mode
    D2 : flipped reverse mode as two DI blocks
    If combine=True one combined D1 + D2 di_block will be returned
    """
    ppars = doprinc(di_block)  # get principle direction
    if combine:
        D3 = []
    D1, D2 = [], []
    for rec in di_block:
        ang = angle([rec[0], rec[1]], [ppars['dec'], ppars['inc']])
        if ang > 90.:
            d, i = (rec[0] - 180.) % 360., -rec[1]
            D2.append([d, i])
            if combine:
                D3.append([d, i])
        else:
            D1.append([rec[0], rec[1]])
            if combine:
                D3.append([rec[0], rec[1]])
    if combine:
        return D3
    else:
        return D1, D2
def doprinc(data):
    """
    Gets principal components from data in form of a list of [dec,inc,int] data.

    Parameters
    ----------
    data : nested list of dec, inc and optionally intensity vectors

    Returns
    -------
    ppars : dictionary with the principal components
        dec : principal direction declination
        inc : principal direction inclination
        V2dec : intermediate eigenvector declination
        V2inc : intermediate eigenvector inclination
        V3dec : minor eigenvector declination
        V3inc : minor eigenvector inclination
        tau1 : major eigenvalue
        tau2 : intermediate eigenvalue
        tau3 : minor eigenvalue
        N  : number of points
    """
    ppars = {}
    rad = np.pi/180.
    X = dir2cart(data)
    # for rec in data:
    #    dir=[]
    #    for c in rec: dir.append(c)
    #    cart= (dir2cart(dir))
    #    X.append(cart)
#   put in T matrix
#
    T = np.array(Tmatrix(X))
#
#   get sorted evals/evects
#
    t, V = tauV(T)
    Pdir = cart2dir(V[0])
    #ppars['Edir'] = cart2dir(V[1])  # elongation direction - this is V2dec!
    dec, inc = doflip(Pdir[0], Pdir[1])
    ppars['dec'] = dec
    ppars['inc'] = inc
    ppars['N'] = len(data)
    ppars['tau1'] = t[0]
    ppars['tau2'] = t[1]
    ppars['tau3'] = t[2]
    Pdir = cart2dir(V[1])
    dec, inc = doflip(Pdir[0], Pdir[1])
    ppars['V2dec'] = dec
    ppars['V2inc'] = inc
    Pdir = cart2dir(V[2])
    dec, inc = doflip(Pdir[0], Pdir[1])
    ppars['V3dec'] = dec
    ppars['V3inc'] = inc
    return ppars


def dir2cart(d):
    """
    Converts a list or array of vector directions in degrees (declination,
    inclination) to an array of the direction in cartesian coordinates (x,y,z).

    Parameters
    ----------
    d : list or array of [dec,inc] or [dec,inc,intensity]

    Returns
    -------
    cart : array of [x,y,z]

    Examples
    --------
    >>> pmag.dir2cart([200,40,1])
    array([-0.71984631, -0.26200263,  0.64278761])
    
    >>> pmag.dir2cart([200,40])
    array([[-0.719846310392954, -0.262002630229385,  0.642787609686539]])
    
    >>> data = np.array([  [16.0,    43.0, 21620.33],
           [30.5,    53.6, 12922.58],
            [6.9,    33.2, 15780.08],
          [352.5,    40.2, 33947.52], 
          [354.2,    45.1, 19725.45]])
    >>> pmag.dir2cart(data)
    array([[15199.574113612794 ,  4358.407742577491 , 14745.029604010038 ],
       [ 6607.405832448041 ,  3892.0594770716   , 10401.304487835589 ],
       [13108.574245285025 ,  1586.3117853121191,  8640.591471770322 ],
       [25707.154931463603 , -3384.411152593326 , 21911.687763162565 ],
       [13852.355235322588 , -1407.0709331498472, 13972.322052043308 ]])
    """
    rad = np.pi/180.
    ints = np.ones(len(d)).transpose()  # get an array of ones to plug into dec,inc pairs
    d = np.array(d).astype('float')
    if len(d.shape) > 1:  # array of vectors
        decs, incs = d[:, 0] * rad, d[:, 1] * rad
        if d.shape[1] == 3:
            ints = d[:, 2]  # take the given lengths
    else:  # single vector
        decs, incs = np.array(float(d[0])) * rad, np.array(float(d[1])) * rad
        if len(d) == 3:
            ints = np.array(d[2])
        else:
            ints = np.array([1.])
    cart = np.array([ints * np.cos(decs) * np.cos(incs), ints *
                     np.sin(decs) * np.cos(incs), ints * np.sin(incs)])
    cart = np.array([ints * np.cos(decs) * np.cos(incs), ints *
                     np.sin(decs) * np.cos(incs), ints * np.sin(incs)]).transpose()
    return cart
def cart2dir(cart):
    """
    Converts a direction in cartesian coordinates into declinations and inclination.

    Parameters
    ----------
    cart : list of [x,y,z] or list of lists [[x1,y1,z1],[x2,y2,z2]...]

    Returns
    -------
    direction_array : array of [declination, inclination, intensity]

    Examples
    --------
    >>> pmag.cart2dir([0,1,0])
    array([ 90.,   0.,   1.])
    """
    cart = np.array(cart)
    rad = np.pi/180.  # constant to convert degrees to radians
    if len(cart.shape) > 1:
        Xs, Ys, Zs = cart[:, 0], cart[:, 1], cart[:, 2]
    else:  # single vector
        Xs, Ys, Zs = cart[0], cart[1], cart[2]
    if np.iscomplexobj(Xs):
        Xs = Xs.real
    if np.iscomplexobj(Ys):
        Ys = Ys.real
    if np.iscomplexobj(Zs):
        Zs = Zs.real
    Rs = np.sqrt(Xs**2 + Ys**2 + Zs**2)  # calculate resultant vector length
    # calculate declination taking care of correct quadrants (arctan2) and
    # making modulo 360.
    Decs = (np.arctan2(Ys, Xs) / rad) % 360.
    try:
        # calculate inclination (converting to degrees) #
        Incs = np.arcsin(Zs / Rs) / rad
    except:
        print('trouble in cart2dir')  # most likely division by zero somewhere
        return np.zeros(3)
    
    direction_array = np.array([Decs, Incs, Rs]).transpose()  # directions list
    
    return direction_array 
def doflip(dec, inc):
    """
    Flips upper hemisphere data to lower hemisphere.
    
    Parameters
    ----------
    dec : float
        declination
    inc : float
        inclination
    
    Returns 
    -------
    tuple
        containing the flipped declination and inclination 
        
    Examples
    -------
    >>> pmag.doflip(30,-45)
    (210.0, 45)
    """
    if inc < 0:
        inc = -inc
        dec = (dec + 180.) % 360.
    return dec, inc
def angle(D1, D2):
    """
    Calculate the angle between two directions.

    Parameters
    ----------
    D1 : Direction 1 as an array of [declination, inclination] pair or pairs
    D2 : Direction 2 as an array of [declination, inclination] pair or pairs

    Returns
    -------
    angle : single-element array 
        angle between the input directions

    Examples
    --------
    >>> pmag.angle([350.0,10.0],[320.0,20.0])
    array([ 30.59060998])
    
    >>> pmag.angle([[350.0,10.0],[320.0,20.0]],[[345,13],[340,14]])
    array([ 5.744522410794302, 20.026413431433475])
    """
    D1 = np.array(D1)
    if len(D1.shape) > 1:
        D1 = D1[:, 0:2]  # strip off intensity
    else:
        D1 = D1[:2]
    D2 = np.array(D2)
    if len(D2.shape) > 1:
        D2 = D2[:, 0:2]  # strip off intensity
    else:
        D2 = D2[:2]
    X1 = dir2cart(D1)  # convert to cartesian from polar
    X2 = dir2cart(D2)
    angles = []  # set up a list for angles
    for k in range(X1.shape[0]):  # single vector
        angle = np.arccos(np.dot(X1[k], X2[k])) * \
            180. / np.pi  # take the dot product
        angle = angle % 360.
        angles.append(angle)
    return np.array(angles)