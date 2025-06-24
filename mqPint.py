import pickle
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
from adjustText import adjust_text
import proplot as pplt
from scipy.stats import gaussian_kde
pplt.rc.cycle = '538'
matplotlib.rcParams['font.size'] = 6
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['lines.markeredgewidth'] = 0.5  
matplotlib.rcParams['lines.linewidth'] = 0.5 
matplotlib.rcParams['lines.markersize'] = 4
"""
based on proplot
"""
#transdataclass
class transspecimen:
    def __init__(self,specimendata):
        self.specimendataframe=specimendata
    def __getitem__(self):
        return self.specimendataframe
class transsite:
    def __init__(self,sitedata):
        self.specimens = {}
        self.sitedataframe=sitedata
    def add_specimen(self,specimen_name,specimendata):
        self.specimens[specimen_name] = transspecimen(specimendata)
    
    def specimenlist(self):
        return list(self.specimens.keys())
    def __getitem__(self, specimen_name):
        return self.specimens[specimen_name] 
class transdata:
    def __init__(self):
        self.sites = {}


    def add_site(self,site_name,sitedata):
        self.sites[site_name] = transsite(sitedata)
    def __getitem__(self, site_name):
        return self.sites[site_name]  
    def sitelist(self):
        return list(self.sites.keys())     


#load data
def loadBICEP(path):
    with open(path, 'rb') as file:
        thellierData = pickle.load(file)
        plt.close()
    return thellierData


def axesnumber(i,coln,start_r,start_c,skip_r,skip_c,base=0,):
    """
    Generate axes number of axis numbers for the plot.
    base: int
        The number of axes to start with.
    i: int
        The number of axes. for give from 0
    fn: int
        The amount of figures
    coln: int
        colum amount
    start_r : int
        Starting row number.
    strat_c : int
        Starting column number.
    skip_r : int
        Row spacing between axes.
    skip_c : int
        Column spacing between axes.
    Returns axes number for axes[]
    """
    if skip_c>coln:
        raise ValueError("Column skiping is larger than column amount.")
    # 每一行有多少个axes
    ax_onerow=math.ceil((coln-start_c)/(1+skip_c))
    #起始数字
    STARTnumber=coln*start_r
    #fillrow第i个axes所在的行数,第i个axes在第几个周期
    #当i<ax_onerow时fillrow=0 addnumber=i+1
    #考虑i被整除时fillrow应该减一
    fillrow, addnumber = divmod(i+1, ax_onerow)
    #最后一行是累加
    if addnumber==0:
        newrowcol=start_c+ax_onerow+(ax_onerow-1)*skip_c
        addrownumber=coln*(fillrow-1)*(skip_r+1)
    else:
        newrowcol=start_c+addnumber+(addnumber-1)*skip_c
        addrownumber=coln*(fillrow)*(skip_r+1)
    
    axn=base+STARTnumber+addrownumber+newrowcol-1
    return axn

#0.axes 网格建立
def mplot(pn=None,rn=None,cn=3,topc='0.7cm',bottomc='0.7cm',rightc='0.7cm',
        leftc='1cm',
        wspacec='1.3cm',
        hspacec='1.3cm', 
        spanc=False, 
        sharec=False):
    if pn==None and rn!=None:
        rown=rn
    elif pn!=None and rn==None:
        rown=math.ceil(pn/cn)
    axesh=18/cn
    fig, axes = pplt.subplots(figsize=(18 / 2.54, rown *axesh / 2.54),ncols=cn, nrows=rown,
        top=topc,
        bottom=bottomc,
        right=rightc,
        left=leftc,
        wspace=wspacec,
        hspace=hspacec, 
        span=spanc, 
        share=sharec
        )
    return fig , axes
def equallime(allo):
    all=np.vstack((allo, [0,0,0]))
    x=np.ravel(allo[:,0])
    y=np.ravel(allo[:,1:3])
    x0=np.ravel(all[:,0])
    y0=np.ravel(all[:,1:3])
    M_array=np.array([np.min(x),np.min(y),np.max(x),np.max(y)])
    M_array0=np.array([np.min(x0),np.min(y0),np.max(x0),np.max(y0)])
    abs_max = np.abs(M_array[np.argmax(np.abs(M_array))])
    X_absmax=x[np.argmax(np.abs(x))]
    Y_absmax=y[np.argmax(np.abs(y))]
    if X_absmax>=0 and Y_absmax>=0: #第一象限
        OUTer=np.abs((abs_max-np.min(M_array0[0:2]))/10)
        ylm=[np.min(M_array0[0:2])-OUTer, abs_max+OUTer]
        xlm=[np.min(M_array0[0:2])-OUTer, abs_max+OUTer]
    elif X_absmax<=0 and Y_absmax>=0: #第二象限
        OUTer=np.abs((abs_max+np.max(np.abs(M_array0[1:3])))/10)
        ylm=[-np.max(np.abs(M_array0[1:3]))-OUTer,abs_max+OUTer]
        xlm=[-abs_max-OUTer,np.max(np.abs(M_array0[1:3]))+OUTer]
    elif X_absmax<=0 and Y_absmax<=0: #第三象限
        OUTer=np.abs((np.max(M_array0[2:4])+abs_max)/10)
        ylm=[-abs_max-OUTer, np.max(M_array0[2:4])+OUTer]
        xlm=[-abs_max-OUTer, np.max(M_array0[2:4])+OUTer]
    elif X_absmax>=0 and Y_absmax<=0: #第四象限
        OUTer=np.abs((max(np.abs(M_array0[0]),np.abs(M_array0[3]))+abs_max)/10)
        ylm=[-abs_max-OUTer,max(np.abs(M_array0[0]),np.abs(M_array0[3]))+OUTer]
        xlm=[-max(np.abs(M_array0[0]),np.abs(M_array0[3]))-OUTer,abs_max+OUTer]
    return ylm , xlm
         
def Eztypeax(ax,allo,xlabel='E',ylabel=r'$N \bullet UP \square$' ):
    """
    Equal X-axis range for z plot
    """
    # This font supports many Unicode characters
    ylm,xlm=equallime(allo)
    
    X_ticks =ax.get_xticks()[1]-ax.get_xticks()[0]
    # ax.yaxis.set_major_locator(MultipleLocator(X_ticks))
    Xformatted_value = "{:.0e}".format(X_ticks)
    Xmantissa, Xexponente = Xformatted_value.split('e')
    Xexponent =int(Xexponente)
    Xlatex_string = r'{}×$10^{{{}}}$'.format(Xmantissa, Xexponent)
    xlimp= xlm[1]
    ylimp= ylm[1]
    ax.text(0,ylimp, ylabel, ha='center', va='bottom', transform=ax.transData)
    ax.text(xlimp+X_ticks*0.15, -X_ticks*0.15,xlabel,transform=ax.transData)
    ax.format(grid=False,xlim=xlm, ylim=ylm,xloc='zero', yloc='zero',xticklabels=False, yticklabels=False, )

    return Xlatex_string

def SingleMTchangeplot(specimendata,ax=None):
    if ax is None:
        fig, ax = pplt.subplots(figsize=(18/2.54, 18/2.54),ncols=1, nrows=1,share=False)
    else:
        pass
    max_value = np.max(specimendata.IZZIdata)
    ox = ax.alty(ylim=(0,max_value*1.1),label='pTRM/NRM$_0$')
    temprange=specimendata.temps[specimendata.tempbool].astype(int)
    ox.plot(specimendata.temps,specimendata.IZZIdata[:,0],'-o',markerfacecolor='w',c='k',markeredgecolor='black',label='Not Used',zorder=0)
    ax.plot(specimendata.temps,specimendata.IZZIdata[:,1],'-o',markerfacecolor='w',c='k',markeredgecolor='black',label='Not Used',zorder=0)
    try:

        ox.plot(temprange,specimendata.IZZITMdata[:,0],'-o',markerfacecolor='#E45F2B',c='#E45F2B',markeredgecolor='black',zorder=10)
        ax.plot(temprange,specimendata.IZZITMdata[:,1],'-o',markerfacecolor='#2A88FA',c='#2A88FA',markeredgecolor='black',zorder=10)
    except:
        pass 

    ax.axvline(x=specimendata.lowerTemp, color='black',alpha=0.5 ,linestyle='--', label='lowerTemp')  # 绘制竖线
    ax.axvline(x=specimendata.upperTemp, color='black',alpha=0.5 ,linestyle='--', label='upperTemp')  # 绘制竖线
    ax.format(
        xlim=(0,max(specimendata.temps)*1.1), 
        ylim=(0,max_value*1.1),
        xlabel='Temperature(°C)', ylabel='NRM/NRM$_0$'
    )        
    return ax
def araitable(specimendata,ax):
    def colorize(value):
        if value is True:
            return '#2A88FA'
        elif value is False:
            return '#E45F2B'
        else:
            return 'black'
    color_row = [colorize(value) for value in specimendata.critic.iloc[1, :]]
    specimendata.critic.loc[2, :] = np.array(color_row).T
    table_data = [
        [r'$n:$', int(specimendata.critic.int_n_measurements[0]), specimendata.critic.int_n_measurements[2]],
        [r'$n_{pTRM}:$', int(specimendata.critic.int_n_ptrm[0]), specimendata.critic.int_n_ptrm[2]],
        [r'$SCAT:$', 'Pass' if specimendata.critic.int_scat[0] else 'Fail', specimendata.critic.int_scat[2]],
        [r'$FRAC:$', '{:.2f}'.format(specimendata.critic.int_frac[0]), specimendata.critic.int_frac[2]],
        [r'$G_{max}:$', '{:.2f}'.format(specimendata.critic.int_gmax[0]), specimendata.critic.int_gmax[2]],
        [r'$β:$', '{:.2f}'.format(specimendata.critic.int_b_beta[0]), specimendata.critic.int_b_beta[2]],
        [r'$MAD:$', '{:.2f}'.format(specimendata.critic.int_mad_free[0]), specimendata.critic.int_mad_free[2]],
        [r'$DANG:$', '{:.2f}'.format(specimendata.critic.int_dang[0]), specimendata.critic.int_dang[2]],
        [r'$\overrightarrow{k}$:', '{:.2f}'.format(specimendata.critic.int_k_prime[0]), specimendata.critic.int_k_prime[2]]
    ]
    #cell_text = [[entry[1] for entry in table_data]]
    row_labels = [entry[0] for entry in table_data]
    #cell_colours = [[entry[2] for entry in table_data]]
    table = ax.table(cellText=[entry[1:2] for entry in table_data],  # 每个子数组代表一行
                    rowLabels=row_labels,
                    cellColours=[entry[2:3] for entry in table_data],
                    edges='open',
                    cellLoc='center', loc='center', bbox=[0.83, 0.55, 0.15, 0.45])
    for key, cell in table.get_celld().items():
        if key[1] == 0:
            cell.get_text().set_color(cell.get_facecolor())
def single_arai(specimendata,ax=None,
                IZZImarker='-o',
                PTRMmarker='^',
                MDmarker='-s',
                IZcolor='#2A88FA',
                ZIcolor='#E45F2B',):
    """
    #mqpint.site.specimen data or same formed df onto the Arai plot.
    
    """

    if ax is None:
        fig, axes = pplt.subplots(figsize=(18/2.54, 18/2.54),ncols=1, nrows=1,share=False)
    else:
        pass

    for i in range(len(specimendata.ptrm_baseline)):
        step_PTRMs=[specimendata.ptrm_baseline.iloc[i].PTRM/specimendata.NRM0,
                    specimendata.Pcheck.iloc[i].PTRM/specimendata.NRM0,
                    specimendata.Pcheck.iloc[i].PTRM/specimendata.NRM0]
        step_NRMs=[specimendata.ptrm_baseline.iloc[i].NRM/specimendata.NRM0,
                    specimendata.ptrm_baseline.iloc[i].NRM/specimendata.NRM0,
                    specimendata.Pcheck.iloc[i].NRM/specimendata.NRM0]
        ax.plot(step_PTRMs,step_NRMs,'k',lw=1,alpha=0.5)
    ax.plot(specimendata.IZZIdata[:,0],specimendata.IZZIdata[:,1],IZZImarker,markerfacecolor='w',c='k',markeredgecolor='black',label='Not Used')
    ax.plot(specimendata.PTRMdata[:,0],specimendata.PTRMdata[:,1],PTRMmarker,markerfacecolor='w',linestyle='none', markeredgecolor='black',label='PTRM Check')
    try:
        ax.plot(specimendata.MDdata[:,0],specimendata.MDdata[:,1],MDmarker,markerfacecolor='w',linestyle='none',markeredgecolor='black',label='MD Check')
    except:
        pass
    try:
        iz_plot=ax.plot(specimendata.IZdata[:,0],specimendata.IZdata[:,1],IZZImarker,
                        markerfacecolor=IZcolor,linestyle='none',markeredgecolor='black',label='I step',zorder=10)
    except:
        pass 
    try:   
        zi_plot=ax.plot(specimendata.ZIdata[:,0],specimendata.ZIdata[:,1],IZZImarker,
                        markerfacecolor=ZIcolor,linestyle='none',markeredgecolor='black',label='Z step',zorder=10)
    except:
        pass 
    try:    
        pcaplot=ax.plot(specimendata.PCAline[:,0],specimendata.PCAline[:,1],'k--',alpha=0.5,zorder=10)
    except:
        pass 
        #---------------------------------------------------------
    #try:
        #if 3*specimendata.lab_B>specimendata.Brang[1][0]:
    for i in range(specimendata.Nxy):
        ax.plot(specimendata.cx[:, i], specimendata.cy[:, i], c='light teal',alpha=0.2,zorder=-1)
    #     else:
    #         pass
    # except:
    #     pass
    tempstext=[]
    for i,temp in enumerate(specimendata.temps):
        tempstext.append(
            ax.text(specimendata.IZZIdata[i,0],specimendata.IZZIdata[i,1],str(int(temp)))
            )
    araitable(specimendata,ax)
    B_anc_sigma=specimendata.B_anc_sigma
    B_anc_up=specimendata.Banc_PCA+B_anc_sigma
    B_anc_down=specimendata.Banc_PCA-B_anc_sigma
    
    try:
        Btext= r'{}{:.2f}$_{{{:.2f}}}^{{{:.2f}}}${}{}{}{:.2f}$_{{{:.2f}}}^{{{:.2f}}}$ {}{}'.format('$B_{anc-PCA}=$',specimendata.Banc_PCA,B_anc_down,B_anc_up,'$μT$','\n',
                                                        '$B_{anc-BiCEP}=$',specimendata.medB[0],specimendata.Brang[0][0],specimendata.Brang[1][0],'$μT$','\n',
                                                        )
        if 3*specimendata.lab_B<specimendata.Brang[1][0]:
            Btext= r'{}{:.2f}$_{{{:.2f}}}^{{{:.2f}}}$ {}{}{}{}'.format('$B_{anc}=$',specimendata.Banc_PCA,B_anc_down,B_anc_up,'$μT$','\n',
                                                        '$B_{anc-BiCEP}$ Failed','\n',
                                                        )
        else:
            Btext= r'{}{:.2f}$_{{{:.2f}}}^{{{:.2f}}}${}{}{}{:.2f}$_{{{:.2f}}}^{{{:.2f}}}$ {}{}'.format('$B_{anc-PCA}=$',specimendata.Banc_PCA,B_anc_down,B_anc_up,'$μT$','\n',
                                                            '$B_{anc-BiCEP}=$',specimendata.medB[0],specimendata.Brang[0][0],specimendata.Brang[1][0],'$μT$','\n',
                                                            )
    except:
        Btext= r'{}{:.2f}$_{{{:.2f}}}^{{{:.2f}}}$ {}{}{}{}'.format('$B_{anc}=$',specimendata.Banc_PCA,B_anc_down,B_anc_up,'$μT$','\n',
                                                        '$B_{anc-BiCEP}$ Failed','\n',
                                                        )
    

    ax.text(0.02,0.01,Btext,ha='left', va='bottom',transform=ax.transAxes)
    axname='Specimen: '+specimendata.specimen_name+'  '+specimendata.comment+'  '+specimendata.quality
    ax.format(
        title=axname,
        grid=False,
        xlim=(0,max(specimendata.IZZIdata[:,0])*1.1), ylim=(0,max(specimendata.IZZIdata[:,1])*1.1),
        xlabel='pTRM/NRM$_0$', ylabel='NRM/NRM$_0$'
    )
    adjust_text(tempstext,ax=ax)
    
def Szplot(specimendata,ax=None,Iunit='A/m'):
    if ax is None:
        fig, axes = pplt.subplots(figsize=(18/2.54, 18/2.54),ncols=1, nrows=1,share=False)
    else:
        pass
    #---------------------------------------------------------------
    ax.plot(specimendata.Zplotdata[:,0],specimendata.Zplotdata[:,1],
            color='black', linestyle='-', marker='o',markerfacecolor='black', markeredgecolor='black',label='N')#N实心圆
    ax.plot(specimendata.Zplotdata[:,0],specimendata.Zplotdata[:,2],
            color='black', linestyle='-', marker='s',markerfacecolor='white', markeredgecolor='black',label='UP')#UP空心正方形
    ax.plot(specimendata.Zplotdatarange[:,0],specimendata.Zplotdatarange[:,1],'o',markerfacecolor='none',markeredgecolor='#E45F2B',markeredgewidth=0.8,)
    ax.plot(specimendata.Zplotdatarange[:,0],specimendata.Zplotdatarange[:,2],'s',markerfacecolor='none',markeredgecolor='#2A88FA',markeredgewidth=0.8,)
    ax.plot(0, 0, alpha=0)
    #--------------------------------------------------------------     


    ax.annotate("",xy=(specimendata.Varrow_point),\
                xytext=(specimendata.Vtext_point),\
                arrowprops=dict(arrowstyle='simple',mutation_scale=5, linewidth=0.5,color='#E45F2B',alpha=0.5),
                annotation_clip=False)
    ax.annotate("",xy=(specimendata.Harrow_point),\
                xytext=(specimendata.Htext_point),\
                arrowprops=dict(arrowstyle='simple',mutation_scale=5,linewidth=0.5,color='#2A88FA',alpha=0.5),
                annotation_clip=False)

    #--------------------------------------------------------------
    # nrm_value = "{:.2e}".format(specimendata.NRM0)
    # NRM_man, NRM_exp = nrm_value.split('e')
    # NRM_str = r'NRM={}×$10^{{{}}}${}'.format(NRM_man, int(NRM_exp),Iunit)
    steps = []
    if len(specimendata.temps) > 5:
        INDIC = np.arange(1, len(specimendata.temps), int(len(specimendata.temps)/5))
    else:
        INDIC = np.arange(1, len(specimendata.temps),1)
    for t in range(len(INDIC)):
        steps.append(ax.text(np.array(specimendata.Zplotdata)[INDIC[t],0],
                            np.array(specimendata.Zplotdata)[INDIC[t],1],
                            str(specimendata.temps[INDIC[t]]),
                            ha='center', va='center',color='gray',fontsize=6))

    Xax_str=Eztypeax(ax,specimendata.Zplotdata,xlabel='E',ylabel='N UP')
    Scalestr="Scale:{}{}".format(Xax_str,Iunit)
    #figtext=[
    #        ax.text(0.1,1.06, NRM_str, fontsize=6,transform=ax.transAxes,  ha='left', va='center'),
    #        ax.text(0.1,1.02, Scalestr, fontsize=6,transform=ax.transAxes,  ha='left', va='center')
    #            ]
    adjust_text(steps,ax=ax)    
    ax.format(title=Scalestr)
    return ax


def singlesiteintplot(sitedata,axes=None):
    pplt.rc.cycle = 'Paired'
    if axes is None:
        fig, axes = pplt.subplots(figsize=(18/2.54, 7/2.54),
            ncols=3, nrows=1, span=True,sharey='limits',
            wspace=(0, 0),   
        )

    else:
        pass
    axcirc=axes[0]
    histaxes=axes[1]
    pointaxes=axes[2]
    specimenlist=sitedata.specimenslist()
    site_name=sitedata.site_name
    Biceplegendlist=[]
    PCAanclegendlist=[]
    sitepcaint=sitedata.int_abs
    sitepcaintup=sitedata.int_abs+sitedata.int_abs_sigma
    sitepcaintdown=sitedata.int_abs-sitedata.int_abs_sigma
    pointaxes.axhline(sitepcaint,linestyle=':',color='gray',linewidth=1,alpha=0.5)
    pointaxes.axhline(sitepcaintdown,linestyle=':',color='gray',linewidth=1,alpha=0.5)
    pointaxes.axhline(sitepcaintup,linestyle=':',color='gray',linewidth=1,alpha=0.5)
    x = np.arange(len(specimenlist))
    for i, specimen in enumerate(specimenlist):
        try:
            krang=sitedata[specimen].krang
            medk=sitedata[specimen].medk.tolist()
            Brang=sitedata[specimen].Brang
            medB=sitedata[specimen].medB.tolist() 
            BICEPancpoint=axcirc.errorbar(medk, medB,xerr=abs(krang-medk), yerr=abs(Brang-medB),fmt='o', linewidth=1, capsize=3,label=specimen)
            Biceplegendlist.append(BICEPancpoint)
        except:
            pass           
        Banc_PCA=sitedata[specimen].Banc_PCA
        B_anc_sigma=sitedata[specimen].B_anc_sigma
        PCAancpoint=pointaxes.errorbar(x[i], Banc_PCA, yerr=B_anc_sigma,fmt='o', linewidth=1, capsize=3,label=specimen) 
        PCAanclegendlist.append(PCAancpoint)
    pointaxes.format(xlim=(-1, len(specimenlist)), xlocator=np.arange(len(specimenlist)),
        xticklabels=specimenlist,
        )
    #pointaxes.set_xticklabels(specimenlist)
    axcirc.axvline(0,linestyle='--',color='k',linewidth=1,alpha=0.5)
    int_sites=sitedata.data.int_sites.dropna().values
    Blimmin,minB,sig1lB,medB,sig1rB,maxB,Blimmax=np.percentile(int_sites,(1,2.5,31.7,50,68.3,97.5,99))
    axcirc.axhline(minB,linestyle='--',color='k',linewidth=1,alpha=0.5)
    axcirc.axhline(sig1lB,linestyle='--',color='k',linewidth=1,alpha=0.5)
    axcirc.axhline(medB,linestyle='--',color='k',linewidth=1,alpha=0.5)
    axcirc.axhline(sig1rB,linestyle='--',color='k',linewidth=1,alpha=0.5)
    axcirc.axhline(maxB,linestyle='--',color='k',linewidth=1,alpha=0.5)
    x1 = np.repeat(sitedata.mink, 100)
    x2 = np.repeat(sitedata.maxk, 100)
    circx= np.vstack((x1,x2))
    circy = np.vstack((sitedata.minB,sitedata.maxB))
    axcirc.plot(circx,
                circy,color='gray',alpha=0.2)
    axcirc.set_xlabel('k')
    axcirc.legend(Biceplegendlist,ncols=4,loc='top')
    pointaxes.legend(PCAanclegendlist,ncols=4,loc='top')

    histaxes.hist(sitedata.data.int_sites,orientation='horizontal',bins=500,density=True,alpha=0.3)
    kde = gaussian_kde(sitedata.data.int_sites)
    y_values = np.linspace(min(sitedata.data.int_sites), max(sitedata.data.int_sites), 100)
    histaxes.plot(kde(y_values), y_values)
    histaxes.axhline(minB,linestyle='--',color='k',linewidth=1,alpha=0.5)
    histaxes.axhline(sig1lB,linestyle='--',color='k',linewidth=1,alpha=0.5)
    histaxes.axhline(medB,linestyle='--',color='k',linewidth=1,alpha=0.5)
    histaxes.axhline(sig1rB,linestyle='--',color='k',linewidth=1,alpha=0.5)
    histaxes.axhline(maxB,linestyle='--',color='k',linewidth=1,alpha=0.5)
    # for specimen in specimenlist:

    BICEPsiteTEXT = r'{}{}{}{}{}{}{}{}{}{}{:.3f}{}{}{}'.format(
        '$Model:$ ',sitedata.sampletype,'\n',
        '$Sn:$ ',int(sitedata.Nsample),'\n',
        '$N_{eff}:$ ',int(sitedata.Neff),'\n',
        '$R_{hat}:$ ',sitedata.Rhat,'\n',
        sitedata.category,'\n',
    )

    histaxes.text(0.95, 0.05, BICEPsiteTEXT, ha='right', va='bottom', transform=histaxes.transAxes)
    histx=histaxes.get_xlim()
    histaxes.text(np.mean(histx)*2/3,minB,f"{minB:.2f}"+'$uT$',ha='center', va='bottom')
    histaxes.text(np.mean(histx),sig1lB,f"{sig1lB:.2f}"+'$uT$',ha='center', va='bottom')
    histaxes.text(np.mean(histx)*2/3,medB,f"{medB:.2f}"+'$uT$',ha='center', va='bottom')
    histaxes.text(np.mean(histx),sig1rB,f"{sig1rB:.2f}"+'$uT$',ha='center', va='bottom')
    histaxes.text(np.mean(histx)*2/3,maxB,f"{maxB:.2f}"+'$uT$',ha='center', va='bottom')
    Btext= r'{}{:.2f}$_{{{:.2f}}}^{{{:.2f}}}$ {}{}'.format(
        '$B_{anc}=$',sitepcaint,sitepcaintdown,sitepcaintup,'$μT$','\n',
                                                        )
    pointaxes.text(0.95, 0.95, Btext, ha='right', va='top', transform=pointaxes.transAxes)

    # Find the maximum of the start points
    minofmax = min(maxB, sitepcaintup)
    # Find the minimum of the end points
    maxofmin = max(minB, sitepcaintdown)

    # Check if there is an intersection
    if minofmax > maxofmin:
        cirx=axcirc.get_xlim()
        
        pointx=pointaxes.get_xlim()
        axcirc.fill_between(cirx, minofmax, maxofmin, color='y', alpha=0.5, label='Filled Area')
        histaxes.fill_between(histx, minofmax, maxofmin, color='y', alpha=0.5, label='Filled Area')
        pointaxes.fill_between(pointx, minofmax, maxofmin, color='y', alpha=0.5, label='Filled Area')
        title='Site: '+site_name+' $B_{anc}$'+' consistent range'+f' {minofmax:.2f}-{maxofmin:.2f}μT'
    else:
        title='Site: '+site_name+' $B_{anc}$'+'without consistent range'
        pass
    #    pointaxes.errorbar(x, y,xerr, yerr, fmt='s', linewidth=2, capsize=6)
    axes.format(
        ylim=(Blimmin, Blimmax),
        grid=False,
        ylabel='$B_{anc} (uT)$',
        suptitle=title,
    )
    axes[1:3].format(ytickloc='neither',)
    fig=axes.get_figure()
    return fig[0],axes

def Multiarai(pintdata,savepath=None,):

    siteslist=pintdata.siteslist()
    for site in siteslist:
        specimenslist=pintdata[site].specimenslist()
        fig,axes=mplot(pn=len(specimenslist),rn=None,cn=2,topc='1cm',bottomc='1cm',rightc='0.7cm',
            leftc='1cm',
            wspacec='1.3cm',
            hspacec='1.3cm', 
            spanc=False, 
            sharec=False)
        for i,specimen in enumerate(specimenslist):
            ax=axes[i]
            single_arai(pintdata[site][specimen],ax=ax,)
        axes.format(
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
                axes[i+1+of].axis('off')    
        if savepath is not None:
            savepdf=savepath+site+'arai.pdf'
            savepng=savepath+site+'arai.png'
            fig.savefig(savepdf)
            fig.savefig(savepng)
        else:
            pass
    return fig,axes

def AZMTplot(pintdata,savepath='./'):
    siteslist=pintdata.siteslist()
    for site in siteslist:
        specimenslist=pintdata[site].specimenslist()
        rn=2*len(specimenslist)
        cn=3
        AX_array = np.zeros((rn, cn))
        for i in range(len(specimenslist)):
            AX_array[2*i:2*i+2,0:2]= 3*i+1
            AX_array[2*i,2]= 3*i+1+1 
            AX_array[2*i+1,2]= 3*i+2+1
        fig, axes = pplt.subplots(AX_array,figsize=(18 / 2.54, 6*rn/ 2.54) ,
                            span=False, 
                            share=False,
                            )
        for i,specimen in enumerate(specimenslist):
            axarai=axes[i*3]
            single_arai(pintdata[site][specimen],ax=axarai,)
            axz=axes[i*3+1]
            Szplot(pintdata[site][specimen],ax=axz,Iunit='A/m')
            axint=axes[i*3+2]
            SingleMTchangeplot(pintdata[site][specimen],ax=axint)
        axes.format(
            suptitle='Site: '+site,  
            )
        if savepath is not None:
            savepdf=savepath+site+'arai_Z_MT.pdf'
            savepng=savepath+site+'arai_Z_MT.png'
            fig.savefig(savepdf)
            fig.savefig(savepng)
        else:
            pass
    return fig,axes








#excel class
# class Specimen:
#     def __init__(self, specimendata,int_reals,ks,Bs,Single_specimens_critic,Single_comment,specimen_criticrule):
#         self.specimen_name=specimendata.specimen.iloc[0]
#         self.data = specimendata
#         NRM0=specimendata.NRM.iloc[0]
#         self.NRM0=NRM0
#         Blab=specimendata.B_lab.unique()[0]
#         self.lab_B=Blab
#         IZZI=specimendata.loc[(specimendata.steptype=='IZ')|(specimendata.steptype=='ZI'),:]
#         Pcheck=specimendata.loc[specimendata.steptype=='P',:]
#         self.Pcheck=Pcheck
#         Tdata=specimendata.loc[specimendata.steptype=='T',:]
#         tempbool=(np.array(IZZI.temps)>=specimendata.lowerTemp.dropna()[0])&(np.array(IZZI.temps)<=specimendata.upperTemp.dropna()[0])
#         self.tempbool=tempbool
#         self.lowerTemp=specimendata.lowerTemp.dropna()[0]
#         self.upperTemp=specimendata.upperTemp.dropna()[0]
#         IZZIintemprange=IZZI.loc[tempbool,:].reset_index(drop=True)
#         IZZIdata=np.array([IZZI.PTRM/NRM0,IZZI.NRM/NRM0]).T
#         self.IZZITMdata=IZZIdata[tempbool,:]
#         PTRMdata=np.array([Pcheck.PTRM/NRM0,Pcheck.NRM/NRM0]).T
#         MDdata=np.array([Tdata.PTRM/NRM0,Tdata.NRM/NRM0]).T
        
# #------------------------------------Zplot-------------------------------------    
#         self.Zplotdata=IZZI.loc[:,'NRM_x':'NRM_z'].values   
#         self.Zplotdatarange=IZZIintemprange.loc[:,'NRM_x':'NRM_z'].values
#         def PCAarrow(x,y):
#             xstart=np.array([x[0], y[0]])#START is arrow position
#             points = np.hstack([x.reshape(-1, 1), y.reshape(-1, 1)])
#             mean_point = np.mean(points, axis=0)
#             pca = PCA(n_components=1)
#             pca.fit(points)
#             principal_component = pca.components_[0]
#             distances = np.linalg.norm(points - mean_point, axis=1)
#             mean_distance = np.max(distances)
#             line_start = mean_point - principal_component * mean_distance 
#             line_end = mean_point + principal_component * mean_distance 
#             distance_AB = np.linalg.norm(xstart - line_start)
#             distance_AC = np.linalg.norm(xstart - line_end)

#             # 找出距离A最近的点
#             arrow_point = line_start if distance_AB < distance_AC else line_end
#             text_pont = line_start if distance_AB >= distance_AC else line_end
            
#             return arrow_point,text_pont
#         self.Varrow_point,self.Vtext_point=PCAarrow(self.Zplotdatarange[:,0],self.Zplotdatarange[:,1])
#         self.Harrow_point,self.Htext_point=PCAarrow(self.Zplotdatarange[:,0],self.Zplotdatarange[:,2])
# #--------------------------------------------------------------------------        
        
        
#         self.ptrm_baseline=IZZI.loc[IZZI.temp_step.isin(Pcheck.baseline_temp),:]
#         self.temps=np.array(IZZI.temps)
#         self.IZZIdata=IZZIdata
#         self.PTRMdata=PTRMdata
#         self.MDdata=MDdata
#         if specimendata.lowerTemp.dropna()[0]==specimendata.upperTemp.dropna()[0]:
#             pass 
#         else:
#             IZ=IZZIintemprange.loc[IZZIintemprange.steptype=='IZ',:]
#             ZI=IZZIintemprange.loc[IZZIintemprange.steptype=='ZI',:]
            

#             self.IZdata=np.array([IZ.PTRM/NRM0,IZ.NRM/NRM0]).T 
#             self.ZIdata=np.array([ZI.PTRM/NRM0,ZI.NRM/NRM0]).T 
#             PCAdata =np.array([IZZIintemprange.PTRM/NRM0,IZZIintemprange.NRM/NRM0]).T 
#             kbanc, self.Banc_PCA, self.B_anc_sigma,mean=calculate_int(PCAdata,Blab)
#             #Single_specimens_critic.int_abs[0]*1000000
#             #Single_specimens_critic.int_abs[0]/Blab
#             y_intercept=mean[1]-kbanc*mean[0]
#             x_intercept =-y_intercept/kbanc
#             self.PCAline = np.array([[0,y_intercept],[x_intercept,0]])
            
#             #self.Banc_PCA=abs(direction_vector[1] / direction_vector[0])*Blab*1000000
#         try:
#             Circl_index = specimendata.columns.get_loc(0)
#             Circledata = np.array(specimendata.iloc[:, Circl_index:])
#             Nxy=int(Circledata.shape[1]/2)
#             self.Nxy=Nxy
#             self.cx = Circledata[:, :Nxy]
#             self.cy = Circledata[:, Nxy:]
#         except:
#             self.cx=[]
#             self.cy=[]
#         self.int_reals=int_reals
#         self.ks=ks
#         self.Bs=Bs
#         self.krang=np.percentile(ks,(2.5,97.5),axis=0)
#         self.medk=np.median(ks,axis=0)
#         self.Brang=np.percentile(Bs,(2.5,97.5),axis=0)
#         self.medB=np.median(Bs,axis=0)
#         #------------------------------------critic-------------------------------------        
#         if specimen_criticrule is None:
#             pass
#         else:           
#             J1=specimen_criticrule.loc[:,'table_column'].str.split('.', expand=True).values[:,1]
#             J2=specimen_criticrule.loc[:,'criterion_operation'].values
#             J3=specimen_criticrule.loc[:,'criterion_value'].values
#             Judge = np.array([J1,J2,J3]).T

#             for i in range(len(Judge)):
#                 item=str(Judge[i,0])
#                 number=Single_specimens_critic.loc[0,item]
#                 symbol=Judge[i,1]
#                 value=Judge[i,2]
#                 if symbol=='=':
#                     symbol='=='
#                 expression = f"{number} {symbol} {value}"
#                 try:
#                     single_result = eval(expression)
#                     Single_specimens_critic.loc[1,item]=single_result
#                 except:
#                     Single_specimens_critic.loc[1,item]=None
#             self.judge_items=Judge[:,0]
#         self.critic=Single_specimens_critic
#         #------------------------------------------------------------
#         if pd.isna(Single_comment.comment.iloc[0]): 
#             self.comment=''
#         else:
#             self.comment=Single_comment.comment.iloc[0]
#         if pd.isna(Single_comment.quality.iloc[0]): 
#             self.quality=''
#         else:
#             self.quality=Single_comment.quality.iloc[0]
# class Site:
#     def __init__(self, site_name,sitedata,single_sites_critic,site_criticrule):
#         self.site_name = site_name
#         self.data=sitedata
#         self.specimens = {}
#         self.int_abs=single_sites_critic.int_abs[0]
#         self.int_abs_sigma=single_sites_critic.int_abs_sigma[0]

#         self.mink=sitedata.loc[:,'mink'].dropna().values
#         self.maxk=sitedata.loc[:,'maxk'].dropna().values
#         self.minB=sitedata.loc[:,'minB'].dropna().values
#         self.maxB=sitedata.loc[:,'maxB'].dropna().values
#         self.Rhat=sitedata.loc[:,'Rhat'].dropna().values[0]
#         self.Neff=sitedata.loc[:,'Neff'].dropna().values[0]
#         self.category=sitedata.loc[:,'category'].dropna().values[0]
#         self.Nsample=sitedata.loc[:,'Nsample'].dropna().values[0]
#         self.sampletype=sitedata.loc[:,'sampletype'].dropna().values[0]

#         if site_criticrule is None:
#             pass
#         else:
#             J1=site_criticrule.loc[:,'table_column'].str.split('.', expand=True).to_numpy()[:,1].reshape(1, -1)
#             J2=site_criticrule.loc[:,'criterion_operation'].to_numpy().reshape(1, -1)
#             J3=site_criticrule.loc[:,'criterion_value'].to_numpy().reshape(1, -1)
#             Judge = np.concatenate((J1,
#                             J2,
#                             J3), axis=0).T
#             for i in range(len(Judge)):
#                 item=str(Judge[i,0])
#                 number=single_sites_critic.loc[0,item]
#                 symbol=Judge[i,1]
#                 value=Judge[i,2]
#                 if symbol=='=':
#                     symbol='=='
#                 expression = f"{number} {symbol} {value}"
#                 try:
#                     single_result = eval(expression)
#                     single_sites_critic.loc[1,item]=single_result
#                 except:
#                     single_sites_critic.loc[1,item]=None
#             self.judge_items=Judge[:,0]
#         self.critic=single_sites_critic
#     def __getitem__(self, specimen_name):
#         return self.specimens[specimen_name]
#     def add_specimen(self, specimen_name,data,sitedata,single_specimens_critic,single_comment,specimen_criticrule):
#         specimenindex=sitedata[sitedata['Specimen'] == specimen_name].index.astype(str)
#         intcol='int_reals'+specimenindex
#         int_reals=sitedata[intcol]
#         kscol='ks'+specimenindex
#         ks=sitedata[kscol]
#         Bscol='Bs'+specimenindex
#         Bs=sitedata[Bscol]
#         self.specimens[specimen_name] = Specimen(data,int_reals,ks,Bs,single_specimens_critic,single_comment,specimen_criticrule)
#     def specimenslist(self):
#         return list(self.specimens.keys())


# class mqpintdata:
#     """
#     _summary_
#     BiCEP  data xlsx contains "thellierData" exported from BiCEP notebook export 
#     site result data xlsx contains "critic" and "site" exported from pamg_gui 
#     specimen result data xlsx contains "critic" and "specimen" exported from pamg_gui
#     comment xlsx contains 'comment' inclueding the comment of each specimen and use or not made by the user
#     crictical rule xlsx contains "criticrule" for site and specimen exported from pamg_gui
#     """
#     def __init__(self, Datapath='./'):
#         self.sites = {}
#         self.load_data(Datapath)

#     def load_data(self, Datapath):
#         Sitefilename = [file for file in os.listdir(Datapath) if file.endswith('.xlsx') and 'thellierData' in file]
#         sites_criticfilename=[file for file in os.listdir(Datapath) if file.endswith('.xlsx') and 'critic' in file and 'site' in file]
#         sites_critic=pd.read_excel(os.path.join(Datapath, sites_criticfilename[0]), skiprows=1).dropna(subset=['int_abs'])     
#         specimens_criticfilename=[file for file in os.listdir(Datapath) if file.endswith('.xlsx') and 'critic' in file and 'specimen' in file]
#         specimens_critic=pd.read_excel(os.path.join(Datapath, specimens_criticfilename[0]), skiprows=1).dropna(subset=['int_n_measurements'])
#         commentfilename=[file for file in os.listdir(Datapath) if file.endswith('.xlsx') and 'comment' in file]
#         comment=pd.read_excel(os.path.join(Datapath, commentfilename[0]))
#         criticrulefilename=[file for file in os.listdir(Datapath) if file.endswith('.xlsx') and 'criticrule' in file]

#         try:
#             criticrule=pd.read_excel(os.path.join(Datapath, criticrulefilename[0]),skiprows=1)
#             specimen_criticrule=criticrule.loc[criticrule.loc[:,'table_column'].str.split('.', expand=True).values[:,0]=='specimens',:]
#             site_criticrule=criticrule.loc[
#                 criticrule.loc[:,'table_column'].str.split('.', expand=True).values[:,0]=='sites',:]
#         except:
#             criticrule=None
#             specimen_criticrule=None
#             site_criticrule=None
        
#         for filename in Sitefilename:
#             file_path = os.path.join(Datapath, filename)
#             sheet_names = pd.ExcelFile(file_path).sheet_names
#             site_name = str(filename.split('_')[1].split('.')[0]) 
#             sitedata=pd.read_excel(file_path, sheet_name=sheet_names[0])
#             single_sites_critic=sites_critic.loc[sites_critic['site']==site_name,:].reset_index(drop=True)
#             site = Site(site_name,sitedata,single_sites_critic,site_criticrule)
#             for sheet_name in sheet_names[1:]:
#                 data = pd.read_excel(file_path, sheet_name=sheet_name)
#                 single_specimens_critic=specimens_critic.loc[specimens_critic.loc[:,'specimen']==sheet_name,:].reset_index(drop=True)
#                 single_comment=comment[comment['specimen']==sheet_name]
#                 site.add_specimen(sheet_name,data,sitedata,single_specimens_critic,single_comment,specimen_criticrule)
#             self.sites[site_name] = site
#     def siteslist(self):
#         return list(self.sites.keys())
#     def __getitem__(self, site_name):
#         return self.sites[site_name]
    
    
# def calculate_int(selectarray,lab_dc_field):
#     n=len(selectarray)
#     mean=np.mean(selectarray, axis=0)
#     x_select=selectarray[:,0]
#     y_select=selectarray[:,1]
#     x_mean=mean[0]
#     y_mean=mean[1]
#     x_err = x_select - x_mean
#     y_err = y_select - y_mean
#     k = -1* np.sqrt(sum(y_err**2)/sum(x_err**2))   # averaged slope

#     b = np.sign(sum(x_err * y_err)) * np.std(y_select, ddof=1)/np.std(x_select, ddof=1) # ddof is degrees of freedom
#     if b == 0:
#         k = 1e-10
#     else:
#         k = b

#     k_sigma= np.sqrt( (2 * sum(y_err**2) - 2*k* sum(x_err*y_err))/ ( (n-2) * sum(x_err**2) )) 
#     if k_sigma == 0: # prevent divide by zero
#         k_sigma = 1e-10

#     B_lab = lab_dc_field 
#     B_anc = abs(k) * B_lab # in microtesla
#     B_anc_sigma = k_sigma * B_lab

#     return k, B_anc, B_anc_sigma,mean