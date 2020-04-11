import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
plt.rcParams['font.size']=16

def initskip_olg(ttype):
    if ttype=='terms': wgrep="'LIST OF TERMS WITH A WEIGHTED MEAN OVER THE FINE STRUCTURE'"
    if ttype=='E1': wgrep="'E1-DATA'"
    if ttype=='E2': wgrep="'E2/M1-DATA'"
    if ttype=='E3': wgrep="'E3/M2-DATA'"
    wsys="grep -n -i "+wgrep+" olg > dum.txt"
    os.system(wsys)
    dumline=pd.read_csv("dum.txt",header=None,sep='\s+')
    nrows=len(dumline)
    iline=dumline.loc[nrows-1][0]
    os.system("rm dum.txt")
    return int(iline[0:-1])

def determine_parity(df,ip):
    df.insert(ip,'P',-1)
    odd=df.index[df.loc[:]['2S+1']<0].tolist()
    even=df.index[df.loc[:]['2S+1']>0].tolist()
    for i in odd: df.at[i,'P']=1
    for i in even: df.at[i,'P']=0
    return

def termsymbol_to_quantnumber(chterm):
    chsq=chterm[0]
    chlq=chterm[1]
    # multiplicity 2S+1
    sq=int(chsq)
    # angular momenta L
    if chlq=='S': lq=0
    if chlq=='P': lq=1
    if chlq=='D': lq=2
    if chlq=='F': lq=3
    if chlq=='G': lq=4
    if chlq=='H': lq=5
    if chlq=='I': lq=6
    if chlq=='K': lq=7
    # parity
    if (lq%2)==0: pq=0
    if (lq%2)!=0: pq=1
    return sq,lq,pq

def quantnumber_to_termsymbol(sq,lq,pq):
    # multiplicity 2S+1
    chsq=str(int(sq))
    # angular momenta L
    if lq==0: chlq='S'
    if lq==1: chlq='P'
    if lq==2: chlq='D'
    if lq==3: chlq='F'
    if lq==4: chlq='G'
    if lq==5: chlq='H'
    if lq==6: chlq='I'
    if lq==7: chlq='K'
    # parity
    if pq==0: chpq=''
    if pq!=0: chpq='*'
    chterm=chsq+chlq+chpq
    return chterm

def conv_Aki2fik(ttype,gi,gk,Eki,Aki):
    #
    # convert Aki to fik using SI units
    #
    #    Eki : rydbergs
    #    Aki : s^{-1} 
    #
    IP=13.6056923           # eV/ry
    h=4.135667731E-15       # eV*s
    h_J=6.62607E-34         # J*s
    c=299792458             # m/s
    hc=h*c                  # eV*m
    me=9.10938E-11          # kg
    e=1.60218E-19           # C
    eps0=8.8541878128E-12   # F/m and F= C^2 J^-1 m^-1 
    convJ2eV=6.24E+18       # eV/J
    convm2A=1.0E+10         # Angstrom
    convm2nm=1.0E+9         # nm
    convEh2J=4.35974417E-18 # J/a.u. (hartree)
    convry2J=convEh2J/2     # J/ry
    lam=hc/(Eki*IP)         # eV*m/eV = m
    lamAng=lam*convm2A
    if ttype=='E1':
        E1_f2A=2*np.pi*e**2/(me*c*eps0)*1.0E20
        fik=lam**2*gk*Aki/(E1_f2A*gi)
    if ttype!='E1':
        fik=None
    return fik

def drop_doublexcited(EXtran,pd_droplev):
    #
    # drop double excited levels except 3p^2
    #
    levsk=pd_droplev.loc[:]['K'].tolist()
    for i in levsk:
        drop_k=EXtran.index[EXtran.loc[:]['K']==i].tolist()
        EXtran.drop(drop_k,axis=0,inplace=True)
    EXtran.reset_index(drop=True,inplace=True)
    return

def absolute_Aki(ttype,EXtran):
    #
    # take absolute value of electric and magnetic transition elements
    # sum electric and magnetic einstein coefficients
    #
    if ttype=='E1': Aki_cols=['A(EK)*SEC']
    if ttype!='E1': Aki_cols=['A(EK)*SEC','A(MK)*SEC']
    for i in Aki_cols: EXtran[i]=abs(EXtran[i])
    EXtran['Aki']=sum(EXtran[i] for i in Aki_cols)
    return

def insert_LVdata(EXtran):
    #
    # use k and kp index to insert lv index
    # include statistic weight 
    # replace Eki with observed values
    #
    colin=['LV_k','gk','LV_i','gi','Eki']
    valin=[-999,-999,-999,-999,-999.999]
    idxin=[3,4,5,6,7]
    for i in range(5): EXtran.insert(idxin[i],colin[i],valin[i])
    ntran=len(EXtran)
    for ii in range(ntran):
        dummy_k=as_levs.loc[as_levs['K']==EXtran.loc[ii]['K']]
        dummy_i=as_levs.loc[as_levs['K']==EXtran.loc[ii]['KP']]
        Ek=dummy_k['Energy'].tolist()[0]
        Ei=dummy_i['Energy'].tolist()[0]
        EXtran.at[ii,'LV_k']=dummy_k['LV'].tolist()[0]
        EXtran.at[ii,'LV_i']=dummy_i['LV'].tolist()[0]
        EXtran.at[ii,'gk']=dummy_k['2J'].tolist()[0]+1
        EXtran.at[ii,'gi']=dummy_i['2J'].tolist()[0]+1
        EXtran.at[ii,'Eki']=Ek-Ei
    return

def compute_fik(ttype,EXtran):
    #
    # compute absortion oscillator strength and weigthed fik for allowed transitions only
    #
    EXtran['fik']=conv_Aki2fik(ttype,EXtran['gi'],EXtran['gk'],EXtran['Eki'],EXtran['Aki'])
    if ttype=='E1': EXtran['gf']=EXtran['gi']*EXtran['fik']
    if ttype!='E1': EXtran['gf']=None
    return

def drop_Ncols(EXtran):
    #
    # drop unnecesary columns
    #
    dropcols=[]
    cols=['E1-DATA','E2/M1-DATA','E3/M2-DATA','K','KP','A(EK)*SEC','A(MK)*SEC','Eki']
    for i in cols:
        if i in EXtran.columns: dropcols.append(i)
    EXtran.drop(dropcols,axis=1,inplace=True)
    return

def transformsymb_levs(ttype,as_X,nist_cfgs):
    db_X=as_X.copy()
    # move ttype column to front
    X_list=db_X.loc[:][ttype].tolist()
    coldrop=[ttype]
    if ttype=='LV': coldrop.append('T')
    db_X.drop(coldrop,axis=1,inplace=True)
    db_X.insert(0,ttype,X_list)
    # insert columns for configuration and term symbols
    db_X.insert(2,'Conf','x')
    db_X.insert(3,'Term','x')
    # fill configuration and term columns with symbols
    ndata=len(db_X)
    for i in range(ndata):
        sq=db_X.iloc[i]['2S+1']
        lq=db_X.iloc[i]['L']
        pq=db_X.iloc[i]['P']
        icf=db_X.iloc[i]['CF']
        term=quantnumber_to_termsymbol(sq,lq,pq)
        cf=nist_cfgs.loc[nist_cfgs.loc[:]['i']==icf]['CFG'].tolist()[0]
        db_X.at[i,'Term']=term
        db_X.at[i,'Conf']=cf
    # drop CF column
    db_X.drop('CF',axis=1,inplace=True)
    return db_X

def transformsymb_radtran(EX_levs,db_levs):
    EX_tranlevs=EX_levs.copy()
    # insert columns for configuration, term, multiplicity, L and parity of final and initial state
    colin_k=['Confk','Termk','Sk','Lk','Pk','Ek(Ry)']
    colin_i=['Confi','Termi','Si','Li','Pi','Ei(Ry)']
    coldb=['Conf','Term','2S+1','L','P','Energy']
    valin=['x','x',-999,-999,-999,-999.999]
    idxin=[1,2,3,4,5,7,9,10,11,12,13,15]
    ncols=len(idxin)
    ncolin=len(colin_k)
    for i in range(ncols):
        ii=idxin[i]
        if ii<8: EX_tranlevs.insert(ii,colin_k[i],valin[i])
        if ii>8: EX_tranlevs.insert(ii,colin_i[i-ncolin],valin[i-ncolin])
    # fill configuration and term symbols of final and initial states
    ntran=len(EX_tranlevs)
    for i in range(ntran):
        lvk=EX_tranlevs.loc[i]['LV_k']
        lvi=EX_tranlevs.loc[i]['LV_i']
        dummyk=db_levs.loc[db_levs[:]['LV']==lvk].squeeze()
        dummyi=db_levs.loc[db_levs[:]['LV']==lvi].squeeze()
        for j in range(ncolin):
            EX_tranlevs.at[i,colin_k[j]]=dummyk[coldb[j]]
            EX_tranlevs.at[i,colin_i[j]]=dummyi[coldb[j]]
    # remove AS's CF and LV index
    EX_tranlevs.drop(['LV_k','LV_i'],axis=1,inplace=True)
    return EX_tranlevs

def include_NIST(ttype,db_EXtran,nist_tranlevs):
    db_EXtran['source']='auto'
    nascols=['Aki','fik']
    nistcols=['Aki_NIST','fik_NIST']
    if ttype=='E1': 
        dummy_nist=nist_tranlevs.loc[nist_tranlevs.loc[:]['fik']>0]
    if ttype!='E1': 
        nascols=[nascols[0]]
        nistcols=[nistcols[0]]
        dummy_nist=nist_tranlevs.loc[nist_tranlevs.loc[:]['fik']<0]
    dummy_nist.reset_index(drop=True,inplace=True)
    nadd=len(nascols)
    ntran=len(dummy_nist)
    for i in range(nadd):
        head=nistcols[i]
        db_EXtran[head]=-9.999
    for i in range(ntran):
        dummy=dummy_nist.loc[i][:]
        idx=db_EXtran.index[(db_EXtran.loc[:]['Confk']==dummy['Confk'])&
                            (db_EXtran.loc[:]['Termk']==dummy['Termk'])&
                            (db_EXtran.loc[:]['gk']==dummy['gk'])&
                            (db_EXtran.loc[:]['Confi']==dummy['Confi'])&
                            (db_EXtran.loc[:]['Termi']==dummy['Termi'])&
                            (db_EXtran.loc[:]['gi']==dummy['gi'])].tolist()
        if len(idx)==0: nist_notfound.append(i)
        if len(idx)!=0: 
            db_EXtran.at[idx[0],'source']='nist'
            for j in range(nadd):
                head=nascols[j]
                nhead=nistcols[j]
                db_EXtran.at[idx[0],nhead]=dummy[head]
    return

def prep_print(ttype,db_EXtran):
    print_EXtran=db_EXtran.copy()
    nascols=['Aki','fik']
    nistcols=['Aki_NIST','fik_NIST']
    if ttype!='E1': 
        nascols=[nascols[0]]
        nistcols=[nistcols[0]]
    ncols=len(nascols)
    ntran=len(db_EXtran)
    for i in range(ntran):
        for j in range(ncols):
            head=nascols[j]
            nhead=nistcols[j]
            dum=print_EXtran.loc[i][nhead]
            if dum>0: print_EXtran.at[i,head]=dum
    if ttype=='E1': 
        print_EXtran['gf']=print_EXtran['gi']*print_EXtran['fik']
        nascols.append('gf')
    print_EXtran.drop(nistcols,axis=1,inplace=True)
    ncols=len(nascols)
    eformat=[]
    for i in range(ncols): eformat.append('{:.3e}')
    format_mapping=dict(zip(nascols,eformat))
    for key,value in format_mapping.items():
        print_EXtran[key]=print_EXtran[key].apply(value.format)
    return print_EXtran

#####################################################################################################
# NIST input: 
#####################################################################################################
##
## - Levels and configuration
##
mainfolder="/home/ale/mg_xtra/high_nl/"
nist_levs=pd.read_csv(mainfolder+"NIST_levels.dat",sep='\s+',skiprows=[0,2],header='infer')
nist_cfgs=pd.read_csv(mainfolder+"NIST_cfgs.dat",sep='\s+',skiprows=[0,2],header='infer')
nist_levs.rename(columns={"Level(Ry)":"NIST(Ryd)"},inplace=True)
nist_levs.rename(columns={"Configuration":"Conf"},inplace=True)
nlevnist=len(nist_levs)
# match configuration with AutoStructure labeling and decode spectroscopic terms to quantum numbers
cf=[]
sq=[]
lq=[]
pq=[]
for i in range(nlevnist):
    dumcf=nist_levs.loc[i]['Conf']
    dumterm=nist_levs.loc[i]['Term']
    sqq,lqq,pqq=termsymbol_to_quantnumber(dumterm)
    icfg=nist_cfgs.loc[nist_cfgs.loc[:]['CFG']==dumcf]['i'].tolist()
    sq.append(sqq)
    lq.append(lqq)
    pq.append(pqq)
    cf.append(icfg[0])
# insert new columns into NIST dataframe
nist_levs.insert(2,'2S+1',sq)
nist_levs.insert(3,'L',lq)
nist_levs.insert(4,'P',pq)
nist_levs.insert(5,'CF',cf)
##
## - Level radiative transition
##
skip=[i for i in range(5)]
cols=[i for i in range(12)]
header=['Wavelength(nm)','Aki','fik','Acc.','Ei(Ry)','Ek(Ry)','Confi','Termi','Ji','Confk','Termk','Jk']
nist_tranlevs=pd.read_csv(mainfolder+"NIST_lines.dat",sep='\s+',skiprows=skip,header=None,usecols=cols)
colnames=dict(zip(cols,header))
nist_tranlevs.rename(columns=colnames,inplace=True)
nist_tranlevs['gi']=nist_tranlevs['Ji']*2+1
nist_tranlevs['gk']=nist_tranlevs['Jk']*2+1
nist_tranlevs.sort_values(by=['Ei(Ry)','Ek(Ry)'],inplace=True)
nist_tranlevs.reset_index(drop=True,inplace=True)
ntran_nist=len(nist_tranlevs)
##
### >> Create NIST-terms dataframe from NIST-levels
##
nist_terms=nist_levs.drop_duplicates(subset=['Conf','Term'],keep='first')
nist_terms.reset_index(drop=True,inplace=True)
# compute weighted energy and J quantum number for each term
ntermnist=len(nist_terms)
for i in range(ntermnist):
    dumterm=nist_terms.loc[i][:]
    dumlev=nist_levs.loc[(nist_levs.loc[:]['Conf']==dumterm['Conf'])
                        &(nist_levs.loc[:]['Term']==dumterm['Term'])]
    dumlev.reset_index(drop=True,inplace=True)
    ndumlev=len(dumlev)
    sum_giei=0.
    sum_gi=0
    for j in range(ndumlev):
        gi=2*dumlev.loc[j]['J']+1
        ei=dumlev.loc[j]['NIST(Ryd)']
        sum_gi=sum_gi+gi
        sum_giei=sum_giei+gi*ei
    eiterm=sum_giei/sum_gi
    jiterm=(sum_gi-1)/2
    nist_terms.at[i,'NIST(Ryd)']=eiterm
    nist_terms.at[i,'J']=jiterm

#####################################################################################################
# AutoStructure input:
#####################################################################################################
##
## - Level data
##
ncfgmax=79
ndum=6
ncfgs=85
nlevs0=339
as_levs=pd.read_csv("oic",sep='\s+',skiprows=ncfgs+ndum,header='infer',nrows=nlevs0)
# drop all levels higher than 3s.20d
pd_droplev=as_levs.loc[as_levs.loc[:]['CF']>ncfgmax]
idx_droplev=pd_droplev.index.tolist()
as_levs.drop(idx_droplev,axis=0,inplace=True)
as_levs.reset_index(drop=True,inplace=True)
# determine parity and take absolute value of multiplicity
determine_parity(as_levs,5)
as_levs.loc[:]['2S+1']=abs(as_levs.loc[:]['2S+1'])
# rename energy column
as_levs.rename(columns={"(EK-E1)/RY":"AS(Ryd)"},inplace=True)
nlevs=len(as_levs)
##
### >> Include NIST energy levels in AutoStructure levels dataframe
##
# include new column to match NIST energy levels
iflag=-1
ncols=len(as_levs.columns)
as_levs.insert(ncols,'NIST(Ryd)',iflag)
as_levs['NIST(Ryd)']=as_levs['NIST(Ryd)'].astype(float)
# match multiplicity, angular momenta, J and configuration between NIST and AutoStructure
for i in range(nlevs):
    dumnist=nist_levs[
           (nist_levs.loc[:]['2S+1']==as_levs.loc[i]['2S+1']) & 
           (nist_levs.loc[:]['L']   ==as_levs.loc[i]['L'])    & 
           (2*nist_levs.loc[:]['J'] ==as_levs.loc[i]['2J'])   & 
           (nist_levs.loc[:]['CF']  ==as_levs.loc[i]['CF'])][:]
    if len(dumnist)==1:
        as_levs.at[i,'NIST(Ryd)']=float(dumnist.iloc[0]['NIST(Ryd)'])
# create new column 'Energy', which combines NIST and AS (shifted) energy levels
as_levs['Energy']=as_levs['NIST(Ryd)']
# copy computed (with ISHFTLS) energy levels:
imiss_levs=as_levs.index[as_levs.loc[:]['Energy']==-1].tolist()
for i in imiss_levs:
    as_levs.at[i,'Energy']=as_levs.loc[i]['AS(Ryd)']
# check if there is any missing energy level
icheck_levs=as_levs.loc[as_levs.loc[:]['Energy']==-1]
if len(icheck_levs)!=0: print("missing: ",icheck_levs)
##
## - Term data
##
nterms0=189
initskip=initskip_olg('terms')
tcols=[i for i in range(8)]
as_terms=pd.read_csv("olg",sep='\s+',skiprows=initskip,header='infer',nrows=nterms0,usecols=tcols)
# drop all terms higher than 3s.20d
pd_dropterm=as_terms.loc[as_terms.loc[:]['CF']>ncfgmax]
idx_dropterm=pd_dropterm.index.tolist()
as_terms.drop(idx_dropterm,axis=0,inplace=True)
as_terms.reset_index(drop=True,inplace=True)
# drop two useless columns
as_terms.drop(['K*CM','WEIGHTS'],axis=1,inplace=True)
# determine parity and take absolute value of multiplicity
determine_parity(as_terms,5)
as_terms.loc[:]['2S+1']=abs(as_terms.loc[:]['2S+1'])
# rename energy column
as_terms.rename(columns={"(EI-E1)/RY":"AS(Ryd)"},inplace=True)
nterms=len(as_terms)
##
### >> Include NIST energy terms in AutoStructure terms dataframe
##
# include new column to match NIST energy terms
iflag=-1
ncols=len(as_terms.columns)
as_terms.insert(ncols,'NIST(Ryd)',iflag)
as_terms['NIST(Ryd)']=as_terms['NIST(Ryd)'].astype(float)
# match multiplicity, angular momenta and configuration between NIST and AutoStructure
for i in range(nterms):
    dumnist=nist_terms[
           (nist_terms.loc[:]['2S+1']==as_terms.loc[i]['2S+1']) & 
           (nist_terms.loc[:]['L']   ==as_terms.loc[i]['L'])    & 
           (nist_terms.loc[:]['CF']  ==as_terms.loc[i]['CF'])][:]
    if len(dumnist)==1:
        as_terms.at[i,'NIST(Ryd)']=float(dumnist.iloc[0]['NIST(Ryd)'])
# copy NIST energy terms to pseudo database dataframe
as_terms['Energy']=as_terms['NIST(Ryd)']
# copy computed (with ISHFTLS) energy terms:
imiss_terms=as_terms.index[as_terms.loc[:]['Energy']==-1].tolist()
for i in imiss_terms:
    as_terms.at[i,'Energy']=as_terms.loc[i]['AS(Ryd)']
icheck_terms=as_terms.loc[as_terms.loc[:]['Energy']==-1]
if len(icheck_terms)!=0: print("missing: ",icheck_terms)
##
## - Level radiative transition 
##
### From olg file:
##
ttype='E1'
initskip=initskip_olg(ttype)
tcols=[i for i in range(4)]
ntranE1_olg=13982
E1tran=pd.read_csv("olg",sep='\s+',skiprows=initskip-1,header='infer',nrows=ntranE1_olg,usecols=tcols)
drop_doublexcited(E1tran,pd_droplev)
absolute_Aki(ttype,E1tran)
insert_LVdata(E1tran)
compute_fik(ttype,E1tran)
drop_Ncols(E1tran)

ttype='E2'
initskip=initskip_olg(ttype)
tcols=[i for i in range(5)]
ntranE2_olg=18159
E2tran=pd.read_csv("olg",sep='\s+',skiprows=initskip-1,header='infer',nrows=ntranE2_olg,usecols=tcols)
drop_doublexcited(E2tran,pd_droplev)
absolute_Aki(ttype,E2tran)
insert_LVdata(E2tran)
# compute_fik(ttype,E2tran)
drop_Ncols(E2tran)

ttype='E3'
initskip=initskip_olg(ttype)
tcols=[i for i in range(5)]
ntranE3_olg=13165
E3tran=pd.read_csv("olg",sep='\s+',skiprows=initskip-1,header='infer',nrows=ntranE3_olg,usecols=tcols)
drop_doublexcited(E3tran,pd_droplev)
absolute_Aki(ttype,E3tran)
insert_LVdata(E3tran)
# compute_fik(ttype,E3tran)
drop_Ncols(E3tran)

#####################################################################################################
# >> Create pseudo Database dataframes:
#####################################################################################################
##
## - Levels
##
db_levs=transformsymb_levs('LV',as_levs,nist_cfgs)
##
## - Terms
##
db_terms=transformsymb_levs('T',as_terms,nist_cfgs)
##
## - Level to level radiative transition 
##
db_E1tran=transformsymb_radtran(E1tran,db_levs)
db_E2tran=transformsymb_radtran(E2tran,db_levs)
db_E3tran=transformsymb_radtran(E3tran,db_levs)
##
### >> Include NIST transition data to radiative transition dataframe
##
include_NIST('E1',db_E1tran,nist_tranlevs)
include_NIST('E2',db_E2tran,nist_tranlevs)
include_NIST('E3',db_E3tran,nist_tranlevs)
##
### Check % relative error of Aki
##
check_aki=db_E1tran.copy()
drop_aki=db_E1tran.index[db_E1tran.loc[:]['Aki_NIST']<0][:]
check_aki.drop(drop_aki,inplace=True)
check_aki['erp_Aki']=(check_aki['Aki']-check_aki['Aki_NIST'])/check_aki['Aki_NIST']
avgerp=sum(check_aki['erp_Aki'])/len(check_aki)
print("average error with NIST =",avgerp)
plt.plot(check_aki['erp_Aki'],'ko')
plt.ylabel(r'$E_r\% (A_{ki})$')
plt.ylim(-10,10)
plt.savefig('erp_Aki.eps')
##
# Prepare radiative transition data to be printed
##
print_E1tran=prep_print('E1',db_E1tran)
print_E2tran=prep_print('E2',db_E2tran)
print_E3tran=prep_print('E3',db_E3tran)
##
## >> Print pseudo Database for level and term energy data
##
db_levs.to_csv('NIST+AS_levels.dat',index=False,sep='\t',header=True,float_format='%.8f')
db_terms.to_csv('NIST+AS_terms.dat',index=False,sep='\t',header=True,float_format='%.8f')
##
## >> Print pseudo Database for level and term radiative transition data
##
print_E1tran.to_csv('NIST+AS_E1tranlevels.dat',index=False,sep='\t',header=True,float_format='%.8f')
print_E2tran.to_csv('NIST+AS_E2tranlevels.dat',index=False,sep='\t',header=True,float_format='%.8f')
print_E3tran.to_csv('NIST+AS_E3tranlevels.dat',index=False,sep='\t',header=True,float_format='%.8f')

