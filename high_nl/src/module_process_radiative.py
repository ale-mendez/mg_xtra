'''
Radiative transition conversion:

For allowed transitions, we transform $A_{ki}$ to $f_{ik}$ using the relation:

A_{ki} = 2\pi e^2/(m_e c \epsilon_0) (1/\lambda^2) (g_i/g_k) (f_{ik})

=> f_{ik} = m_e c \epsilon_0/(2\pi e^2) \lambda^2 (g_k/g_i) A_{ki}

'''
import os
import pandas as pd
import numpy as np
import src.module_miscellaneous as misc


def read_NIST_radiative(folder):
    skip = [i for i in range(5)]
    cols = [i for i in range(12)]
    header = ['Wavelength(nm)', 'Aki', 'fik' ,'Acc.', 'Ei(Ry)', 'Ek(Ry)', 'Confi', 'Termi', 'Ji', 'Confk', 'Termk', 'Jk']
    nistpath = os.path.join(folder, "NIST_lines.dat")
    nist_tranlevs = pd.read_csv(nistpath, sep='\s+', skiprows=skip, header=None, usecols=cols)
    colnames = dict(zip(cols, header))
    nist_tranlevs.rename(columns=colnames, inplace=True)
    nist_tranlevs['gi'] = nist_tranlevs['Ji'] * 2 + 1
    nist_tranlevs['gk'] = nist_tranlevs['Jk'] * 2 + 1
    nist_tranlevs.sort_values(by=['Ei(Ry)', 'Ek(Ry)'], inplace=True)
    nist_tranlevs.reset_index(drop=True, inplace=True)
    return nist_tranlevs


def read_transition_from_olg(keyword, olgpath, trantype, ntran, as_ener):
    """Read transition data from olg file

    Parameters
    ----------
    keyword : str
        Keyword to tell the code to read term-to-term or level-to-level transition calculations. Accepted values: 'terms', 'levels'
        
    olgpath : str
        Path file to olg datafile
        
    type : str
        Transition type. Accepted values: 'E1', 'E2', 'E3'.
        
    ntran : int
        Number of transitions to be read
        
    as_ener : Dataframe
        Dataframe with autostructure energy calculations

    Returns
    -------
    Dataframe
        Dataframe with transition values

    Raises
    ------
    ValueError
        if keyword/trantype is incorrect, it prints out the only accepted values.
        
    """
    if keyword not in ['terms', 'levels']:
        raise ValueError(f"{keyword} can only be 'terms' or 'levels'")
    if trantype not in ['E1', 'E2', 'E3']:
        raise ValueError(f"{trantype} can only be 'E1', 'E2' or 'E3'")
    
    Ecolname = 'database(NIST)'
    initskip = misc.initskip_olg(olgpath, trantype, keyword)
    ncols = 5
    if trantype == 'E1':
        ncols = 4
    tcols = [i for i in range(ncols)]
    tran = pd.read_csv(olgpath, sep='\s+', skiprows=initskip - 1, header='infer', nrows=ntran, usecols=tcols)
    if keyword == 'terms':
        tran.rename(columns={'I': 'K', 'IP': 'KP'}, inplace=True)
        kw = 'T'
        trantype = 'E1' # E2 and E3 transitions (in LS) are printed out in olg as E1 transition
    elif keyword == 'levels':
        kw = 'LV'
    absolute_Aki(trantype, tran)
    insert_energy_from_dataframe(kw, tran, as_ener, Ecolname)
    if trantype == 'E1':
        compute_fik(trantype, tran)
    # drop NaN valued rows and unnecesary columns
    tran = tran[tran['Eki'].notnull()]
    drop_Ncols(tran)
    return tran


def conv_Aki2fik(ttype, gi, gk, Eki, Aki):
    #
    # convert Aki to fik using SI units
    #
    #    Eki : rydbergs
    #    Aki : s^{-1} 
    #
    if Eki is None:
        return None
    IP = 13.6056923           # eV/ry
    h = 4.135667731E-15       # eV*s
    h_J = 6.62607E-34         # J*s
    c = 299792458             # m/s
    hc = h * c                  # eV*m
    me = 9.10938E-11          # kg
    e = 1.60218E-19           # C
    eps0 = 8.8541878128E-12   # F/m and F= C^2 J^-1 m^-1 
    convJ2eV = 6.24E+18       # eV/J
    convm2A = 1.0E+10         # Angstrom
    convm2nm = 1.0E+9         # nm
    convEh2J = 4.35974417E-18 # J/a.u. (hartree)
    convry2J = convEh2J / 2     # J/ry
    lam = hc / (Eki * IP)         # eV*m/eV = m
    lamAng = lam * convm2A
    if ttype == 'E1':
        E1_f2A = 2 * np.pi * e ** 2 / (me * c * eps0) * 1.0E20
        fik = lam ** 2 * gk * Aki / (E1_f2A * gi)
    if ttype != 'E1':
        fik = None
    return fik


def drop_doublexcited(EXtran, pd_droplev):
    #
    # drop double excited levels except 3p^2
    #
    levsk = pd_droplev.loc[:]['K'].tolist()
    for i in levsk:
        drop_k = EXtran.index[EXtran.loc[:]['K'] == i].tolist()
        EXtran.drop(drop_k, axis=0, inplace=True)
    EXtran.reset_index(drop=True, inplace=True)


def absolute_Aki(ttype, EXtran):
    #
    # take absolute value of electric and magnetic transition elements
    # sum electric and magnetic einstein coefficients
    #
    
    if ttype == 'E1': 
        Aki_cols = ['A(EK)*SEC']
    else:
        Aki_cols = ['A(EK)*SEC', 'A(MK)*SEC']
        
    for i in Aki_cols: 
        EXtran[i] = abs(EXtran[i])
        
    EXtran['Aki'] = sum(EXtran[i] for i in Aki_cols)


def insert_energy_from_dataframe(key, EXtran, as_data, Ecolname):
    #
    # use k and kp index (from as_data dataframe) to insert lv index
    # include statistic weight 
    # replace Eki with observed values (given in Ecolname column)
    #
    colin = [key + '_k', 'gk', key + '_i', 'gi', 'Eki']
    valin = [-999, -999, -999, -999, -999.999]
    idxin = [3, 4, 5, 6, 7]
    for i in range(5): 
        EXtran.insert(idxin[i], colin[i], valin[i])
    ntran = len(EXtran)
    for ii in range(ntran):
        dummy_k = as_data.loc[as_data['K'] == EXtran.loc[ii]['K']]
        dummy_i = as_data.loc[as_data['K'] == EXtran.loc[ii]['KP']]
        # si encuentra ese level/term, calculo el salto de energia
        if (len(dummy_i) > 0) and (len(dummy_k) > 0):
            Ek = dummy_k[Ecolname].tolist()[0]
            Ei = dummy_i[Ecolname].tolist()[0]
            EXtran.at[ii, key + '_k'] = dummy_k[key].tolist()[0]
            EXtran.at[ii, key + '_i'] = dummy_i[key].tolist()[0]
            if key == 'LV':
                EXtran.at[ii,'gk'] = 2 * dummy_k['J'].tolist()[0] + 1
                EXtran.at[ii,'gi'] = 2 * dummy_i['J'].tolist()[0] + 1
            if key == 'T':
                EXtran.at[ii, 'gk'] = dummy_k['2S+1'].tolist()[0] * (2 * dummy_k['L'].tolist()[0] + 1)
                EXtran.at[ii, 'gi'] = dummy_i['2S+1'].tolist()[0] * (2 * dummy_i['L'].tolist()[0] + 1)
            EXtran.at[ii, 'Eki'] = Ek - Ei
        else: # sino flag con None
            EXtran.at[ii, 'Eki'] = None


def compute_fik(ttype, EXtran):
    #
    # compute absortion oscillator strength and weigthed fik for allowed transitions only
    #
    EXtran['fik'] = conv_Aki2fik(ttype, EXtran['gi'], EXtran['gk'], EXtran['Eki'], EXtran['Aki'])
    if ttype == 'E1': EXtran['gf'] = EXtran['gi'] * EXtran['fik']
    if ttype != 'E1': EXtran['gf'] = None


def drop_Ncols(EXtran):
    #
    # drop unnecesary columns
    #
    dropcols = []
    cols = ['E1-DATA', 'E2/M1-DATA', 'E3/M2-DATA', 'E2-DATA', 'E3-DATA', 'K', 'KP', 'A(EK)*SEC', 'A(MK)*SEC', 'Eki', 'S']
    for i in cols:
        if i in EXtran.columns: dropcols.append(i)
    EXtran.drop(dropcols, axis=1, inplace=True)


def transformsymb_levs(ttype, as_X, nist_cfgs):
    db_X = as_X.copy()
    # move ttype column to front
    X_list = db_X.loc[:][ttype].tolist()
    coldrop = [ttype]
    if ttype == 'LV': 
        coldrop.append('T')
    db_X.drop(coldrop, axis=1, inplace=True)
    db_X.insert(0, ttype, X_list)
    # insert columns for configuration and term symbols
    db_X.insert(2, 'Conf', 'x')
    db_X.insert(3, 'Term', 'x')
    # fill configuration and term columns with symbols
    ndata = len(db_X)
    for i in range(ndata):
        sq = db_X.iloc[i]['2S+1']
        lq = db_X.iloc[i]['L']
        pq = db_X.iloc[i]['P']
        icf = db_X.iloc[i]['CFG']
        term = misc.quantnumber_to_termsymbol(sq, lq, pq)
        cf = nist_cfgs.loc[nist_cfgs.loc[:]['i'] == icf]['CFG'].tolist()[0]
        db_X.at[i, 'Term'] = term
        db_X.at[i, 'Conf'] = cf
    # drop CF column
    db_X.drop('CFG', axis=1, inplace=True)
    return db_X


def transformsymb_radtran(key, EX_levs, db_levs):
    EX_tranlevs = EX_levs.copy()
    # insert columns for configuration, term, multiplicity, L and parity of final and initial state
    colin_k = ['k', 'Confk', 'Termk', 'Sk', 'Lk', 'Pk', 'Ek(Ry)']
    colin_i = ['i', 'Confi', 'Termi', 'Si', 'Li', 'Pi', 'Ei(Ry)']
    coldb = ['K', 'CFG', 'Term', '2S+1', 'L', 'P', 'database(NIST)']
    valin = [None, 'x', 'x', -999, -999, -999, -999.999]
    idxin = [1, 2, 3, 4, 5, 6, 8, 10, 11, 12, 13, 14, 15, 17]
    ncols = len(idxin)
    ncolin = len(colin_k)
    for i in range(ncols):
        ii = idxin[i]
        if ii < 9:
            EX_tranlevs.insert(ii, colin_k[i], valin[i])
        elif ii > 9:
            EX_tranlevs.insert(ii, colin_i[i - ncolin], valin[i - ncolin])
    # fill configuration and term symbols of final and initial states
    for i in EX_tranlevs.index:
        lvk = EX_tranlevs.loc[i][key + '_k']
        lvi = EX_tranlevs.loc[i][key + '_i']
        dummyk = db_levs.loc[db_levs[:][key] == lvk].squeeze()
        dummyi = db_levs.loc[db_levs[:][key] == lvi].squeeze()
        for j in range(ncolin):
            EX_tranlevs.at[i, colin_k[j]] = dummyk[coldb[j]]
            EX_tranlevs.at[i, colin_i[j]] = dummyi[coldb[j]]
    # remove AS's CF and LV index
    EX_tranlevs.drop([key + '_k', key + '_i'], axis=1, inplace=True)
    return EX_tranlevs


def include_NIST(ttype, db_tran, nist_tranlevs):
    db_EXtran = db_tran.copy()
    # nascols = ['Aki', 'fik']
    # nistcols = ['Aki(NIST)', 'fik(NIST)']
    nascols = ['Aki']
    nistcols = ['Aki(NIST)']
    # if ttype == 'E1':
    #     dummy_nist = nist_tranlevs.loc[nist_tranlevs.loc[:]['fik'] > 0]
    # if ttype != 'E1': 
    #     nascols = [nascols[0]]
    #     nistcols = [nistcols[0]]
    #     dummy_nist = nist_tranlevs.loc[nist_tranlevs.loc[:]['fik'] < 0]
    nist_tranlevs.reset_index(drop=True, inplace=True)
    ntran = len(nist_tranlevs)
    for col in nistcols:
        db_EXtran[col] = None
        
    db_EXtran['source'] = 'auto'

    for i in range(ntran):
        dummy = nist_tranlevs.loc[i][:]
        idx = db_EXtran.index[
            (db_EXtran.loc[:]['Confk']==dummy['Confk']) &
            (db_EXtran.loc[:]['Termk']==dummy['Termk']) &
            (db_EXtran.loc[:]['gk']==dummy['gk']) &
            (db_EXtran.loc[:]['Confi']==dummy['Confi']) &
            (db_EXtran.loc[:]['Termi']==dummy['Termi']) &
            (db_EXtran.loc[:]['gi']==dummy['gi'])
        ].tolist()
        if len(idx) != 0: 
            db_EXtran.at[idx[0],'source'] = 'nist'
            for col, nist_col in zip(nascols, nistcols):
                db_EXtran.at[idx[0], nist_col] = dummy[col]
    return db_EXtran


def prep_print(ttype, db_EXtran):
    print_EXtran = db_EXtran.copy()
    nascols = ['Aki', 'fik']
    nistcols = ['Aki_NIST', 'fik_NIST']
    if ttype != 'E1': 
        nascols = [nascols[0]]
        nistcols = [nistcols[0]]
    ncols = len(nascols)
    ntran = len(db_EXtran)
    for i in range(ntran):
        for j in range(ncols):
            head = nascols[j]
            nhead = nistcols[j]
            dum = print_EXtran.loc[i][nhead]
            if dum is not None:
                print_EXtran.at[i, head] = dum
    if ttype == 'E1': 
        print_EXtran['gf'] = print_EXtran['gi'] * print_EXtran['fik']
        nascols.append('gf')
    print_EXtran.drop(nistcols, axis=1, inplace=True)
    ncols = len(nascols)
    eformat = []
    for i in range(ncols): 
        eformat.append('{:.3e}')
    format_mapping = dict(zip(nascols,eformat))
    for key, value in format_mapping.items():
        print_EXtran[key] = print_EXtran[key].apply(value.format)
    return print_EXtran


def prep_print_tranterms(ttype, db_EXtran_terms):
    print_EXtran_terms = db_EXtran_terms.copy()
    nascols = [['Aki']]
    if ttype == 'E1': nascols.append(['fik', 'gf'])
    nascols = [item for i in nascols for item in i ]
    ncols = len(nascols)
    eformat = []
    for i in range(ncols): eformat.append('{:.3e}')
    format_mapping = dict(zip(nascols, eformat))
    for key, value in format_mapping.items():
        print_EXtran_terms[key] = print_EXtran_terms[key].apply(value.format)
    return print_EXtran_terms
