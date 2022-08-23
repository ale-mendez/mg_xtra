import pandas as pd
import os
import src.module_miscellaneous as misc

kB=8.6173324E-05 # eV/K
convRyd2eV=13.6057 # eV/Ryd
convcm2Ryd=1./109737.26 


# def initskip_olg():
#     os.system("grep -n -i 'LIST OF TERMS WITH A WEIGHTED MEAN OVER THE FINE STRUCTURE' olg > dum.txt")
#     dumline = pd.read_csv("dum.txt", header=None, sep='\s+')
#     nrows = len(dumline)
#     iline = dumline.loc[nrows - 1][0]
#     os.system("rm dum.txt")
#     return int(iline[0:-1])


def read_AS_TERMS(termpath, nterms, colname):
    # input level data
    terms = pd.read_csv(termpath, header=None, skiprows=[0], nrows=nterms, usecols=[0,1,2,3,4])
    colsname = ['CFG', 'S', 'L', 'J', 'E']
    terms.columns = colsname
    for j in range(nterms):
        s = terms.loc[j]['CFG']
        sn = " ".join(s.split()[1:])
        terms.at[j,'CFG'] = sn
    terms.index = terms.index + 1
    print('checking energy conversion:')
    print(terms.loc[2]['CFG'], terms.loc[2]['E'], 'cm^2 => ', terms.loc[2]['E']*convcm2Ryd, 'Ryd\n')
    terms['E'] = terms['E'] * convcm2Ryd
    terms = terms.rename(columns={'E': colname, 'S': '2S+1'})
    terms.insert(0, 'idx', terms.index.copy())
    return terms


def read_NIST_levels(levpath, cfgpath):
    nist_levs = pd.read_csv(levpath, sep='\s+', skiprows=[0,2], header='infer')
    nist_cfgs = pd.read_csv(cfgpath, sep='\s+', skiprows=[0,2], header='infer')
    nist_levs.rename(columns={"Level(Ry)": "NIST(Ryd)"}, inplace=True)
    nist_levs.rename(columns={"Configuration": "Conf"}, inplace=True)
    nlevnist = len(nist_levs)
    # match configuration with AutoStructure labeling and decode spectroscopic terms to quantum numbers
    cf = []
    sq = []
    lq = []
    pq = []
    for i in range(nlevnist):
        dumcf = nist_levs.loc[i]['Conf']
        dumterm = nist_levs.loc[i]['Term']
        sqq, lqq, pqq = misc.termsymbol_to_quantnumber(dumterm)
        icfg = nist_cfgs.loc[nist_cfgs.loc[:]['CFG'] == dumcf]['i'].tolist()
        sq.append(sqq)
        lq.append(lqq)
        pq.append(pqq)
        cf.append(icfg[0])
    # insert new columns into NIST dataframe
    nist_levs.insert(2, 'S', sq)
    nist_levs.insert(3, 'L', lq)
    nist_levs.insert(4, 'P', pq)
    nist_levs.insert(5, 'CFG', cf)
    return nist_levs


def determine_NIST_terms(nist_levs):
    nist_terms = nist_levs.drop_duplicates(subset=['Conf','Term'], keep='first')
    nist_terms.reset_index(drop=True, inplace=True)
    # compute weighted energy and J quantum number for each term
    ntermnist = len(nist_terms)
    for i in range(ntermnist):
        dumterm = nist_terms.loc[i][:]
        dumlev = nist_levs.loc[
            (nist_levs.loc[:]['Conf'] == dumterm['Conf']) & 
            (nist_levs.loc[:]['Term'] == dumterm['Term'])
        ]
        dumlev.reset_index(drop=True, inplace=True)
        ndumlev = len(dumlev)
        sum_giei = 0.
        sum_gi = 0
        for j in range(ndumlev):
            gi = 2 * dumlev.loc[j]['J'] + 1
            ei = dumlev.loc[j]['NIST(Ryd)']
            sum_gi = sum_gi + gi
            sum_giei = sum_giei + gi * ei
        eiterm = sum_giei / sum_gi
        jiterm = (sum_gi - 1) / 2
        nist_terms.at[i,'NIST(Ryd)'] = eiterm
        nist_terms.at[i,'J'] = jiterm
    nist_terms.insert(7, 'gi', [2 * j + 1 for j in nist_terms.J])
    nist_terms.index = nist_terms.index + 1
    # nist_terms = nist_terms.rename(columns={'Energy': 'Energy(Ryd)'})
    return nist_terms


def read_database_energies(dbenerpath):
    mod_levels = pd.read_csv(dbenerpath, sep="\s+", header='infer', usecols=[3,6,7,8,9,11,12])
    mod_levels.rename(columns={'2S': 'S'}, inplace=True)
    mod_levels['S'] = mod_levels['S'] + 1
    mod_levels['Energy(Ryd)'] = mod_levels['ExcitationWaven'] * convcm2Ryd
    mod_levels = mod_levels.drop(['LevelNumber', 'ExcitationWaven'], axis=1)
    mod_levels = mod_levels.rename({'ElectronConfig': "Conf", 'LevelWeight': "gi"}, axis=1)
    mod_levels.index = mod_levels.index + 1
    return mod_levels


def read_NIST_terms(termspath):
    terms = pd.read_csv(termspath, sep='\t', usecols=[1,2,3,4,5,6,9])
    terms = terms.rename(columns={'2S+1': 'S', 'CF': 'CFG'})
    terms.insert(6, 'gi', [s*(2*l+1) for s, l in zip(terms['S'], terms['L'])])
    terms.index = terms.index + 1
    terms = terms.rename(columns={'Energy': 'Energy(Ryd)'})
    return terms


def insert_energy_column(terms, colname, termdic, df):
    terms[colname] = None
    for idx_as, idx in termdic.items():
        terms.at[idx_as, colname] = df.loc[idx]['Energy(Ryd)']


def combine_energy_values(terms, colnameA, colnameB):
    combvalues = []
    for i in terms.index:
        ei = terms.loc[i][colnameA]
        if ei is None: ei = terms.loc[i][colnameB]
        combvalues.append(ei)
    return combvalues


def read_levels_from_oic(oicpath, ncfgs, nlevs, colname, ncfgmax=None):
    """Read 

    Parameters
    ----------
    oicpath : _type_
        _description_
    ncfgs : _type_
        _description_
    nlevs : _type_
        _description_
    colname : _type_
        _description_
    ncfgmax : _type_, optional
        _description_, by default None

    Returns
    -------
    _type_
        _description_
    """
    ndum = 6
    as_levs = pd.read_csv(oicpath, sep='\s+', skiprows=ncfgs + ndum, header='infer', nrows=nlevs)
    # drop all levels higher than ncfgmax
    if ncfgmax is not None: 
        pd_droplev = as_levs.loc[as_levs.loc[:]['CF'] > ncfgmax]
        idx_droplev = pd_droplev.index.tolist()
        as_levs.drop(idx_droplev, axis=0, inplace=True)
        as_levs.reset_index(drop=True, inplace=True)
    # determine parity and take absolute value of multiplicity
    misc.determine_parity(as_levs, 5)
    as_levs['2S+1'] = abs(as_levs['2S+1'])
    as_levs['2J'] = as_levs['2J'] / 2
    # include string format configurations
    cfgs = [misc.config(i) for i in as_levs['CF']]
    as_levs.insert(2, 'CFG', cfgs)
    terms = [misc.quantnumber_to_termsymbol(sq, lq, pq) for sq, lq, pq in zip(as_levs['2S+1'], as_levs['L'], as_levs['P'])]
    as_levs.insert(3, 'Term', terms)
    as_levs.drop(axis=1, columns=['CF'], inplace=True)
    # rename columns and index
    as_levs.rename(columns={"(EK-E1)/RY": colname, "2J": "J"}, inplace=True)
    as_levs.index = as_levs.index + 1
    return as_levs


def read_terms_from_olg(olgpath, nterms, colname, ncfgmax=None):
    """Read term energy values from olg file

    Parameters
    ----------
    olgpath : str
        Path to olg file
        
    nterms : int
        Number of terms to be read
        
    colname : str
        Name for energy column to be given to calculation
        
    ncfgmax : int, optional
        Number of maximum configuration value, by default None 
        Values with configurations greater than ncfgmax will be dropped

    Returns
    -------
    Dataframe
        Dataframe with energy calculation values
        
    """
    initskip = misc.initskip_olg(olgpath, 'terms', 'terms')
    tcols = [i for i in range(8)]
    as_terms = pd.read_csv(olgpath, sep='\s+', skiprows=initskip, header='infer', nrows=nterms, usecols=tcols)
    # drop all terms higher than 3s.20d
    if ncfgmax is not None: 
        pd_dropterm = as_terms.loc[as_terms.loc[:]['CF'] > ncfgmax]
        idx_dropterm = pd_dropterm.index.tolist()
        as_terms.drop(idx_dropterm, axis=0, inplace=True)
        as_terms.reset_index(drop=True, inplace=True)
    # drop two useless columns
    as_terms.drop(['K*CM', 'WEIGHTS'], axis=1, inplace=True)
    # determine parity and take absolute value of multiplicity
    misc.determine_parity(as_terms, 5)
    as_terms['2S+1'] = abs(as_terms['2S+1'])
    gi = as_terms['2S+1'] * (as_terms['L'] * 2 + 1)
    as_terms.insert(4, "J", (gi - 1) / 2)
    # include string format configurations
    cfgs = [misc.config(i) for i in as_terms['CF']]
    as_terms.insert(2, 'CFG', cfgs)
    as_terms.drop(axis=1, columns=['CF'], inplace=True)
    terms = [misc.quantnumber_to_termsymbol(sq, lq, pq) for sq, lq, pq in zip(as_terms['2S+1'], as_terms['L'], as_terms['P'])]
    as_terms.insert(3, 'Term', terms)
    # rename columns and index
    # as_terms.rename(columns={'I': 'K', "(EI-E1)/RY": colname, '2S+1': 'S'}, inplace=True)
    as_terms.rename(columns={'I': 'K', "(EI-E1)/RY": colname}, inplace=True)
    as_terms.index = as_terms.index + 1
    return as_terms
