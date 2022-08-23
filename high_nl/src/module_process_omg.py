import pandas as pd

def read_AS_omg(omgpath, ener, nskip, ntran):
    nener = len(ener)
    # input collision strength data
    omgcols = [i for i in range(nener + 4)]
    omg = pd.read_csv(omgpath, sep="\s+", header=None,skiprows=nskip, nrows=ntran, usecols=omgcols)
    colsname = ['k', 'i', 'aki']
    for i in range(nener):
        colsname.append(str(ener[i]))
    colsname.append("inf")
    omg.columns = colsname
    return omg

def sort_transition_type(omg, terms, debug=False):
    print('sorting transition type...')
    ntran = len(omg)
    columns = terms.columns
    fbig = 0.01
    fzero = 1.0E-04
    ntype = []
    for tran in range(ntran):
        row = omg.loc[tran]
        k, i, Aki = int(row['k']), int(row['i']), row['aki']
        # check if transition terms are found in energy dataframe
        try:
            jk, ji = terms.loc[k]['J'], terms.loc[i]['J']
        except:
            ntype.append(0)
            continue
        gk, gi = 2*jk+1, 2*ji+1
        # read term energy from database (NIST), if not found consider NIST+AS
        ek, ei = terms.loc[k]['Energy'], terms.loc[i]['Energy']
        # compute oscillator strength
        eik = abs(ek-ei)
        if eik == 0: 
            if debug: print(f'tran: {tran:>6d}, {i:>3d} => {k:>3d} | ei, ek = {ei}')
            ntype.append(0)
            continue
        S = 3.73491E-10*gk*Aki/eik**3
        fij = eik*S/(3.*gi)
        # check delta S
        if terms.loc[k]['2S+1'] == terms.loc[i]['2S+1']:
            # check delta L
            if (abs(terms.loc[k]['L']-terms.loc[i]['L'])<=1) & (fij>=fbig):
                ntype.append(1)
            else:
                # check oscillator strength value
                if (fij>fzero) & (fij<fbig):
                    ntype.append(4)
                else:
                    ntype.append(2)
        else:
            if (fij>fzero) & (fij<fbig):
                ntype.append(4)
            else:
                ntype.append(3)
    omg['type'] = ntype
    print('... all done!')
    return omg