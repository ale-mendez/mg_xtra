import os
import sys
import numpy as np
import pandas as pd
import src.module_om2ups as om2ups
import src.module_process_energies as proc_ener
import src.module_process_omg as proc_omg
import src.module_miscellaneous as misc


kB = 8.6173324E-05 # eV/K
ener = misc.default_energy_grid()
temperature = misc.default_temperature_grid()
ntemp = len(temperature)
kT = kB*np.array(temperature)



def energy_data(folder, file_om, colname, nterms, compterms = None):

    adf_path = os.path.join(folder, file_om)
    adfproc_path = os.path.join(folder, "processed_" + file_om)
    
    # process adf04 filetype
    if not os.path.isfile(adfproc_path): 
        misc.findexp_replace(adf_path, adfproc_path)

    # read AS term data
    terms = proc_ener.read_AS_TERMS(adfproc_path, nterms, colname)

    if compterms is None:
        terms.rename({colname: 'Energy'})
        idx_dict = dict(zip(terms.index, terms.index))
    else:
        if len(compterms) > 1: 
            print('error: there are more than one term dataframe to compare energies with. stop')
            sys.exit()
        # create dictionaries with index-labeling correspondance 
        #   keys: AutoStructure index 
        #   values: DataBase/NIST index
        for colname, vdict in compterms.items():
            idx_dict = misc.dict_levels_AS('T', terms, vdict['terms'], nsuplev = vdict['suplev'])
            proc_ener.insert_energy_column(terms, colname, idx_dict, vdict['terms'])

        # drop all terms without NIST+AS energy values
        terms = terms[terms[colname].notnull()]
        # create a new column with energy values that will be used throughout the ECS calculation
        terms['Energy'] = terms[colname]

    return terms, idx_dict


def collisional_data(folder, file_om, terms, nterms, ntran, debug=False):

    adfproc_path = os.path.join(folder, "processed_" + file_om)
    upsfpath = os.path.join(folder, "upsAS_unformatted.dat")

    # read computed collisions strengths (omg)
    nskip = nterms+3
    omg = proc_omg.read_AS_omg(adfproc_path, ener, nskip, ntran)

    # sort transition type
    omg = proc_omg.sort_transition_type(omg, terms)
    nomg1 = len(omg)

    # drop all transitions with terms not included in database
    omg = omg[omg['type'] > 0]
    omg = omg.reset_index(drop=True)
    nomg2 = len(omg)
    print(f' ==> {nomg1-nomg2} transition were dropped because eik==0\n')

    # read or compute ECS 
    if os.path.exists(upsfpath):
        ups = pd.read_csv(upsfpath)
    else:
        ups = om2ups.compute_ECS(terms, omg, ener, temperature, debug=debug)
        ups.to_csv(upsfpath, index=False)

    return omg, ups


def rename_transition_levels(df, levdict):

    print('processing index transformation..')

    r_df = df.copy()
    kp, ip = [], []
    for tran in r_df.index:
        k, i = r_df.loc[tran]['k'], r_df.loc[tran]['i']
        if k != levdict[k]: k = levdict[k]
        if i != levdict[i]: i = levdict[i]
        kp.append(int(k))
        ip.append(int(i))
    r_df.insert(loc=2, column='kp', value=kp)
    r_df.insert(loc=3, column='ip', value=ip)
    print("... all done!")

    return r_df


def compute_superlevels_ECS(terms, ups, nlevmax, nlevmin, debug=False):

    ups_suplev = om2ups.initialize_upsilon_dataframe(temperature, nlevmax=nlevmax, nlevmin=nlevmin)

    iempty = []
    for tran in ups_suplev.index:
        gL, gU, eL, eU = 0, 0, 0, 0
        k_slv, i_slv = ups_suplev.loc[tran, 'k'], ups_suplev.loc[tran, 'i']
        dummy = ups.loc[(ups['kp'] == k_slv) & (ups['ip'] == i_slv)]
        if len(dummy) > 1: # compute superlevels
            if debug:
                print(k_slv, "=>", i_slv, "//", dummy['k'].tolist(), "=>", dummy['i'].tolist())
            tot_ups = np.zeros(ntemp)
            for slv_tran in dummy.index:
                k, i = dummy.loc[slv_tran]['k'], dummy.loc[slv_tran]['i']
                if k < i: i, k = k, i
                gk, gi = 2 * terms.loc[k, 'J'] + 1, 2 * terms.loc[i, 'J'] + 1
                ek, ei = terms.loc[k, 'Energy'], terms.loc[i, 'Energy']
                eik = ei - ek
                gL += gi
                eL += gi * ei
                gU += gk
                eU += gk * ek
                tran_ups = dummy.loc[slv_tran][5:].values
                tot_ups += np.exp(-eik/kT) * tran_ups/gi
            eLU = eL/gL - eU/gU
            ups_suplev.iloc[tran, 2:] = gL * np.exp(eLU/kT) * tot_ups
        elif len(dummy) == 1: # no superlevels
            ups_suplev.iloc[tran, 2:] = dummy.iloc[0][5:].values
        elif len(dummy) == 0: # transitions not found
            iempty.append((k_slv, i_slv))
    iempty_rows = ups_suplev.loc[ups_suplev[3000.0] == 0].index
    if len(iempty) != len(iempty_rows): print('empty rows inconsistent!!')
    ups_suplev.drop(iempty_rows, inplace=True)

    return ups_suplev


def compute_superlevels_aki(terms, ups, nlevmax, nlevmin, debug=False):

    ups_suplev = om2ups.initialize_upsilon_dataframe(temperature, nlevmax=nlevmax, nlevmin=nlevmin)

    iempty = []
    for tran in ups_suplev.index:
        gL, gU, eL, eU = 0, 0, 0, 0
        k_slv, i_slv = ups_suplev.loc[tran, 'k'], ups_suplev.loc[tran, 'i']
        dummy = ups.loc[(ups['kp'] == k_slv) & (ups['ip'] == i_slv)]
    #     if len(dummy) > 1: # compute superlevels
    #         if debug:
    #             print(k_slv, "=>", i_slv, "//", dummy['k'].tolist(), "=>", dummy['i'].tolist())
    #         tot_ups = np.zeros(ntemp)
    #         for slv_tran in dummy.index:
    #             k, i = dummy.loc[slv_tran]['k'], dummy.loc[slv_tran]['i']
    #             if k < i: i, k = k, i
    #             gk, gi = 2*terms.loc[k, 'J'] + 1, 2*terms.loc[i, 'J'] + 1
    #             ek, ei = terms.loc[k, 'Energy'], terms.loc[i, 'Energy']
    #             eik = ei - ek
    #             gL += gi
    #             eL += gi*ei
    #             gU += gk
    #             eU += gk*ek
    #             tran_ups = dummy.loc[slv_tran][5:].values
    #             tot_ups += np.exp(-eik/kT) * tran_ups/gi
    #         eLU = eL/gL - eU/gU
    #         ups_suplev.iloc[tran, 2:] = gL*np.exp(eLU/kT)*tot_ups
    #     elif len(dummy) == 1: # no superlevels
    #         ups_suplev.iloc[tran, 2:] = dummy.iloc[0][5:].values
    #     elif len(dummy) == 0: # transitions not found
    #         iempty.append((k_slv, i_slv))
    # iempty_rows = ups_suplev.loc[ups_suplev[3000.0] == 0].index
    # if len(iempty) != len(iempty_rows): print('empty rows inconsistent!!')
    # ups_suplev.drop(iempty_rows, inplace=True)
    # return ups_suplev


def process_superlevels(terms, ups, termdic_DB, nlevmax, nlevmin, debug=False):

    print('processing index transformation... and',
          ' \ncomputing ECS of superlevels at the same time...')

    # create effective collision strength dataframe
    ups_suplev1 = om2ups.initialize_upsilon_dataframe(temperature, nlevmax, nlevmin)
    ups_suplev2 = om2ups.initialize_upsilon_dataframe(temperature, nlevmax, nlevmin)

    iempty = []
    for tran in ups_suplev1.index:
        gL, eL = 0, 0
        klv, ilv = ups_suplev1.loc[tran, 'k'], ups_suplev1.loc[tran, 'i']
        i_as = [idx for idx, kdb in termdic_DB.items() if kdb == ilv]
        k_as = [idx for idx, kdb in termdic_DB.items() if kdb == klv]
        if debug: 
            print(f'tran: {tran:>6d}, {ilv:<3d} => {klv:>3d}  //  {i_as} => {k_as}')
        tot_ups1 = np.zeros(ntemp)
        tot_ups2 = np.zeros(ntemp)
        for i in i_as:
            gi, ei = 2*terms.loc[i, 'J'] + 1, terms.loc[i, 'Energy']
            gL += gi
            eL += gi*ei
            for k in k_as:
                ek = terms.loc[k, 'Energy']
                eik = ei - ek
                if k < i: i, k = k, i
                if debug: print(15*' ', f'{i:>3d} => {k:>3d}')
                if eik == 0: 
                    iempty.append((tran, i, k))
                umatch = ups[(ups['k'] == k) & (ups['i'] == i)]
                if umatch.empty: continue
                tran_ups = umatch.iloc[0][3:].values
                # statistic average (see notes)
                tot_ups1 += np.exp(-eik/kT) * tran_ups/gi
                tot_ups2 += gi * tran_ups
        gU = sum([2*terms.loc[k, 'J'] + 1 for k in k_as])
        eU = sum([(2*terms.loc[k]['J']+1)*terms.loc[k]['Energy'] for k in k_as])
        eLU = eL/gL - eU/gU
        # copy results in df
        ups_suplev1.iloc[tran,2:] = gL*np.exp(eLU/kT)*tot_ups1
        ups_suplev2.iloc[tran,2:] = tot_ups2

    print('\nlooking for empty rows...')
    # finally, check if there is any transition in the superlevel ECS 
    # dataframe that is empty
    ierr = False
    izero = 0
    for tran in ups_suplev1.index:
        if ups_suplev1.loc[tran][2] == 0: 
            ifound = 0
            izero += 1
            for tupla in iempty:
                if tran == tupla[0]: 
                    ifound = 1
                    break
            if ifound == 0:
                ierr = True
                print(f'\touch! {tran} is empty and eik != 0.') 
    print('==> number of transition not found', izero)

    if ierr: 
        print('=> hold on! empty rows may be due to transitions',
              'not included in current calculation. check it!') 

    print("\nyou're all done! :)")
    return ups_suplev1, ups_suplev2
