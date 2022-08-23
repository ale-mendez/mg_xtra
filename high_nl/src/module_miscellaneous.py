import os
import sys
# from src.module_configurations import cfgs
import json
import pandas as pd


cfgfolder = '/home/ale/magnesio/mg_xtra/high_nl/src'
with open(cfgfolder + '/configurations.json') as json_file:
    data = json.load(json_file)
    cfgs_AS = data["autostructure"]
    cfgs_DB = data["database"]

cfgs = cfgs_DB.copy()
cfgs.update(cfgs_AS)
df_config = pd.DataFrame(data=cfgs_DB.values(), index=cfgs_DB.keys())
df_config = df_config.drop_duplicates(keep='last')
df_config = df_config.reset_index()
df_config.rename(columns={'index': 'CFG', 0: 'i'}, inplace=True)


def default_energy_grid():
    ener = [1.00E-02,2.20E-01,4.89E-01,8.17E-01,1.22E+00,
            1.71E+00,2.30E+00,3.03E+00,3.92E+00,5.00E+00]
    return ener


def default_temperature_grid():
    temp = [3.0E+03, 4.0E+03, 5.0E+03, 6.0E+03,7.0E+03, 8.0E+03, 9.0E+03,
            1.0E+04, 1.1E+04, 1.2E+04, 1.3E+04,1.4E+04, 1.5E+04, 1.6E+04]
    return temp


def findexp_replace(filepath, procpath):
    print("converting " + filepath)
    try:
        os.system("cp " + filepath + " " + procpath)
        os.system("sed -i 's/ /  /g' " + procpath)
        os.system("sed -i 's/+/E+/g' " + procpath)
        os.system("sed -i 's/-/E-/g' " + procpath)
        os.system("sed -i 's/(/,/g' " + procpath)
        os.system("sed -i 's/)/,/g' " + procpath)
        for i in range(99):
            if i<10: os.system("sed -i 's/+0" + str(i) + "E/+0" + str(i) + " /g' " + procpath)
            if i<10: os.system("sed -i 's/-0" + str(i) + "E/-0" + str(i) + " /g' " + procpath)
            if i>=10: os.system("sed -i 's/+" + str(i) + "E/+" + str(i) + " /g' " + procpath)
            if i>=10: os.system("sed -i 's/-" + str(i) + "E/-" + str(i) + " /g' " + procpath)
    except:
        print(f'Error in {filepath}')


def dict_levels_AS(keyword, levels_AS, levels, nsuplev=None, debug=False):
    levdic = {idx_as: None for idx_as in levels_AS.index}
    if nsuplev is not None: 
        print('=> warning: using superlevels predefined')
    for i in levels_AS.index:
        flag = 0
        value1 = levels_AS.loc[i]['CFG']
        # cfg1, ncfg1 = value1, nconfig(value1)
        if type(value1) == str:
            cfg1, ncfg1 = value1, nconfig(value1)
        else:
            cfg1, ncfg1 = config(value1), value1
        s1 = levels_AS.loc[i]['2S+1']
        l1 = levels_AS.loc[i]['L']
        g1 = 2 * levels_AS.loc[i]['J'] + 1
        if ncfg1 is None: 
            print(f'{i}: {value1} not found')
            sys.exit()
        for j in levels.index:
            cfg2 = levels.loc[j]['Conf']
            ncfg2 = nconfig(cfg2)
            s2 = levels.loc[j]['S']
            l2 = levels.loc[j]['L']
            g2 = levels.loc[j]['gi'] 
            if (ncfg1 == ncfg2) & (s1 == s2) & (l1 == l2) & (g1 == g2): 
                levdic[i] = j
                flag = 1
                if debug:
                    print(i, cfg1, ncfg1, s1, l1, g1)
                    print(j, cfg2, ncfg2, s2, l2, g2)
                    print(' ')
                break
            elif nsuplev is not None and i >= nsuplev:
                flag = 1
                isuplev = super_terms(cfg1) if keyword == 'T' else super_levels(cfg1)
                levdic[i] = isuplev
        if flag == 0:
            if debug:
                print(f'{i}, {cfg1} not found\n')
    print('... all done!\n')
    levdic = {k: v for k, v in levdic.items() if v is not None}
    return levdic


def config(number):
    # return [key for key, val in cfgs.items() if val == number][-1]
    return [key for key, val in cfgs_DB.items() if val == number][-1]


def nconfig(name):
    try:
        ncfg = cfgs[name]
    except:
        ncfg = None
    return ncfg


def super_terms(name):
    nlev = None
    if (name == '3s.9s') or (name == '3S1 9S1'): nlev = 55
    if (name == '3s.8d') or (name == '3S1 8D1'): nlev = 56
    if (name == '3s.8f') or (name == '3S1 8F1'): nlev = 56
    if (name == '3s.8g') or (name == '3S1 8G1'): nlev = 56
    if (name == '3s.8h') or (name == '3S1 8H1'): nlev = 56
    if (name == '3s.9p') or (name == '3S1 9P1'): nlev = 57
    if (name == '3s.9d') or (name == '3S1 9D1'): nlev = 59
    if (name == '3s.9f') or (name == '3S1 9F1'): nlev = 59
    if (name == '3s.9g') or (name == '3S1 9G1'): nlev = 59
    if (name == '3s.9h') or (name == '3S1 9H1'): nlev = 59
    if (name == '3s.9i') or (name == '3S1 9I1'): nlev = 59
    if (name == '3s.9k') or (name == '3S1 9K1'): nlev = 59
    if (name == '3s.10s') or (name == '3S110S1'): nlev = 58
    if (name == '3s.10d') or (name == '3S110D1'): nlev = 62
    if (name == '3s.10f') or (name == '3S110F1'): nlev = 62
    if (name == '3s.11d') or (name == '3S111D1'): nlev = 65
    if (name == '3s.11f') or (name == '3S111F1'): nlev = 65
    if (name == '3s.12d') or (name == '3S112D1'): nlev = 68
    if (name == '3s.12f') or (name == '3S112F1'): nlev = 68
    if (name == '3s.13d') or (name == '3S113D1'): nlev = 71
    if (name == '3s.13f') or (name == '3S113F1'): nlev = 71
    if (name == '3s.14d') or (name == '3S114D1'): nlev = 73
    if (name == '3s.14f') or (name == '3S114F1'): nlev = 73
    if (name == '3s.15d') or (name == '3S115D1'): nlev = 75
    if (name == '3s.15f') or (name == '3S115F1'): nlev = 75
    return nlev


def super_levels(name):
    nlev = None
    if (name == '3s.9s') or (name == '3S1 9S1'): nlev = 99
    if (name == '3s.8d') or (name == '3S1 8D1'): nlev = 100
    if (name == '3s.8f') or (name == '3S1 8F1'): nlev = 100
    if (name == '3s.8g') or (name == '3S1 8G1'): nlev = 100
    if (name == '3s.8h') or (name == '3S1 8H1'): nlev = 100
    if (name == '3s.9p') or (name == '3S1 9P1'): nlev = 101
    if (name == '3s.9d') or (name == '3S1 9D1'): nlev = 103
    if (name == '3s.9f') or (name == '3S1 9F1'): nlev = 103
    if (name == '3s.9g') or (name == '3S1 9G1'): nlev = 103
    if (name == '3s.9h') or (name == '3S1 9H1'): nlev = 103
    if (name == '3s.9i') or (name == '3S1 9I1'): nlev = 103
    if (name == '3s.9k') or (name == '3S1 9K1'): nlev = 103
    if (name == '3s.10s') or (name == '3S110S1'): nlev = 102
    if (name == '3s.10d') or (name == '3S110D1'): nlev = 106
    if (name == '3s.10f') or (name == '3S110F1'): nlev = 106
    if (name == '3s.11d') or (name == '3S111D1'): nlev = 109
    if (name == '3s.11f') or (name == '3S111F1'): nlev = 109
    if (name == '3s.12d') or (name == '3S112D1'): nlev = 112
    if (name == '3s.12f') or (name == '3S112F1'): nlev = 112
    if (name == '3s.13d') or (name == '3S113D1'): nlev = 115
    if (name == '3s.13f') or (name == '3S113F1'): nlev = 115
    if (name == '3s.14d') or (name == '3S114D1'): nlev = 117
    if (name == '3s.14f') or (name == '3S114F1'): nlev = 117
    if (name == '3s.15d') or (name == '3S115D1'): nlev = 119
    if (name == '3s.15f') or (name == '3S115F1'): nlev = 119
    return nlev
    


def print_unformatted(fpath, ups):
    ups.to_csv(fpath, index=False)


def print_formatted(fout, ups):
    print('print calculation...')
    fp = open(fout, "w+")
    zn = 12
    ionch = 0
    modidx = 1
    source = 999
    header = ["AtomicNumber", "IonCharge", "ModelIndex", "LowerLevel", "UpperLevel",
              "Temperature", "CollisionStrength", "Source"]
    print(*header, sep="\t", file=fp)
    temperature = default_temperature_grid()
    ntemp = len(temperature)
    for tran in ups.index:
        row = ups.loc[tran]
        i, k = int(row["i"]), int(row["k"])
        rups = row[2:].values
        for temp, value in zip(temperature, rups):
            if i < k: 
                    #print(zn,ionch,modidx,ni,nk,tjj,upsjj,source,sep="\t")
                print(zn, ionch, modidx, i, k, temp, value, source, sep = "\t", file = fp)
            if i > k: 
                print(i, '>', k)
                #print(zn,ionch,modidx,nk,ni,tjj,upsjj,source,sep="\t")
                print(zn, ionch, modidx, k, i, temp, value, source, sep = "\t", file = fp)

    print("... all print ok")
    return


def determine_parity(df, ip):
    df.insert(ip, 'P', -1)
    odd = df.index[df.loc[:]['2S+1'] < 0].tolist()
    even = df.index[df.loc[:]['2S+1'] > 0].tolist()
    for i in odd: df.at[i,'P'] = 1
    for i in even: df.at[i,'P'] = 0


def initskip_olg(olg_pathfile, trantype, key):
    if trantype == 'terms': wgrep = "'LIST OF TERMS WITH A WEIGHTED MEAN OVER THE FINE STRUCTURE' "
    if trantype == 'E1': wgrep = "'E1-DATA' "
    if trantype == 'E2': wgrep = "'E2' "
    if trantype == 'E3': wgrep = "'E3' "
    wsys = "grep -n -i " + wgrep + olg_pathfile + " > dum.txt"
    os.system(wsys)
    dumline = pd.read_csv("dum.txt", header=None, sep='\s+')
    nrows = len(dumline)
    if key == 'terms': iline = dumline.loc[0][0]
    if key == 'levels': iline = dumline.loc[nrows - 1][0]
    os.system("rm dum.txt")
    return int(iline[0:-1])


def quantnumber_to_termsymbol(sq, lq, pq):
    # multiplicity 2S+1
    chsq = str(int(sq))
    # angular momenta L
    if lq == 0: chlq = 'S'
    if lq == 1: chlq = 'P'
    if lq == 2: chlq = 'D'
    if lq == 3: chlq = 'F'
    if lq == 4: chlq = 'G'
    if lq == 5: chlq = 'H'
    if lq == 6: chlq = 'I'
    if lq == 7: chlq = 'K'
    # parity
    if pq == 0: chpq=''
    if pq != 0: chpq='*'
    chterm = chsq + chlq + chpq
    return chterm


def termsymbol_to_quantnumber(chterm):
    chsq = chterm[0]
    chlq = chterm[1]
    # multiplicity 2S+1
    sq = int(chsq)
    # angular momenta L
    if chlq == 'S': lq = 0
    if chlq == 'P': lq = 1
    if chlq == 'D': lq = 2
    if chlq == 'F': lq = 3
    if chlq == 'G': lq = 4
    if chlq == 'H': lq = 5
    if chlq == 'I': lq = 6
    if chlq == 'K': lq = 7
    # parity
    if (lq % 2) == 0: pq = 0
    if (lq % 2) != 0: pq = 1
    return sq, lq, pq
