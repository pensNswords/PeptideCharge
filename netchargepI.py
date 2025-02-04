import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import NetCharge

'''
The isoelectic point (pI) of any peptide is the pH at which the netc harge of the protein is 0.

As we know each amino acid has a charge associated with it and theoretically,
summation of the individual charges of the amino acids should equal to the charge of the peptide. 

But practically, each amio acid, depeding on its nature, the number of free ions, its natures, etc. has additional effect(s)
on the neighbouring amoino acid making the theoretical calculation of the amono acid dfficult. 

Moreover, the pH of the condition, chages the charges and therefore the effect(s) the amino acid(s) have on eachother.

The iterative net-charge calculation as described first by Bjellqvist, using a modified form of the Henderson-Hasselbach equation
is to sum all the positive and negative charges of the sequence at any given pH.

By iteratively calculating the net-charge of a peptide at differnet pH conditions, 
the pH which the net-charge is 0 can be determined to be the theoretical pI of the peptide. 

'''

def charge_plot(charge_data1, method=''):
    fig = plt.figure(figsize = (10,5))

    fig.suptitle('pH Vs Charge of Protein - ' + str(method))
    ax1 = fig.add_subplot(1, 1,1)

    plt.axhline(y=0, color = 'black', alpha = 0.3, linestyle = '--')
    plt.axvline(x=7, color = 'black', alpha = 0.3, linestyle = '--')

    ax1 = sns.lineplot(data = charge_data1, x= 'pH', y = 'Net-Charge', color = 'tab:blue')
    plt.xticks(np.arange(0, 14+1, 1.0))

    x1 = charge_data1['pH'].min()
    y1 = charge_data1.loc[charge_data1['pH'] == charge_data1['pH'].min(), 'Net-Charge'].iloc[0]
    plt.text(x1, y1+0.2, round(y1, 2), fontdict=dict(color='red',size=7))

    x2 = charge_data1['pH'].max()
    y2 = charge_data1.loc[charge_data1['pH'] == charge_data1['pH'].max(), 'Net-Charge'].iloc[0]
    plt.text(x2-0.5, y2+0.2, round(y2, 2), fontdict=dict(color='red',size=7))

    ch_7_1 = charge_data1.loc[(charge_data1['pH']-7).abs().argsort()[:1],"Net-Charge"].iloc[0]
    plt.text(7,ch_7_1+0.6 , round(ch_7_1,2), fontdict=dict(color='red',size=7))

    plt.tight_layout()
    plt.show()



def pH_sim():
    ini = [7]
    t = 7
    ls = []

    for i in range(0,10):
        for p in ini:
            ph1 = p + t/2
            if ph1 not in ini:
                ls.append(ph1)
                print(ph1)
            ph2 = p - t/2
            if ph2 not in ini:
                ls.append(ph2)
                print(ph2)
        ini = ini + list(set(ls) - set(ini))
        t = t/2
        i = i+1
    ini.sort()
    return(ini)


def ncpI(aa, pI_plot = True):
    ph_list = pH_sim()
    nc = []
    for ph in ph_list:
        sum_chg = NetCharge.chain_net_charge(aa, ph)
        nc.append(sum_chg)
    
    df = pd.DataFrame(ph_list, nc, columns = ["pH", "Net-Charge"])

    theo_pI = np.mean(df.iloc[(df['Net-Charge']-0).abs().argsort()[:2]]['pH'])

    if pI_plot == True:
        charge_plot(charge_data1=df, method='Iterative Net-Charge Method')
    return (theo_pI)

    