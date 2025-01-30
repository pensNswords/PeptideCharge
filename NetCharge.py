import numpy as np


# Master Amino Acid Dictionary
'''
Based on Handbook of Chemistry and Physics
Ref: Reference: D.R. Lide, Handbook of Chemistry and Physics, 72nd Edition, CRC Press, Boca Raton, FL, 1991.
'''
aa_names = {"A": "Alanine", "R": "Arginine", "N": "Aspargine", "D": "Aspartic Acid", "C": "Cysteine", "Q": "Glutamine", "E": "Glutamic Acid", 
            "G": "Glycine", "H": "Histidine", "I": "Isoleucine", "L": "Leucine", "K": "Lysine", "M": "Methionine", "F": "Phenylalanine",
            "P": "Proline", "S": "Serine", "T": "Threonine", "W": "Tryptophan", "Y": "Tyrosine", "V": "Valine"}

aa_pI = {"A": 6.00, "R": 10.76, "N": 5.41, "D": 2.77, "C": 5.07, "Q": 5.65, "E": 3.22, "G": 5.97, "H": 7.56, "I": 6.02, "L": 5.98, "K": 9.74,
        "M": 5.74, "F": 5.48, "P": 6.30, "S": 5.68, "T": 5.60, "W": 5.89, "Y": 5.66, "V": 5.96}

pka_n_avg = 9.69
pka_c_avg = 2.34

# Values Based on D.R. Lide, Handbook of Chemistry and Physics, 72nd Edition, CRC Press, Boca Raton, FL, 1991.
pka_f = {"C": 8.18, "D": 3.65, "E": 4.25, "H": 6.00, "K": 10.53, "R": 12.48, "Y": 10.07}

pka_c = {"A": 2.34, "R": 2.17, "N": 2.02, "D": 1.88, "C": 1.96, "Q": 2.17, "E": 2.19, "G": 2.34, "H": 1.82, "I": 2.36, "L": 2.36, "K": 2.18,
         "M": 2.28, "F": 1.83, "P": 1.99, "S": 2.21, "T": 2.09, "W": 2.83, "Y": 2.20, "V": 2.32}

pka_n = {"A": 9.69, "R": 9.04, "N": 8.80, "D": 9.60, "C": 10.28, "Q": 9.13, "E": 9.67, "G": 9.60, "H": 9.17, "I": 9.60, "L": 9.60, "K": 8.95,
         "M": 9.21, "F": 9.13, "P": 10.60, "S": 9.15, "T": 9.10, "W": 9.39, "Y": 9.11, "V": 9.62}

aa_mass = {"A": 89.094, "R": 174.203, "N": 132.119, "D": 113.104, "C": 121.154, 
            "Q": 146.146, "E": 147.131, "G": 75.067, "H": 155.156, "I": 131.175, 
            "L": 131.175, "K": 146.189, "M": 165.192, "F": 165.192, "P": 115.132, 
            "S": 105.093, "T": 119.119, "W": 204.228, "Y": 181.191, "V": 117.148}

aa_3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", 
         "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
         "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", 
         "Y": "Tyr", "V": "Val"}


aa_type = {"A": "Non-Polar, Aliphatic, Hydrophobic", 
           "R": "Polar, Basic", 
           "N": "Polar", 
           "D": "Polar, Acidic", 
           "C": "Polar", 
           "Q": "Polar", 
           "E": "Polar, Acidic", 
            "G": "Non-Polar, Aliphatic", 
            "H": "Polar, Basic, Essential", 
            "I": "Non-Polar, Hydrophobic, Aliphatic, Essential", 
            "L": "Non-Polar, Aliphatic, Essential", 
            "K": "Polar, Basic, Essential", 
            "M": "Non-Polar, Hydrophobic, Essential", 
            "F": "Non-Polar, Hydrophobic, Aromatic, Essential",
            "P": "Non-Polar, Cyclic", 
            "S": "Polar, Hydroxyl", 
            "T": "Polar, Hydroxyl, Essential", 
            "W": "Non-Polar, Hydrophobic, Aromatic, Essential", 
            "Y": "Polar, Aromatic, Hydrophobic", 
            "V": "Non-Polar, Hydrophobic, Aliphatic, Essential"
           }



def aa_info(amino_acid):    
    if len(amino_acid) > 1:
        try:
            aone = aa_names()[aa_names.values().index(amino_acid)]
        except Exception as e1:
            print("Enter a valid natural aminoa acid")
            print(e1)
        aful = amino_acid
    else:
        aone = amino_acid
        try:
            aful = aa_names[amino_acid]
        except Exception as e2:
            print("Enter a valid natural aminoa acid")
            print(e2)

    info = {"Amino Acid": aful,
            "Abbreviation": aa_3[aone], 
            "One-Letter Code": aone,
            "Molar Mass": aa_mass[aone],
            "Amino Acid pI": aa_pI[aone],
            "General Attributes": aa_type[aone] }
    
    return info



# Amino Acid Charge Calculator
'''
Formula Ref: https://www.bachem.com/knowledge-center/peptide-calculator/
The Formula is a simplified version of the Henderson-Hasselbach equatiion [A staple in chemistry]

# Simplified Formula
# Ref: Cargile B.J., et al. (2008) Calculation of the isoelectric point of tryptic peptides in the pH 3.5–4.5 range based on adjacent amino acid effects. Electrophoresis, 29, 2768–2778

The Net-Charge method involves summing the N-Terminal and C-Terminal charges alsong with the charges contributed by the Funcional groups of the acidic anf basic amino acids

aa: Amino Acid
ph: ph (change is calculated for a given pH)
n: the pKa value of the N-terminal amino acid
c: the pKa value of the C-terminal aminoa acid
f: the pKa calue of the amino acid with other functional group

for any individual aminoa acid, there will only be either n, c or f value where others are set to 0
'''

def charge_cal(aa, ph, n, c, f):
    sum = 0
    term_n = 1 / (np.power(10, (ph-n)) + 1) if n !=0 else 0

    term_c = - 1 / ((np.power(10, (c-ph))) + 1) if c !=0 else 0

    
    if f != 0:
        if aa in ["R", "K", "H"]:
            term_f = (1) / ((np.power(10, (ph-f))) + 1)
        elif aa in ["D", "E", "C", "Y"]:
            term_f = - (1) / ((np.power(10, (f-ph))) + 1)
        else:
            term_f = 0
    else:
        term_f = 0
    
    sum = term_n + term_c + term_f
    return (sum)


#Net Charge of aa Peptide (an amino acid chain);
'''
While there are amany methods to calculate the peptide net charge, the most common and the easiest one is to sum the individual amino acid charges. 

chain: amino acid chain. Needs to be in the form of a string of Amino acids or a list with the one-letter code of each amino acid being separated by a comma
ph: pH at which the charge needs to be calculated
'''

def chain_net_charge(chain, ph):

    try:

        if type(chain) not in [str, list]:
            raise ValueError("The input needs to be in the form of string of amino acids or a lits of Amino Acid sequence")

        if isinstance(chain, list) == False:
            chain = list(chain)

        if set(chain) - set(list(aa_names.keys())) != set():
            raise ValueError("The input needs to have a valid natural amino acid")

        start = chain[0]
        end = chain[-1]
        chain_n = pka_n[start]
        chain_c = pka_c[end]
        chain_f_start = pka_f[start] if start in ["R", "K", "H","D", "E", "C", "Y"] else 0
        chain_f_end = pka_f[end] if end in ["R", "K", "H","D", "E", "C", "Y"] else 0

        start_chg = charge_cal(aa = start, ph=ph, n= chain_n,c=0,f=chain_f_start)
        end_chg = charge_cal(aa=end, ph=ph, n=0, c=chain_c, f=chain_f_end)

        sum_chg = start_chg + end_chg

        for i in range (1,(len(chain)-1)):
            if chain[i] in ["R", "K", "H","D", "E", "C", "Y"]:
                chain_f = pka_f[chain[i]]
                chg = charge_cal(aa = chain[i], ph = ph, n = 0, c = 0, f = chain_f)
            else:
                chg = 0

            sum_chg = sum_chg + chg
        return sum_chg
    
    except Exception as e:
        print("The error is: ", e)