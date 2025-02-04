# PeptideCharge

Calculate the net charge of a small peptide or an amino acid sequence

The net-charge of the amino acid or the amino acid sequence at any given pH is calculated using the Henderson-Hasselback equation. 

The isoelectic point (pI) of any peptide is the pH at which the net-charge of the protein is 0.

As we know each amino acid has a charge associated with it and theoretically, summation of the individual charges of the amino acids should equal to the charge of the peptide. 

But practically, each amio acid, depeding on its nature, the number of free ions, its natures, etc. has additional effect(s) on the neighbouring amoino acid making the theoretical calculation of the amono acid dfficult. 

Moreover, the pH of the condition, chages the charges and therefore the effect(s) the amino acid(s) have on eachother.

The iterative net-charge calculation as described first by Bjellqvist, using a modified form of the Henderson-Hasselbach equation is to sum all the positive and negative charges of the sequence at any given pH.

By iteratively calculating the net-charge of a peptide at differnet pH conditions, the pH which the net-charge is 0 can be determined to be the theoretical pI of the peptide. 

NetCharge.py
1. Get general info about the amino acid. 
2. Calculate the charge of the amino acid or the net charge of an amino acid sequence.

netchargePI.py:
1. Identify the theoretical pI of any peptide

The example usage for both functions can be found in the test notebook. 