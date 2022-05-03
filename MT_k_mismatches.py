##
from complement_generator import mismatchfinder
from MT_Calculator import Tm_NN
import csv
import numpy as np
import pandas as pd

##
Probe = 'ACACCTTAATCACCGCTTCA'
Max_mismatches = 5

##
for i in range(1, Max_mismatches + 1):
    with open(f'MT_mismatch_results/{i}mismatches.csv', 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        mismathes = mismatchfinder(Probe, i)
        csvwriter.writerows(mismathes)

##
for i in range(1, Max_mismatches + 1):
    df = pd.read_csv(f'MT_mismatch_results/{i}mismatches.csv', sep='\n', header=None)
    df.columns = ['sequence']
    df['MT'] = " "
    for index in df.index:
        try:
            temp = Tm_NN(seq=Probe, c_seq=df['sequence'][index])
            df['MT'][index] = temp
        except:
            df['MT'][index] = np.nan
    df = df.dropna()
    df.to_csv(f'MT_mismatch_results/MT for {i}mismatches.csv')

