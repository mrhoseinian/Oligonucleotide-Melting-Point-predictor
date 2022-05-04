##
import itertools
from MT_Calculator import Tm_NN
import pandas as pd
##
K = 11
Probe1 = 'GGACCTGATCG'
Probe2 = 'GCACCAGATCG'
Probe3 = 'GCACCTGATCG'
Probe4 = 'GCACCTGATTG'
Bases = 'ACGT'


##
# Generate the library
def generate_library(k=K):
    library = []
    for item in itertools.product(Bases, repeat=k):
        library.append(''.join(item))
    return library


##
Lib = generate_library()
print(len(Lib))
##
captured1, captured2, captured3, captured4 = {}, {}, {}, {}
for curr in Lib:
    try:
        temp = Tm_NN(seq=Probe1, c_seq=curr)
        if temp >= 27:
            captured1[curr] = int(temp)
    except:
        pass
    try:
        temp = Tm_NN(seq=Probe2, c_seq=curr)
        if temp >= 27:
            captured2[curr] = int(temp)
    except:
        pass
    try:
        temp = Tm_NN(seq=Probe3, c_seq=curr)
        if temp >= 27:
            captured3[curr] = int(temp)
    except:
        pass
    try:
        temp = Tm_NN(seq=Probe4, c_seq=curr)
        if temp >= 27:
            captured4[curr] = int(temp)
    except:
        pass


##
df1 = pd.DataFrame(data=captured1, index=["Melting Point"]).T
df1 = df1[df1["Melting Point"] < 100]
df1.to_excel("./Serial Imprinting/probe1.xlsx")

df2 = pd.DataFrame(data=captured2, index=["Melting Point"]).T
df2 = df2[df2["Melting Point"] < 100]
df2.to_excel("./Serial Imprinting/probe2.xlsx")

df3 = pd.DataFrame(data=captured3, index=["Melting Point"]).T
df3 = df3[df3["Melting Point"] < 100]
df3.to_excel("./Serial Imprinting/probe3.xlsx")

df4 = pd.DataFrame(data=captured4, index=["Melting Point"]).T
df4 = df4[df4["Melting Point"] < 100]
df4.to_excel("./Serial Imprinting/probe4.xlsx")

##
cap1_2 = captured1.keys() & captured2.keys()
cap1_3 = captured1.keys() & captured3.keys()
cap1_4 = captured1.keys() & captured4.keys()
cap2_3 = captured2.keys() & captured3.keys()
cap2_4 = captured2.keys() & captured4.keys()
cap3_4 = captured1.keys() & captured2.keys()

cap1_2_3 = captured1.keys() & captured2.keys() & captured3.keys()
cap1_2_4 = captured1.keys() & captured2.keys() & captured4.keys()
cap1_3_4 = captured1.keys() & captured3.keys() & captured4.keys()
cap2_3_4 = captured2.keys() & captured3.keys() & captured4.keys()

cap1_2_3_4 = captured1.keys() & captured2.keys() & captured3.keys() & captured4.keys()







##

