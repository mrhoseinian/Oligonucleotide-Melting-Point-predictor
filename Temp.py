##
import random
# import itertools
# import sys
from math import ceil
from math import floor
from Bio.SeqUtils import MeltingTemp as MT
import pandas as pd

# from Bio.Seq import Seq

# import numpy as np

##
# K = int(input("K=?"))
K = 20
# K = int(sys.argv[1])
"""For DNA"""
Bases = 'ACGT'
"""For RNA"""

# Bases = 'ACGU'


# Bases = '0123'


##
# # generating the library
# def generate_library(k):
#     library = []
#     for item in itertools.product(Bases, repeat=k):
#         library.append(''.join(item))
#     return library


##
# generate a random dna
# def random_dna_and_rna(length):
#     DNA = ""
#     for count in range(length):
#         DNA += random.choice(Bases)
#     return DNA
DNA = "ACACCTTAATCACCGCTTCA"


##
# a function for adding 1-edit distance sequences
def add_sequences(dna):
    sequences = []
    current_sequence = list(dna)
    """WITHOUT TERMINAL MISMATCH"""
    # for char in range(1, len(dna) - 1):
    """WITH TERMINAL MISMATCH"""
    for char in range(len(dna)):
        temp_list = list(current_sequence)
        for base in Bases:
            temp_list[char] = base
            sequences.append(''.join(temp_list))
        sequences.remove(dna)

    return sequences


##
# a function for adding 2-edit distance sequences
def add_sequences_again(dna, offset):
    sequences = []
    current_sequence = list(dna)
    """WITHOUT TERMINAL MISMATCH"""
    # for char in range(offset + 1, len(dna) - 1):
    """WITH TERMINAL MISMATCH"""
    for char in range(offset + 1, len(dna)):
        temp_list = list(current_sequence)
        for base in Bases:
            temp_list[char] = base
            sequences.append(''.join(temp_list))
        sequences.remove(dna)

    return sequences


##
# finding the complementary sequence
def comp_seq(candidate):
    """For DNA"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    """For RNA"""
    # complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
    bases = list(candidate)
    bases = [complement[base] for base in bases]
    return ''.join(bases)


##
# def melting_temp(sequences):
#     temperatures = []
#     for sequence in sequences:
#         temperature = mt.Tm_NN('GGCCA', c_seq=sequence)
#         temperatures.append(temperature)
#     return temperatures
#

##
# Probe = random_dna_and_rna(K)
Probe = DNA
# 'CCGGCGCATACGTCTCGGTC'
Target = comp_seq(Probe)
##
One_edit_seq = add_sequences(Target)
# Temp_list = melting_temp(temp)

##
Two_edit_seq = []
"""without terminal mismatch"""
# for index, item in enumerate(One_edit_seq):
#     if ceil((index + 1) / 3) == (len(item) - 2):
#         continue
#     Two_edit_seq.append(add_sequences_again(item, ceil((index + 1) / 3)))

"""with terminal mismatch"""
for index, item in enumerate(One_edit_seq):
    if ceil((index + 1) / 3) == (len(item)):
        continue
    Two_edit_seq.append(add_sequences_again(item, floor(index / 3)))

# ##
# """WRITING TO A FILE"""
# with open("E:/EDU/Research/Coding/R/R/data/list.txt", "w") as file:
#     file.write(Candidate + "\n")
#     file.write(Comp_candidate + "\n")
#     for item in One_edit_seq:
#         file.write(item + "\n")
# ##
# with open("E:/EDU/Research/Coding/R/R/data/list_two.txt", "w") as file:
#     for list_item in Two_edit_seq:
#         for item in list_item:
#             file.write(item + "\n")

##
"""COMPUTING THE MELTING TEMPERATURE"""
MELTING_TEMP = {}
MELTING_TEMP_2MISS = {}
MELTING_TEMP[Target] = MT.Tm_NN(seq=Probe, c_seq=Target, nn_table=MT.DNA_NN4, dnac1=100000, dnac2=100000,
                                Mg=1, Tris=1, Na=0, K=50)
for one_seq in One_edit_seq:
    MELTING_TEMP[one_seq] = MT.Tm_NN(seq=Probe, c_seq=one_seq, nn_table=MT.DNA_NN4, dnac1=100000, dnac2=100000,
                                     Mg=1, Tris=1, Na=0, K=50)
##
for two_seq_list in Two_edit_seq:
    for two_seq in two_seq_list:
        try:
            MELTING_TEMP_2MISS[two_seq] = MT.Tm_NN(seq=Probe, c_seq=two_seq, nn_table=MT.DNA_NN4, dnac1=100000,
                                                   dnac2=100000, Mg=1, Tris=1, Na=0, K=50)

        except:
            MELTING_TEMP_2MISS[two_seq] = None
##
""""WRITE TO A FILE"""
df = pd.DataFrame(data=MELTING_TEMP, index=[0])
df = df.T
df = df.set_axis(["Melting Point"], axis=1)
df.to_excel("~/EDU/Research/Coding/Melting point result/Winter 22/One_Miss.xlsx")
df = pd.DataFrame(data=MELTING_TEMP_2MISS, index=[0])
df = df.T
df = df.set_axis(["Melting Point"], axis=1)
df.to_excel("~/EDU/Research/Coding/Melting point result/Winter 22/Two_Miss.xlsx")
##
# # Test here!
# print('%f' % MT.Tm_NN('CCGGCGCATACGTCTCGGT', c_seq='GGCCGCGTATGCAGAGCCA', nn_table=MT.DNA_NN4, dnac1=250,
#                       dnac2=250, Mg=1, saltcorr=1))

##
# Let's make one list of all!
temp_dict = {Target: MT.Tm_NN(seq=Probe, c_seq=Target, nn_table=MT.DNA_NN4, dnac1=100000,
                              dnac2=100000, Mg=1, Tris=1, Na=0, K=50)}
for one_seq in One_edit_seq:
    temp_dict[one_seq] = MT.Tm_NN(seq=Probe, c_seq=one_seq, nn_table=MT.DNA_NN4, dnac1=100000,
                                  dnac2=100000, Mg=1, Tris=1, Na=0, K=50)
for two_seq_list in Two_edit_seq:
    for two_seq in two_seq_list:
        try:
            temp_dict[two_seq] = MT.Tm_NN(seq=Probe, c_seq=two_seq, nn_table=MT.DNA_NN4, dnac1=100000,
                                          dnac2=100000, Mg=1, Tris=1, Na=0, K=50)

        except:
            temp_dict[two_seq] = None
##
temp_3miss = []
for two_miss_list in Two_edit_seq:
    for sequnce_of_list in two_miss_list:
        temp_3miss.append(add_sequences(sequnce_of_list))
for list_item in temp_3miss:
    for item in list_item:
        if item not in temp_dict:
            try:
                temp_dict[item] = MT.Tm_NN(seq=Probe, c_seq=item, nn_table=MT.DNA_NN4, dnac1=100000,
                                           dnac2=100000, Mg=1, Tris=1, Na=0, K=50)
            except:
                temp_dict[item] = None
##
df = pd.DataFrame(data=temp_dict, index=[0])
df = df.T
df = df.set_axis(["Melting Point"], axis=1)
df.to_excel("~/EDU/Research/Coding/Melting point result/Winter 22/MeltingPoint.xlsx")

## Four mismatches for one specific sequence: TTTGAAATTAGTGGCGAAAT
Four_edit_seq = add_sequences('TTTGAAATTAGTGGCGAAAT')
MELTING_TEMP_4MISS = {}
for seq in Four_edit_seq:
    try:
        MELTING_TEMP_4MISS[seq] = MT.Tm_NN(seq=Probe, c_seq=seq, nn_table=MT.DNA_NN4, dnac1=100000,
                                           dnac2=100000, Mg=1, Tris=1, Na=0, K=50)

    except:
        MELTING_TEMP_4MISS[seq] = None
df = pd.DataFrame(data=MELTING_TEMP_4MISS, index=[0])
df = df.T
df = df.set_axis(["Melting Point"], axis=1)
df.to_excel("~/EDU/Research/Coding/Melting point result/Winter 22/Four_Miss_for_one_sequence.xlsx")

##

