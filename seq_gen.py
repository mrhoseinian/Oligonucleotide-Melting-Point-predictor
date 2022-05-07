import math
import Bio
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp

## Name: base4
## Parameters: integer i and integer seq_length
## Function: Given i and seq_length, converts i to base 4 representation of seq_length digits long and returns as a string
def base4(i, seq_length) :
    dna_string = []
    count = 0
    while i > 0:
        dna_string.append(i % 4) ## Euclidean division to convert bases
        i = math.floor(i / 4) 
        count = count + 1 

    for i in range(0, seq_length - count) : ## Bit extension by appending zero so length of dna_string + count = seq_length
        dna_string.append(0)

    return "".join(map(str, dna_string))

## Name: dna_string_convert
## Parameters: A string of an integer in base4 representation
## Output: String nucleotide base equivalent of given base4 string, ex 0123 --> ATCG
def dna_string_convert(dna_string) : 
    dna_string = dna_string.replace("0", "A")
    dna_string = dna_string.replace("1", "T")
    dna_string = dna_string.replace("2", "C")
    dna_string = dna_string.replace("3", "G")
    return dna_string

## Name: generate_sequences
## Parameters: An integer seq_length
## Output: CSV format of all possible DNA sequences of length seq_length, and corresponding melting temperature as predicted by nearest neighbor thermodynamics
def generate_sequences(seq_length) :
    for i in range(0, pow(4, seq_length)) :
        sequence = Seq(dna_string_convert(base4(i, seq_length))) ## Map integer to DNA sequence, construct Sequence Object

        print('%s,%3.4f' %(sequence, MeltingTemp.Tm_NN(sequence))) ## Predict melting temp using nearest neighbor method, print sequence and temp in CSV format
    return