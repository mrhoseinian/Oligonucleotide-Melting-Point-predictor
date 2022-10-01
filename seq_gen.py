import math
import Bio
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp
from numpy import equal

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
def generate_sequences(seq_length, format) :
    f = open("sequences.txt", "w")

    for i in range(0, pow(4, seq_length)) :
        sequence = Seq(dna_string_convert(base4(i, seq_length))) ## Map integer to DNA sequence, construct Sequence Object
        
        if format == "CSV":
            f.write('%s,%3.4f\n' %(sequence, MeltingTemp.Tm_NN(seq=sequence, nn_table=MeltingTemp.DNA_NN4, dnac1=100000, dnac2=100000, Mg=1, Tris=1, Na=0, K=50))) ## Predict melting temp using nearest neighbor method, print sequence and temp in CSV format
        elif format == "FASTA":
            f.write('>Seq%s\n%s\n' %(i, sequence))
            if(i < pow(4,seq_length) - 1):
                f.write('\n')
        else:
            f.write('%s' %(sequence))

    f.close()
    return