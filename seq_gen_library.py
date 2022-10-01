import math
import MT_Calculator as mt
from Bio.Seq import Seq
from numpy import equal
import random


## Name: base4
## Parameters: integer i and integer seq_length
## Function: Given i and seq_length, converts i to base 4 representation of seq_length digits long and returns as a string
def base4(i, seq_length):
    dna_string = []
    count = 0
    while i > 0:
        dna_string.append(i % 4)  ## Euclidean division to convert bases
        i = math.floor(i / 4)
        count = count + 1

    for i in range(0, seq_length - count):  ## Bit extension by appending zero so length of dna_string + count = seq_length
        dna_string.append(0)

    return "".join(map(str, dna_string))


## Name: base10
## Parameters: sequence
## Function: Given i and seq_length, converts i to base 4 representation of seq_length digits long and returns as a string
def base10(seq):
    sum = 0
    for i in range(len(seq)):
        sum += int(seq[i]) * (4**i)
    return sum


## Name: dna_string_convert
## Parameters: A string of an integer in base4 representation
## Output: String nucleotide base equivalent of given base4 string, ex 0123 --> ATCG
def dna_string_convert(dna_string):
    dna_string = dna_string.replace("0", "A")
    dna_string = dna_string.replace("1", "T")
    dna_string = dna_string.replace("2", "C")
    dna_string = dna_string.replace("3", "G")
    return dna_string


## Name: dna_string_invert
## Parameters: A string of an DNA sequence
## Output: Base 4 equivalent of given sequence string, ex ATCG --> 0123
def dna_string_invert(dna_string):
    dna_string = dna_string.replace("A", "0")
    dna_string = dna_string.replace("T", "1")
    dna_string = dna_string.replace("C", "2")
    dna_string = dna_string.replace("G", "3")
    return dna_string


def get_stack(string):
    if string == "GC":
        return 13
    if string == "CC" or string == "GG":
        return 11
    if string == "CG" or string == "AC" or string == "GT":
        return 10
    if string == "TC" or string == "CT" or string == "AG" or string == "GA":
        return 8
    if string == "TG" or string == "CA" or string == "AT":
        return 7
    if string == "TT" or string == "AA":
        return 5
    if string == "TA":
        return 4


def Tm_stacking(seq, dna, con):
    stack_param = 0
    for i in range(len(seq) - 1):
        stack_param += get_stack(seq[i : i + 2])
    return 7.35 * (float(stack_param) / len(seq)) + (17.34 * math.log(len(seq))) + (4.96 * math.log(con)) + (0.89 * math.log(dna)) - 25.42


## Name: generate_sequences
## Parameters: An integer seq_length
## Output: Array of all possible DNA sequences of length seq_length, and corresponding melting temperature as predicted by nearest neighbor thermodynamics
def gen_seq(seq_length):
    f = []
    for i in range(0, pow(4, seq_length)):
        seq = Seq(dna_string_convert(base4(i, seq_length)))
        f.append((seq, mt.Tm_NN(seq)))
    return f


## Name: random_gen_seq
## Parameters: An integer seq_length, and integer n
## Output: Array of length n of random sequences of length seq_length, and corresponding melting temperature as predicted by nearest neighbor thermodynamics
def random_gen_seq(seq_length, num):
    seq_array = []
    for i in range(0, num):
        seq = Seq("".join(random.choices(["A", "T", "C", "G"], k=seq_length)))
        seq_array.append((seq, mt.Tm_NN(seq)))
    return seq_array
