import math
import Bio
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp

def base4(i, seq_length) :
    dna_string = []
    count = 0
    while i > 0:
        dna_string.append(i % 4)
        i = math.floor(i / 4)
        count = count + 1

    for i in range(0, seq_length - count) :
        dna_string.append(0)

    dna_string.reverse()

    return "".join(map(str, dna_string))

def dna_string_convert(dna_string) :
    dna_string = dna_string.replace("0", "A")
    dna_string = dna_string.replace("1", "T")
    dna_string = dna_string.replace("2", "C")
    dna_string = dna_string.replace("3", "G")
    return dna_string

seq_length = int(input("Enter sequence length: "))

print(seq_length)


for i in range(0, pow(4, seq_length)) :
    sequence = Seq(dna_string_convert(base4(i, seq_length)))

    print('%s,%3.4f' %(sequence, MeltingTemp.Tm_NN(sequence)))
