##
import csv
from Bio import Seq


def mismatchfinder(Probe, changes):
    answer = []
    Bases = 'ACGT'

    def generate_complements(dna, change, index):
        if not change:
            answer.append(dna)
            return
        change -= 1
        for i in range(index, len(dna) - change):
            curr = dna[i]
            for char in Bases:
                if char == curr:
                    continue
                seq = dna[:i] + char + dna[i + 1:]
                generate_complements(seq, change, i + 1)

    c_seq = Seq.Seq(Probe).complement()
    generate_complements(c_seq, changes, 0)
    return answer