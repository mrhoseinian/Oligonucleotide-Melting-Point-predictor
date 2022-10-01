import numpy
import os.path
import biovec
from biovec.utils import split_ngrams

pv = biovec.models.ProtVec("sequence_data/new_p2v_6.fasta", corpus_fname="sequence_data/split_8_0.txt", n=2, size=100)

pv.save("saved_models/prot2vec/new_p2v_8.model")
