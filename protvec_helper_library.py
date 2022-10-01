import numpy
import os.path
import biovec
from biovec.utils import split_ngrams


def classify(num, t1, t2):
    if num <= t1:
        return 0
    elif num <= t2:
        return 1
    else:
        return 2


def to_vecs(pv, n=3, seq="ATATATATAT", size=15):
    ngram_patterns = split_ngrams(seq, n)
    protvecs = numpy.zeros([2, size])
    for num in range(1):
        ngram_vecs = []
        for ngram in ngram_patterns[num]:
            try:
                ngram_vecs.append(pv[ngram])
            except:
                raise Exception("Model has never trained this n-gram: " + ngram)
        protvecs[
            num,
        ] = sum(ngram_vecs)
    return protvecs


def multi_to_vecs(pv, seq_array, n=3, size=15):
    vec_array = []
    for seq in seq_array:
        vec_array.append(numpy.concatenate(to_vecs(pv, n, seq, size)))
    return vec_array


def data_from_csv(data_file="10mer_data.csv", pv_model_file="saved_model.model"):
    pv = []
    seq_file = []
    if os.path.exists("saved_models/prot2vec/" + pv_model_file):
        print("P2v File found, using " + pv_model_file)
        pv = biovec.models.load_protvec("saved_models/prot2vec/" + pv_model_file)
    else:
        print("P2V Model File not found, using saved_model.model")
        pv = biovec.models.load_protvec("saved_models/prot2vec/saved_model.model")

    if os.path.exists("sequence_data/" + data_file):
        print("Data File found, using " + data_file)
        seq_file = open("sequence_data/" + data_file, "r")
    else:
        print("Data File not found, using 10mer_data.csv")
        seq_file = open("sequence_data/10mer_data.csv", "r")

    data_X = []

    for line in seq_file:
        arr = line.split(",")
        vec = numpy.concatenate(to_vecs(pv, 2, arr[0], 15))
        data_X.append((arr[0], vec, float(arr[1])))

    seq_file.close()
    return numpy.array(data_X)


def data_from_array(data_array, pv_model_file="saved_model.model"):
    pv = []
    if os.path.exists("saved_models/prot2vec/" + pv_model_file):
        print("P2v File found, using " + pv_model_file)
        pv = biovec.models.load_protvec("saved_models/prot2vec/" + pv_model_file)
    else:
        print("P2V Model File not found, using saved_model.model")
        pv = biovec.models.load_protvec("saved_models/prot2vec/saved_model.model")

    data_X = []

    for item in data_array:
        vec = numpy.concatenate(to_vecs(pv, 2, item[0], 15))
        data_X.append((item[0], vec, float(item[1])))
    return numpy.array(data_X)


def data_to_file(data_array, test_name, format):
    if format == "CSV":
        test_name = test_name + ".csv"
    elif format == "FASTA":
        test_name = test_name + ".fasta"
    else:
        test_name = test_name + ".txt"

    f = open("sequence_data/" + test_name, "a")

    count = 0
    for i in data_array:
        if format == "CSV":
            f.write("%s, %3.4f\n" % (i[0], i[1]))
        elif format == "FASTA":
            f.write(">Seq%s\n%s\n" % (count, i[0]))
        else:
            f.write(i)
        count += 1
    f.close()
