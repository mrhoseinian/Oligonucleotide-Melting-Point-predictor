import biovec

pv = biovec.models.ProtVec('sequences.txt', corpus_fname='split_sequences.txt', n=2)

pv.save('p2v_model.model')







