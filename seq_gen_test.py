import seq_gen

num = int(input("Seq length? "))
format = input("Format? (CSV or FASTA)")

seq_gen.generate_sequences(num, format)