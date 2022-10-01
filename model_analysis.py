import protvec_helper_library as p2v
import model_helper_library as tester
import seq_gen_library as gen
import numpy


test_name = str(input("Specify test name: "))

new_data = str(input("New Seq Data Gen? [Y or N]: "))


if new_data == "Y":
    format = input("Format? (CSV or FASTA)")
    variation = str(input("Varied Seq length? [Y or N]: "))
    if variation == "Y":
        seq_length1 = int(input("Lower Bound for Seq length: "))
        seq_length2 = int(input("Upper Bound for Seq length: "))
        num = int(input("Number of each sequence: "))
        for i in range(seq_length1, seq_length2):
            p2v.data_to_file(gen.random_gen_seq(i, num), test_name, format)
    else:
        seq_length = int(input("Seq length: "))
        p2v.data_to_file(gen.gen_seq(seq_length), test_name, format)

data_file = input("Data File Name? ")


X = p2v.data_from_csv(data_file)

# t1 = float(input("Specify low region cutoff temp: "))
# t2 = float(input("Specify high region cutoff temp:  "))

# tester.train_predict_with_splitter(test_name, X, t1, t2, test_name)
model, X_test, y_test = tester.train_and_score_reg(test_name, X, test_name)

for i in range(1, 8):
    X = p2v.data_from_array(gen.gen_seq(i))

    y_predict = model.predict([x for x in X[:, 1]])

    tester.model_performance(test_name, [x for x in X[:, 2]], y_predict, test_name)

for i in range(8, 200):
    X = p2v.data_from_array(gen.random_gen_seq(i, 3000))

    y_predict = model.predict([x for x in X[:, 1]])

    tester.model_performance(test_name, [x for x in X[:, 2]], y_predict, test_name)

for i in range(300, 500):
    X = p2v.data_from_array(gen.random_gen_seq(i, 1000))

    y_predict = model.predict([x for x in X[:, 1]])

    tester.model_performance(test_name, [x for x in X[:, 2]], y_predict, test_name)
