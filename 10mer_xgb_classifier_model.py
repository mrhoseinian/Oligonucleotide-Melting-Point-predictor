import numpy
import biovec
import random
import xgboost as xgb
import warnings
import protvec_helper_library as p2v
import seq_gen_library as gen

from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn import metrics as met

warnings.filterwarnings("ignore")

pv = biovec.models.load_protvec("saved_models/prot2vec/saved_model.model")
seq_file = open("sequence_data/10mer_data.csv", "r")
data_X = []

print("Preparing data_X from file...")

for line in seq_file:
    arr = line.split(",")
    vec = numpy.concatenate(p2v.to_vecs(pv, 2, arr[0], 15))
    data_X.append((arr[0], vec, float(arr[1])))

seq_file.close()

print("Finished converting seq to feat vectors, classifying data_X...")

X = numpy.array(data_X)
y_class = []
t1 = 37.3
t2 = 67.6

for x in X[:, 2]:
    y_class.append(p2v.classify(x, t1, t2))

X_class_train, X_class_test, y_class_train, y_class_test = train_test_split([x for x in X[:, 1]], numpy.array(y_class), random_state=11)

classifier = KNeighborsClassifier()

classifier.fit(X_class_train, y_class_train)

X_low = []
X_high = []
X_lin = []
for i in X:
    if p2v.classify(i[2], t1, t2) == 0:
        X_low.append(i)
    if p2v.classify(i[2], t1, t2) == 1:
        X_lin.append(i)
    if p2v.classify(i[2], t1, t2) == 2:
        X_high.append(i)

print("Classifying data_X complete, splitting data set... ")

X_low_train, X_low_test, y_low_train, y_low_test = train_test_split(numpy.array(X_low), numpy.array(X_low)[:, 2], test_size=0.25, random_state=11)
X_lin_train, X_lin_test, y_lin_train, y_lin_test = train_test_split(numpy.array(X_lin), numpy.array(X_lin)[:, 2], test_size=0.25, random_state=12)
X_high_train, X_high_test, y_high_train, y_high_test = train_test_split(numpy.array(X_high), numpy.array(X_high)[:, 2], test_size=0.25, random_state=13)

X_red_low_train = [x for x in X_low_train[:, 1]]
X_red_low_test = [x for x in X_low_test[:, 1]]

X_red_lin_train = [x for x in X_lin_train[:, 1]]
X_red_lin_test = [x for x in X_lin_test[:, 1]]

X_red_high_train = [x for x in X_high_train[:, 1]]
X_red_high_test = [x for x in X_high_test[:, 1]]

print("Fitting model...")

low = xgb.XGBRegressor(random_state=11, tree_method="hist")
lin = xgb.XGBRegressor(random_state=12, tree_method="hist")
high = xgb.XGBRegressor(random_state=13, tree_method="hist")

low.fit(X_red_low_train, y_low_train)
lin.fit(X_red_lin_train, y_lin_train)
high.fit(X_red_high_train, y_high_train)

low.save_model("10mer_low_region.model")
lin.save_model("10mer_lin_region.model")
high.save_model("10mer_high_region.model")

low_predict = low.predict(X_red_low_test)
lin_predict = lin.predict(X_red_lin_test)
high_predict = high.predict(X_red_high_test)


print("Low Region R^2: %3.4f" % (met.r2_score(y_low_test, low_predict)))
print("Lin Region R^2: %3.4f" % (met.r2_score(y_lin_test, lin_predict)))
print("High Region R^2: %3.4f" % (met.r2_score(y_high_test, high_predict)))

print("Low Region MAPE: %3.4f" % (met.mean_absolute_percentage_error(y_low_test, low_predict)))
print("Lin Region MAPE: %3.4f" % (met.mean_absolute_percentage_error(y_lin_test, lin_predict)))
print("High Region MAPE: %3.4f" % (met.mean_absolute_percentage_error(y_high_test, high_predict)))

print("Low Region MAE: %3.4f" % (met.mean_absolute_error(y_low_test, low_predict)))
print("Lin Region MAE: %3.4f" % (met.mean_absolute_error(y_lin_test, lin_predict)))
print("High Region MAE: %3.4f" % (met.mean_absolute_error(y_high_test, high_predict)))

print("Low Region RMSE: %3.4f" % ((met.mean_squared_error(y_low_test, low_predict) ** 0.5)))
print("Lin Region RMSE: %3.4f" % ((met.mean_squared_error(y_lin_test, lin_predict) ** 0.5)))
print("High Region RMSE: %3.4f" % ((met.mean_squared_error(y_high_test, high_predict) ** 0.5)))
