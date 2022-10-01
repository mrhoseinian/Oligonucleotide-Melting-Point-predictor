from cgi import test
import numpy
import xgboost as xgb
import biovec
import warnings
import protvec_helper_library as p2v

from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn import metrics as met
import joblib


def train_and_score_reg(test_name, X, result_file):

    X_train, X_test, y_train, y_test = train_test_split([x for x in numpy.array(X)[:, 1]], numpy.array(X)[:, 2], test_size=0.25, random_state=11)
    print("Train Regressor")

    model = xgb.XGBRegressor(random_state=11, tree_method="hist")

    model.fit(X_train, y_train)

    model.save_model("saved_models/" + test_name + ".model")

    model_performance(test_name, y_test, model.predict(X_test), result_file)
    print("Finished Train Regressor")

    return model, X_test, y_test


def train_tri_region_split(test_name, X, t1, t2, result_file):
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

    model1 = train_and_score_reg(test_name + "_low_region", X_low, result_file)
    model2 = train_and_score_reg(test_name + "_lin_region", X_lin, result_file)
    model3 = train_and_score_reg(test_name + "_high_region", X_high, result_file)

    return model1, model2, model3


def train_predict_with_splitter(test_name, X, t1, t2, result_file):

    print("Splitting")

    X_train, X_test, y_train, y_test = train_test_split(numpy.array(X), numpy.array(X)[:, 2], test_size=0.25, random_state=11)

    clf = train_and_score_clf(test_name, X, t1, t2, result_file)

    model1, model2, model3 = train_tri_region_split(test_name, X_train, t1, t2, result_file)

    X_low = []
    y_low = []
    X_high = []
    y_high = []
    X_lin = []
    y_lin = []

    X_predict = clf.predict([x for x in X_test[:, 1]])
    print("Sort into bins")
    count = 0

    for i in X_test:
        if X_predict[count] == 0:
            X_low.append(i[1])
            y_low.append(i[2])
        if X_predict[count] == 1:
            X_lin.append(i[1])
            y_lin.append(i[2])
        if X_predict[count] == 2:
            X_high.append(i[1])
            y_high.append(i[2])
        count += 1

    print("Eval performance")
    model_performance(test_name + "_low_region_clf_predict", y_low, model1[0].predict(X_low), result_file)
    model_performance(test_name + "_lin_region_clf_predict", y_lin, model2[0].predict(X_lin), result_file)
    model_performance(test_name + "_high_region_clf_predict", y_high, model3[0].predict(X_high), result_file)
    return


def train_and_score_clf(test_name, X, t1, t2, result_file):

    print("Train Classifier")
    y = []
    for x in X[:, 2]:
        y.append(p2v.classify(x, t1, t2))

    X_train, X_test, y_train, y_test = train_test_split([x for x in X[:, 1]], numpy.array(y), test_size=0.25, random_state=11)

    classifier = KNeighborsClassifier()

    classifier.fit(numpy.array(X_train), numpy.array(y_train))

    joblib.dump(classifier, "saved_models/" + test_name + "_" + str(t1).replace(".", "_") + "_" + str(t2).replace(".", "_") + "_classifier.model")

    model_performance(test_name + "_classifier", y_test, classifier.predict(X_test), result_file)

    print("Finished Train Classifier")

    return classifier


def model_performance(test_name, y_test, y_predict, result_file):
    r2_score = met.r2_score(y_test, y_predict)
    mape_score = met.mean_absolute_percentage_error(y_test, y_predict)
    mae_score = met.mean_absolute_error(y_test, y_predict)
    rmse_score = met.mean_squared_error(y_test, y_predict)

    results_to_file(test_name, r2_score, mape_score, mae_score, rmse_score, result_file)


def results_to_file(test_name, r2, mape, mae, rmse, result_file):
    file_results = open("results/" + result_file + "_results.csv", "a")
    file_results.write("%s, %3.5f, %3.5F, %3.5F, %3.5f\n" % (test_name, r2, mape, mae, rmse))
    file_results.close()
