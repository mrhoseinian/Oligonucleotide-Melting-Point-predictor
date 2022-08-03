import numpy
import pandas
import sklearn

from sklearn.preprocessing import OrdinalEncoder

encoder = OrdinalEncoder()

data_file = open("split_sequences.txt", "r")

unprocessed_input = []

count = 0

for line in data_file:
    for i in range(0, 3):
        composite_string = composite_string + data_file.readline() 
    
    composite_string.replace('\n','')
    composite_string.replace(' ', '')
    
    unprocessed_input.append(composite_string)
    
    composite_string = ""

encoder.fit_transform(unprocessed_input)

