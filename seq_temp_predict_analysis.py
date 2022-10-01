data_file = open("sequences.csv", "r")

value_dict = {}
count = 0
for line in data_file:
    data_string = line.split(',')

    value_dict[data_string[1]] = data_string[0]


print(len(value_dict))