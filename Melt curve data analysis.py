##
import string

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()

##
Copies = 3
Concentrations = 3
Pair = 5
alphabet = list(string.ascii_uppercase)
alphabet = alphabet[:Pair]
Well_plate = []
for i in alphabet:
    for j in range(1, Copies + 1):
        Well_plate.append(i + str(j))
##
melting_data = pd.read_csv('./data/Melt Curve RFU Results.csv', header=None)
melting_data.dropna(inplace=True)
print(melting_data.head())

##
Count = len(melting_data.columns) - 1
Temperature = melting_data[0].values
number_of_plots = int(Count / Copies)
##
fig, axes = plt.subplots(number_of_plots, 2, figsize=(15 * 2, 8 * number_of_plots))
curr_ax = 0
melting_points = {}
for i in range(1, Count, Copies):
    plot_data = melting_data[melting_data.columns[i:i + Copies]]
    for j in plot_data.columns:
        curr_data = plot_data[j]
        derivative = -np.gradient(curr_data, 0.5)
        # melting_points[Well_plate[curr_ax]] = melting_points.get(Well_plate[curr_ax], 0) + Temperature[np.argmax(derivative[20:])]
        rfu = sns.lineplot(ax=axes[curr_ax, 0], y=curr_data, x=Temperature)
        der_rfu = sns.lineplot(ax=axes[curr_ax, 1], y=derivative, x=Temperature)
        rfu.set(xlabel='Temperature')
        rfu.set(ylabel=Well_plate[curr_ax])
        der_rfu.set(xlabel='Temperature')
    # der_rfu.set(ylabel=f"Melting Point = {melting_points[Well_plate[curr_ax]] / Copies}")
    curr_ax += 1
plt.show()
# plt.savefig('./figures/result.png')
##
