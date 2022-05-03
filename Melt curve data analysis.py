##
import collections
import string
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()

##
Copies = 6
# Concentrations = 1
Pair = 6
alphabet = list(string.ascii_uppercase)
alphabet = alphabet[:(Pair * Copies) // 12]
# alphabet = alphabet[:(Pair*2 * Copies) // 12]
Well_plate = []
for i in alphabet:
    # for j in range(1, 12 // Copies + 1):
    for j in range(1, (12*2) // Copies + 1):
        Well_plate.append(i + str(j))
##
melting_data = pd.read_csv('./data/Melt Curve RFU Results.csv', header=None)
melting_data = melting_data.drop(columns=melting_data.columns[0])
melting_data = melting_data.drop(labels=0, axis=0)
melting_data.dropna(axis=0, inplace=True)
print(melting_data.head())

##
Count = len(melting_data.columns) - 1
Temperature = melting_data[1].values
number_of_plots = int(Count / Copies)
##
curr_ax = 0
melting_points = collections.defaultdict(list)
# melting_points2 = collections.defaultdict(list)
for i in range(1, Count, 2 * Copies):
    fig, axes = plt.subplots(2, 2, figsize=(15 * 2, 15))
    plot_data = melting_data[melting_data.columns[i:i + Copies]]
    plot_data2 = melting_data[melting_data.columns[i + Copies:i + 2 * Copies]]
    for j, k in zip(plot_data.columns, plot_data2.columns):
        curr_data, curr_data2 = plot_data[j], plot_data2[k]
        derivative, derivative2 = -np.gradient(curr_data, 0.5), -np.gradient(curr_data2, 0.5)
        melting_points[Well_plate[curr_ax]].append(Temperature[20 + np.argmax(derivative[20:])])
        melting_points[Well_plate[curr_ax+1]].append(Temperature[20 + np.argmax(derivative2[20:])])
        # melting_points2[Well_plate[curr_ax]].append(Temperature[20 + np.argmax(derivative2[20:])])
        rfu = sns.lineplot(ax=axes[0, 0], y=curr_data, x=Temperature)
        der_rfu = sns.lineplot(ax=axes[0, 1], y=derivative, x=Temperature)
        rfu2 = sns.lineplot(ax=axes[1, 0], y=curr_data2, x=Temperature)
        der_rfu2 = sns.lineplot(ax=axes[1, 1], y=derivative2, x=Temperature)
    rfu.set(xlabel='Temperature')
    rfu.set(ylabel=Well_plate[curr_ax])
    der_rfu.set(xlabel='Temperature')
    der_rfu.set(
        ylabel=f"Melting Point = {np.mean(melting_points[Well_plate[curr_ax]]):.2f} and std = {np.std(melting_points[Well_plate[curr_ax]]):.2f}")
    rfu2.set(xlabel='Temperature')
    rfu2.set(ylabel=Well_plate[curr_ax+1])
    der_rfu2.set(xlabel='Temperature')
    der_rfu2.set(
        ylabel=f"Melting Point = {np.mean(melting_points[Well_plate[curr_ax+1]]):.2f} and std = {np.std(melting_points[Well_plate[curr_ax+1]]):.2f}")
    plt.savefig('./figures/'+Well_plate[curr_ax]+'.png')
    curr_ax += 2
    plt.show()

##
