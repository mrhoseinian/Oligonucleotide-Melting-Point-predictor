##
# Libraries
# import math
import random

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

##
# Adjusting the model
K = 10
# Initial_concentration = 1
Number_of_Reads = 1000000
Mean_yeild = 0.002
Sd_yeild = 0.0005
Subset_fraction = 0.002
Edit_parameter = 2 - Subset_fraction
Subset_size = Subset_fraction * 4 ** K
Initial_count = 10000
Wash_elution = 1
Number_of_extra_molecules = 1
Interference_fraction = Subset_fraction * Number_of_extra_molecules


# Noise = 0.3


##
# Generate Library
def generate_library():
    # library = np.zeros(4 ** K)
    # # count = 0
    # for i in range(4 ** K):
    #     yield_fraction = random.gauss(Mean_yeild, Sd_yeild)
    #     if yield_fraction < 0:
    #         yield_fraction = 0
    #     if yield_fraction > 1:
    #         yield_fraction = 1
    #     library[i] = yield_fraction
    # return library
    library = np.ones(4 ** K) * Initial_count
    return library


##
# simulate selection
def selection(library):
    library_size = 4 ** K
    random_indices = random.sample(population=range(library_size), k=round(Subset_size * Edit_parameter))
    dna_counts = np.zeros(4 ** K, dtype=int)
    selected_counts = np.zeros(len(random_indices), dtype=int)
    random_indices.sort()
    total = 0
    for i, j in enumerate(random_indices):
        yield_fraction = random.gauss(Mean_yeild, Sd_yeild)
        if yield_fraction < 0:
            yield_fraction = 0
        if yield_fraction > 1:
            yield_fraction = 1
        after_selection_count = round(yield_fraction * library[j])
        dna_counts[j] = after_selection_count * Wash_elution
        selected_counts[i] = after_selection_count * Wash_elution
        total += dna_counts[j]

    return random_indices, dna_counts, total, selected_counts


# def selection():
#     library_size = 4 ** K
#     random_indices = random.sample(population=range(library_size), k=round(Subset_size * Edit_parameter))
#     random_indices.sort()
#     dna_counts = np.zeros(len(random_indices), dtype=int)
#     total = 0
#     for i in range(len(random_indices)):
#         yield_fraction = random.gauss(Mean_yeild, Sd_yeild)
#         if yield_fraction < 0:
#             yield_fraction = 0
#         if yield_fraction > 1:
#             yield_fraction = 1
#         after_selection_count = round(yield_fraction * Initial_count)
#         dna_counts[i] = after_selection_count * Wash_elution
#         total += dna_counts[i]
#     return random_indices, dna_counts, total# def selection():


##
def negative_selection(library):
    library_size = 4 ** K
    random_indices = random.sample(population=range(library_size), k=round(Subset_size * Edit_parameter))
    dna_counts = np.copy(library)
    # selected_counts = np.zeros(len(random_indices), dtype=int)
    random_indices.sort()
    total = 0
    for j in random_indices:
        yield_fraction = random.gauss(Mean_yeild, Sd_yeild)
        if yield_fraction < 0:
            yield_fraction = 0
        if yield_fraction > 1:
            yield_fraction = 1
        after_selection_count = round(yield_fraction * library[j])
        dna_counts[j] -= after_selection_count * Wash_elution
        total += dna_counts[j]
    return random_indices, dna_counts, total


##
#
def interference(random_indices):
    interference_size = round((Interference_fraction ** 2) * (4 ** K) * (Edit_parameter ** 2))
    interference_indices = random.sample(population=random_indices, k=interference_size)
    # interference_DNAs = np.zeros(len(interference_indices), dtype=int)
    interference_DNAs = np.zeros(4 ** K, dtype=int)
    # interference_indices.sort()
    for j in interference_indices:
        yield_fraction = random.gauss(Mean_yeild, Sd_yeild)
        interference_count = round(yield_fraction * Initial_count)
        interference_DNAs[j] = interference_count

    return interference_DNAs, interference_indices


##
# Readout
def readout(selection_counts, selection_indices):
    # X = 0
    # selection_size = len(list_count)
    # list_readout = np.zeros(len(selection_counts), dtype=int)  # list_read_inter =
    list_readout = np.zeros(4 ** K, dtype=int)
    list_noise = np.zeros(4 ** K, dtype=int)
    total_noise = 0
    total_read = 0
    noise_weights = np.ones(4 ** K, dtype=int)
    for i in selection_indices:
        noise_weights[i] = 0

    read_dna_indices = random.choices(population=selection_indices, weights=selection_counts,
                                      k=round(Number_of_Reads * (1 - Noise)))
    read_noise_indices = random.choices(population=range(4 ** K), weights=noise_weights,
                                        k=round(Number_of_Reads * Noise))
    # read_dna_indices.sort()
    # read_noise_indices.sort()
    for i in read_dna_indices:
        # if list_count[j] != 0:
        # if list_inter[i] == 0:
        total_read += 1
        # if i in dict_readout:
        # dict_readout[i] += 1
        list_readout[i] += 1
    # else:
    #     dict_readout[i] = 1

    for i in read_noise_indices:
        total_noise += 1
        # if i in dict_noise:
        #     dict_noise[i] += 1
        # else:
        #     dict_noise[i] = 1
        list_noise[i] += 1
    return list_readout, list_noise, total_read, total_noise  # , list_inter


##
# Plotting
def rawplot(list_read, list_noise, list_inter):
    plt.figure(num=None, figsize=(25, 20), dpi=150, facecolor='w', edgecolor='k')
    x = range(4 ** K)
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    plt.scatter(x, list_read, c='b', label='signal')
    plt.scatter(x, list_noise, c='r', label='noise')
    # plt.scatter(x, list_inter, c='y', label='interference')
    plt.legend(loc='upper left')
    plt.xlabel('DNA Library')
    plt.ylabel('Frequency')
    ymin, ymax = plt.ylim()
    plt.ylim([0, ymax + 50])
    plt.savefig(
        'figures/Noise=' + str(Noise) + '   Mean_yeild=' + str(Mean_yeild) + '   Interference_fraction=' + str(
            Interference_fraction) + '.png')
    plt.close()


##
# Plotting for negative
def negative_lib_plotter(list_read):
    plt.figure(num=None, figsize=(60, 30), dpi=150, facecolor='w', edgecolor='k')
    x = range(4 ** K)
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    plt.scatter(x, list_read, c='b', s=0.05, label='signal')
    # plt.scatter(x, list_noise, c='r', label='noise')
    # plt.scatter(x, list_inter, c='y', label='interference')
    plt.legend(loc='upper left')
    plt.xlabel('DNA Library')
    plt.ylabel('Frequency')
    ymin, ymax = plt.ylim()
    plt.ylim([0, ymax + 50])
    plt.savefig(
        'negative_lib_figures/Noise=' + str(Noise) + '   Mean_yeild=' + str(Mean_yeild) + '   Interference_fraction='
        + str(Interference_fraction) + '.png')
    plt.close()


##
# Plotting
def negative_rawplot(list_read, list_noise, selection_indices):
    plt.figure(num=None, figsize=(25, 20), dpi=150, facecolor='w', edgecolor='k')
    x = range(4 ** K)
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    plt.scatter(x, list_read, c='b', label='signal')
    plt.scatter(x, list_noise, c='r', label='noise')
    for i in selection_indices:
        plt.axvline(i, linewidth=0.01)
    plt.legend(loc='upper left')
    plt.xlabel('DNA Library')
    plt.ylabel('Frequency')
    ymin, ymax = plt.ylim()
    plt.ylim([0, ymax + 100])
    plt.savefig(
        'negative_figures/Noise=' + str(Noise) + '   Mean_yeild=' + str(Mean_yeild) + '   Interference_fraction=' + str(
            Interference_fraction) + '.png')
    plt.close()


##
# data
def data_extractor(list_read, list_noise, list_inter, selection_indices):
    selected_sequences = np.count_nonzero(list_read)
    noisy_sequences = np.count_nonzero(list_noise)
    interfered_sequences = np.count_nonzero(list_inter)
    # selection_average = np.average(list_read[selection_indices])
    selection_average = np.average(list_read)
    noise_average = np.average(list_noise)
    # interference_average = np.average(list_inter[interference_indices])
    with open("data/statistics.txt", "a") as file:
        file.write("For Noise =%f and Mean_yeild=%f and interference_fraction = %f: \n" % (
            Noise, Mean_yeild, Interference_fraction))
        file.write("# read sequences = %d and average = %d  and max = %d\n" % (selected_sequences, selection_average,
                                                                               np.max(list_read)))
        file.write("# noise sequences = %d and average = %f  and max = %d\n" % (noisy_sequences, noise_average,
                                                                                np.max(list_noise)))
        file.write("# non-interfered sequences = %d \n" % (selected_sequences - interfered_sequences))
        file.write("Average SNR = %f \n\n\n" % (selection_average / noise_average))


##
def histo_plotter(list_read, filename):
    list_copy = list_read[np.where(list_read != 0)]
    mu, std = norm.fit(list_copy)
    plt.hist(list_copy, bins=20, density=True, color='g')
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
    plt.plot(x, p, 'k', linewidth=2)
    title = "mu = %.2f,  std = %.2f" % (mu, std)
    plt.title(title)
    right_lim = mu + 3 * std
    left_lim = mu - 3 * std
    plt.axvline(right_lim, linewidth=0.1, c='r')
    plt.axvline(left_lim, linewidth=0.1, c='r')
    plt.savefig('histo/histo_' + filename + '.png')
    plt.close()
    # if filename == 'neg':
    with open("histo/histostat_"+ filename+ ".txt", "a") as file:
        file.write(
            "number of reads less than mu-3sd=%d\n" % (
                        np.count_nonzero(list_read < left_lim) - np.count_nonzero(list_read == 0)))
        file.write("number of reads greater than mu+3sd=%d" % (np.count_nonzero(list_read > right_lim)))


##
# Tester
Noise = 0.1
Library = generate_library()
Selection_indices, Selection_Counts, Selected_total, Selected_summary = selection(Library)
Interference_counts, Interference_indices = interference(Selection_indices)
List_read, List_noise, Total_read, Total_noise = readout(Selected_summary, Selection_indices)

##
'''Noise'''
Library = generate_library()
Selection_indices, Selection_Counts, Selected_total, Selected_summary = selection(Library)
Interference_counts, Interference_indices = interference(Selection_indices)
Noise_fraction = [0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99]
for Noise in Noise_fraction:
    List_read, List_noise, Total_read, Total_noise = readout(Selected_summary, Selection_indices)
    rawplot(List_read, List_noise, Interference_counts)
    data_extractor(List_read, List_noise, Interference_counts, Selection_indices)

##
Noise = 0.5
''''Mean Yield'''
Mean_yeilds = [0.001, 0.0005, 0.0001]
for Mean_yeild in Mean_yeilds:
    Library = generate_library()
    Selection_indices, Selection_Counts, Selected_total, Selected_summary = selection(Library)
    Interference_counts, Interference_indices = interference(Selection_indices)
    List_read, List_noise, Total_read, Total_noise = readout(Selected_summary, Selection_indices)
    rawplot(List_read, List_noise, Interference_counts)
    data_extractor(List_read, List_noise, Interference_counts, Selection_indices)

##
'''Interference'''
Noise = 0.5
Mean_yeild = 0.002
Interference_fractions = [0.002, 0.004, 0.006, 0.01]
Library = generate_library()
Selection_indices, Selection_Counts, Selected_total, Selected_summary = selection(Library)
for Interference_fraction in Interference_fractions:
    Interference_counts, Interference_indices = interference(Selection_indices)
    List_read, List_noise, Total_read, Total_noise = readout(Selected_summary, Selection_indices)
    rawplot(List_read, List_noise, Interference_counts)
    data_extractor(List_read, List_noise, Interference_counts, Selection_indices)

##
''' Negative part'''
# plotting the positive

Noise = 0.1
Mean_yeild = 0.5
Sd_yeild = 0.15
Subset_fraction = 0.02
Number_of_Reads = 10000000
Edit_parameter = 2 - Subset_fraction
Subset_size = Subset_fraction * 4 ** K
Interference_fraction = 0.002
Library = generate_library()
Selection_indices, Selection_Counts, Selected_total, Selected_summary = selection(Library)
Interference_counts, Interference_indices = interference(Selection_indices)
List_read, List_noise, Total_read, Total_noise = readout(Selected_summary, Selection_indices)
histo_plotter(List_read, 'pos')

rawplot(List_read, List_noise, Interference_counts)
# data_extractor(List_read, List_noise, Interference_counts, Selection_indices)

# not A => B #
N_Selection_indices, N_Selection_Counts, N_Selected_total = negative_selection(Library)
negative_lib_plotter(N_Selection_Counts)
Selection_indices, Selection_Counts, Selected_total, Selected_summary = selection(N_Selection_Counts)
Interference_counts, Interference_indices = interference(Selection_indices)
List_read, List_noise, Total_read, Total_noise = readout(Selected_summary, Selection_indices)
histo_plotter(List_read, 'neg')
negative_rawplot(List_read, List_noise, N_Selection_indices)
# data_extractor(List_read, List_noise, Interference_counts, Selection_indices)


##
