import math
import random

##
# Part 1
# This script takes a list of 'counts' from sequencing data and calculates theoretical #
# equilibrium association constants for each library member. Assumes 3 rds of Selection. #

# required values
total_reads = 190363  # total number of usable sequence reads
diversity = 7077888.0  # total number of distinct molecules in the library
del_conc_initial = 1.7 * 10 ** (-5)  # initial concentration of DEL with desired small-molecule moiety
protein_conc = float(.16 * 10 ** (-6))  # concentration active protein
del_conc_final = 24 * (10 ** (-12))  # Estimated assuming a DNA concentration of 1-10 pM ##prior to PRC amplification
conc_ligand_initial = del_conc_initial / diversity  # starting concentration of each ligand


def getConc(Content, Total_reads, Del_conc_final):
    """takes list of sequence counts and returns list of
    final concentrations for each library member"""
    list1 = []
    for n in Content:
        counts = float(n)
        concitem = (counts / Total_reads) * Del_conc_final
        list1.append(concitem)
    return list1


def getKa(B, A, List_conc_final):
    """takes floats A, B, and list of library member final concentrations
    and returns list of calculated association equilibrium constants"""
    list1 = []
    for i in List_conc_final:
        a = (B ** 3) * i - (B ** 3) * A
        b = 3 * (B ** 2) * i
        c = 3 * B * i
        d = i
        delta = 18 * a * b * c * d - 4 * (b ** 3) * d + (b ** 2) * c ** 2 - 4 * a * c ** 3 - 27 * (a ** 2) * (d ** 2)
        delta0 = b * b - 3 * a * c
        delta1 = 2 * (b ** 3) - 9 * a * b * c + 27 * (a ** 2) * d
        delta2 = -27 * (a ** 2) * delta
        u1 = 1
        # u2 = (-1 + 1.0j*(3**(1.0/2)))/2
        # u3 = (-1 - 1.0j*(3**(1.0/2)))/2
        C = ((delta1 + (delta2 ** (1.0 / 2))) / 2) ** (1.0 / 3)
        K1 = (-1 / (3 * a)) * (b + u1 * C + delta0 / (u1 * C))
        # K2 = (-1/(3*a))*(b + u2*C + delta0/(u2*C))
        # K3 = (-1/(3*a))*(b + u3*C + delta0/(u3*C))
        list1.append(K1)
    return list1


def outputFile(list1, name):
    """takes list and str filename and writes list to file"""
    File1 = open(name, 'w')
    for X in list1:
        File1.write(str(X) + "\n")
    File1.close()


def readFile(name):
    """takes filename as string and reads the file row by row to list. Returns list"""
    with open(name) as f:
        Content = f.readlines()
    f.close()
    return Content


content = readFile('listofrawcounts.txt')
list_conc_final = getConc(content, total_reads, del_conc_final)
list_Ka = getKa(protein_conc, conc_ligand_initial, list_conc_final)
outputFile(list_Ka, 'Kvalues.txt')

##
# part 2
# Python script for sampling of association constant values with initial library concentrations,
# final library concentrations, and protein concentrations treated as random variables

read_count = 205  # sequence reads for a particular library member (ligand)
total_reads = 190363  # total number of usable sequence reads
diversity = 7077888.0  # total number of distinct molecules in the library

list_Ka_range = []
count = 0
while count < 100000:
    # initial concentration of DEL.
    del_conc_initial = random.gauss(1.7 * 10 ** (-5), 0.5 * 10 ** (-5))
    # concentration of active protein
    protein_conc = random.gauss(0.16 * 10 ** (-6), 0.05 * 10 ** (-6))
    # see section S3 for explanation of del_conc_final
    del_conc_final = random.gauss(33.0 * (10 ** (-12)), 30.0 * (10 ** (-12)))
    if protein_conc > 0 and del_conc_initial > 0 and del_conc_final > 0:
        conc_ligand_initial = del_conc_initial / diversity
        list_conc_final = getConc([read_count], total_reads, del_conc_final)
        list_Ka = getKa(protein_conc, conc_ligand_initial, list_conc_final)
        if list_Ka[0] > 0:
            list_Ka_range.append(list_Ka[0])
            count += 1

    outputFile(list_Ka_range, 'sampling_Kvalues.txt')

##
# part 3
# Python script generating test_library
input_kvalues = readFile('Ka_values_from_Clark_etal.txt')
log_values = []
for x in input_kvalues:
    log_values.append(math.log10(float(x)))

count = 0
while count < 7070788 - len(input_kvalues):
    value = abs(random.gauss(0, 2))
    if value < 5.8:
        log_values.append(value)
    count += 1
log_values = sorted(log_values, reverse=True)

file1 = open('example_Ka_values.txt', 'w')

for x in log_values:
    file1.write(str(x) + "\n")

file1.close()

##
# part 4
# takes test_library generated in section S5 and simulates n cycles of selection and sequencing.

# The values below can be adjusted to mimic various selection conditions
volume = 100.0 * 10 ** (-6)  # total volume of the selection
diversity = 707788.0  # total # of unique sequences (structures) in library
del_conc_initial = 5 * 10 ** (-5)  # conc of DEL during incubation with protein
protein_conc = 5.0 * 10 ** (-7)  # conc of protein used during a selection
reads = 190363  # sequence reads
mean_yield = 45  # avg yield individual library member as percentage
yield_sd = 15  # standard deviation of mean_yield

ligand_conc = del_conc_initial / diversity
umolesbegin = float((del_conc_initial * volume) * 1000000)


def readFile(name):
    """takes name and reads the file row by row to list. Returns list"""
    with open(name) as f:
        Content = f.readlines()
    f.close()
    return Content


def assignStartConc(List_gen_Kavalues, mean, sd, Ligand_conc):
    """Takes a list of Ka values, and 3 floats (mean, sd, ligand_conc).
    Generates a random 'yield' form a normal distribution
    and assigns a starting ligand concentration for each library member.
    Returns a list of concentrations (float).
    Also returns total micromoles of DNA before selection."""
    total = 0.0
    List_concentrations = []
    for X in range(0, len(List_gen_Kavalues)):
        concAB = random.gauss(mean, sd) * 0.01 * Ligand_conc
        if concAB < 0:
            concAB = 0
        List_concentrations.append(concAB)
        total = total + concAB
    Umolesend = (total * volume) * 1000000.0
    return List_concentrations, Umolesend


def simulate_n_RoundsSelection(List_concentrations, List_gen_Kavalues, n):
    """takes two lists and an integer. Simulates n cycles of selections
    Returns list of concentrations for each library member (floats).
    Also returns total micromoles of DNA after selection is complete."""
    Umolesend = 0
    for X in range(0, n):
        total = 0.0
        for i in range(0, len(List_gen_Kavalues)):
            ka = 10.0 ** (float(List_gen_Kavalues[i]))
            List_concentrations[i] = (ka * protein_conc * List_concentrations[i] / (1 + ka * protein_conc))
            total = total + List_concentrations[i] * 0.5
        Umolesend = (total * volume) * 1000000.0
    return List_concentrations, Umolesend


def calcEnrichment(List_concentrations, Ligand_conc, Del_conc_initial):
    """takes list of floats and returns list of floats. Converts concentrations
    to enrichment ratio."""
    list_enrichment = []
    for X in range(0, len(List_concentrations)):
        list_enrichment.append((List_concentrations[X] / Ligand_conc) / (Ligand_conc / Del_conc_initial))
    return list_enrichment


def convertEnrichmentToInteger(list_enrichment, Mult_factor):
    """takes list of floats and a float. Enrichment values are adjusted to prevent
    int boxes_total from exceeding 10**8. Returns list and int"""
    Boxes_total = 10 ** 9
    List_boxsize = []
    lessthan1 = 0.0

    while Boxes_total > 10 ** 8:
        List_boxsize = []
        lessthan1 = 0.0
        Mult_factor *= 2
        Boxes_total = 0
        for X in list_enrichment:
            X = float(X / Mult_factor)
            if X >= 1:
                Boxes_total += int(X)
                List_boxsize.append(int(X))
                lessthan1 += X - int(X)
            if X < 1:
                List_boxsize.append(0)
                lessthan1 += X

    List_boxsize.append(int(lessthan1))
    Boxes_total += int(lessthan1)
    return List_boxsize, Boxes_total


def randomChoiceBoxes(Boxes_total):
    """takes int. Returns list with length equal to boxes_total.
    Indices of the list are populated randomly to simulate
    sequencing experiment."""
    List_pickboxes = []
    for X in range(0, Boxes_total):
        List_pickboxes.append(0)

    for X in range(0, reads):
        r = random.randrange(0, Boxes_total)
        List_pickboxes[r] = List_pickboxes[r] + 1
    return List_pickboxes


def assignRandomChoicesToLigands(List_boxsize, List_pickboxes, Mult_factor):
    """takes two lists and a float. Indices populated in list_pickboxes are
    assigned to particular library members according to the information stored
    in list_boxsize. Counts and enrichment for each library member is stored
    in listofcount and list_enrich. Both lists are returned. """

    total = 0
    item = 0
    Listofcounts = []
    List_enrich = []
    for X in List_boxsize:
        total = total + X
        Listofcounts.append(0)
        for i in range(total - X, total):
            Listofcounts[item] = Listofcounts[item] + List_pickboxes[i]
        List_enrich.append(X * Mult_factor)  # enrichment
        item += 1
    return List_enrich, Listofcounts


def printToFile(numberlines, List_enrich, Listofcounts, List_gen_Kavalues):
    """takes an int and 3 lists. Writes to a file information from the 3 lists.
    Writes the first number of lines dictated by the int numberline"""
    if numberlines > len(List_enrich) - 1:
        numberlines = len(List_enrich) - 1
    File1 = open('counts_ryields.txt', 'w')
    for X in range(0, numberlines - 1):
        File1.write(str(List_enrich[X]) + " " + str(Listofcounts[X]) + " " + str(X + 1) + " " +
                    str(List_gen_Kavalues[X]) + "\n")
    File1.close()


list_gen_Kavalues = readFile("example_kavalues.txt")
list_concentrations, umolesend = assignStartConc(list_gen_Kavalues, mean_yield, yield_sd, ligand_conc)
print(umolesend)
list_concentrations, umolesend = simulate_n_RoundsSelection(list_concentrations, list_gen_Kavalues, 3)
list_enrichment_adjusted = calcEnrichment(list_concentrations, ligand_conc, del_conc_initial)

# sequencing
mult_factor = 0.5
list_boxsize, boxes_total = convertEnrichmentToInteger(list_enrichment_adjusted, mult_factor)
list_pickboxes = randomChoiceBoxes(boxes_total)
list_enrich, listofcounts = assignRandomChoicesToLigands(list_boxsize, list_pickboxes, mult_factor)
printToFile(40000, list_enrich, listofcounts, list_gen_Kavalues)
