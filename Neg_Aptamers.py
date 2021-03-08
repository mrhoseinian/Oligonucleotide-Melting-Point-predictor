##
import Aptamers
import random

import matplotlib.pyplot as plt
import numpy as np

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





