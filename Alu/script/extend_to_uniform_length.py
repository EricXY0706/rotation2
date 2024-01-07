import pyBigWig
from collections import Counter
import pandas as pd
import numpy as np
from utils import extend_to_uniform_length

BW1 = pyBigWig.open(r'bw/SRP9_EN_rep1_sorted_cpm_fd.bw')
BW2 = pyBigWig.open(r'bw/SRP9_EN_rep1_sorted_cpm_rv.bw')
file = pd.read_csv(r'intervals.csv')
strand_file = {'+': 'fd', '-': 'rv'}
max_length = np.max(np.array(file['total_length'].tolist()))
for index, row in file.iterrows():
    BW_file = pyBigWig.open(rf'bw/SRP9_EN_rep1_sorted_cpm_{strand_file[row.strand]}.bw')
    counts = BW_file.values(row.chromosome, row.start_base + row.start_interval, row.stop_base + row.stop_interval)
    counts = extend_to_uniform_length(counts, max_length)
    print(counts)