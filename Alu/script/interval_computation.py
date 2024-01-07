import pyBigWig
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os
import re

strand_file = {'+': 'fd', '-': 'rv'}
SRP = open(r'SRP9&14/SRP14_copy.txt').readlines()[:-1]
results = []
for line in SRP:
    strand = line.split(',')[0][-1]
    pattern = '|'.join(map(re.escape, ['_', ',', ':', '-']))
    line = re.split(pattern, line)
    chrom, start_base, stop_base = line[-3], int(line[-2]), int(line[-1][:-1])
    bw1 = pyBigWig.open(rf'bw/SRP14_EN_rep1_sorted_cpm_{strand_file[strand]}.bw')
    bw2 = pyBigWig.open(rf'bw/SRP14_EN_rep2_sorted_cpm_{strand_file[strand]}.bw')
    start_base_Alu, stop_base_Alu = start_base, stop_base
    bw = bw1
    if bw.values(chrom, start_base_Alu, stop_base_Alu)[0] == 0.0 and bw.values(chrom, start_base_Alu, stop_base_Alu)[-1] == 0.0:
        results.append([strand, chrom, start_base, start_base_Alu, stop_base, stop_base_Alu])
    elif strand == '+':
        while bw.values(chrom, start_base_Alu, stop_base_Alu)[-1] != 0.0:
            stop_base_Alu += 1
        results.append([strand, chrom, start_base, start_base_Alu, stop_base, stop_base_Alu])
    elif strand == '-':
        while bw.values(chrom, start_base_Alu, stop_base_Alu)[0] != 0.0:
            start_base_Alu -= 1
        results.append([strand, chrom, start_base, start_base_Alu, stop_base, stop_base_Alu])
interval = []
for result in results:
    if result[2] == result[3] and result[5] == result[4]:
        interval.append([result[0], result[1], result[2], result[4], result[2] - result[3], result[5] - result[4]])
    elif result[2] != result[3] and result[5] == result[4]:
        interval.append([result[0], result[1], result[2], result[4], result[3] - result[2] + 1, result[5] - result[4]])
    elif result[5] != result[4] and result[2] == result[3]:
        interval.append([result[0], result[1], result[2], result[4], result[2] - result[3], result[5] - result[4] - 1])
df = pd.DataFrame(interval)
df.columns = ['strand', 'chromosome', 'start_base', 'stop_base', 'start_interval', 'stop_interval']
lengths = []
for index, row in df.iterrows():
    lengths.append(row.stop_base + row.stop_interval - (row.start_base + row.start_interval))
df['total_length'] = pd.Series(lengths)
df.to_csv(r'/analysis2/xuy/SRP9&14_intervals/SRP9_rep1_intervals.csv', index=False)

# BW1 = pyBigWig.open(r'bw/SRP9_EN_rep1_sorted_cpm_fd.bw')
# BW2 = pyBigWig.open(r'bw/SRP9_EN_rep1_sorted_cpm_rv.bw')
# print('--------------------------------------------------------------------')
# print(BW1.values('chr14', 54959482, 54959796))
# print('--------------------------------------------------------------------')
# print(BW2.values('chr7', 23387936, 23388236))
    # plt.figure()
    # os.system("export DISPLAY=:0.0")
    # plt.switch_backend('WebAgg')
    # plt.bar(x=list(range(start_base-100, stop_base+100)), height=bw1_values_background, color='blue', label='background')
    # plt.bar(x=list(range(start_base, stop_base)), height=bw1_values, color='red', label='Alu')
    # plt.legend()
    # plt.show()