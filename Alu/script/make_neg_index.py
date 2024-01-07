import pybedtools
import pyBigWig
from tqdm import *
import pandas as pd
file = pybedtools.BedTool(r'/analysis2/xuy/Alu/ALU.strand.bed')
strand_file = {'+': 'fd', '-': 'rv'}
neg = []
for line in file:
    line = list(line)
    chrom, start_base, stop_base, strand = line[0], int(line[1]), int(line[2]), line[5]
    neg.append((chrom, start_base, stop_base, strand))
pos = []
SRP_9_file = open('/analysis2/xuy/Alu/SRP9&14/SRP9_copy.txt').readlines()[:-1]
for line in SRP_9_file:
    line = line.split(',')
    strand, chrom = line[0][-1], line[1].split(':')[0]
    start_base, stop_base = int(line[1].split(':')[1].split('-')[0]), int(line[1].split(':')[1].split('-')[1][:-1])
    pos.append((chrom, start_base-1, stop_base, strand))
SRP_14_file = open('/analysis2/xuy/Alu/SRP9&14/SRP14_copy.txt').readlines()[:-1]
for line in SRP_14_file:
    line = line.split(',')
    strand, chrom = line[0][-1], line[1].split(':')[0]
    start_base, stop_base = int(line[1].split(':')[1].split('-')[0]), int(line[1].split(':')[1].split('-')[1][:-1])
    pos.append((chrom, start_base-1, stop_base, strand))
print(len(neg), len(pos))
neg = [list(item) for item in (set(neg) - set(pos)) if '.' not in list(item)[0]]
neg.sort()
print(len(neg))
df = pd.DataFrame(neg)
df.columns = ['chrom', 'start', 'stop', 'strand']
df.to_csv(r'/analysis2/xuy/Alu/training_data/negative_total.csv', index=False)