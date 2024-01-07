import pyBigWig
import pybedtools
import pandas as pd
import numpy as np
import random
import os
from tqdm import *

def set_seed(seed=42):
    random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    np.random.seed(seed)
def extend_to_uniform_length(original_array, uniform_length):
    original_length = len(original_array)
    multiplier = uniform_length // original_length
    extended_array = []
    for element in original_array:
        extended_array.extend([element] * multiplier)
    extended_array = [extended_array[0]] * int(((uniform_length - len(extended_array)) / 2)) + extended_array
    extended_array.extend([extended_array[-1]] * (uniform_length - len(extended_array)))
    return extended_array
def make_positive_dataset(uniform_length=3432):
    strand_file = {'+': 'fd', '-': 'rv'}
    for SRP in [9, 14]:
        for rep in [1, 2]:
            pos = []
            locus_file = pd.read_csv(rf'/analysis2/xuy/Alu/SRP9&14_intervals/SRP{SRP}_rep{rep}_intervals.csv')
            for index, row in locus_file.iterrows():
                strand, chrom = row.strand, row.chromosome
                start, stop = row.start_base + row.start_interval, row.stop_base + row.stop_interval
                bw_file = pyBigWig.open(rf'/analysis2/xuy/Alu/bw/SRP{SRP}_EN_rep{rep}_sorted_cpm_{strand_file[strand]}.bw')
                counts = bw_file.values(chrom, start, stop)
                extended_counts = extend_to_uniform_length(counts, uniform_length)
                pos.append(extended_counts)
            pos = np.array(pos)
            np.save(rf'/analysis2/xuy/Alu/training_data/pos/SRP_{SRP}_{rep}.npy', pos)
def make_negative_dataset(uniform_length=3432):
    positive_size = np.load(r'/analysis2/xuy/Alu/training_data/pos/SRP_9_14_1.npy').shape[0] * 2
    neg_total = pd.read_csv(r'/analysis2/xuy/Alu/negative_total.csv')
    strand_file = {'+': 'fd', '-': 'rv'}
    neg = []
    for index, row in tqdm(neg_total.iterrows()):
        chrom, start_base, stop_base, strand = row.chrom, row.start, row.stop, row.strand
        for SRP in [9, 14]:
            for rep in [1, 2]:
                bw_file = pyBigWig.open(rf'/analysis2/xuy/Alu/bw/SRP{SRP}_EN_rep{rep}_sorted_cpm_{strand_file[strand]}.bw')
                counts = bw_file.values(chrom, start_base, stop_base)
                if counts != [0.0] * len(counts):
                    extended_counts = extend_to_uniform_length(counts, uniform_length)
                    neg.append(extended_counts)
    neg = np.array(neg)
    neg_index = random.sample(range(0, neg.shape[0]), positive_size)
    neg_use = neg[neg_index, :]
    print(neg.shape, neg_use.shape)
    np.save(r'/analysis2/xuy/Alu/training_data/neg/neg.npy', neg_use)
def make_independent_dataset(uniform_length=3432):
    strand_file = {'+': 'fd', '-': 'rv'}
    idp = []
    for srp in [9, 14]:
        for rep in [1, 2]:
            srp9 = pd.read_csv(rf'/analysis2/xuy/Alu/nonsoloAlu_SRP_9_14/Alu_nonsolo_SRP{srp}.csv')
            for index, row in srp9.iterrows():
                chrom, start, stop, strand = row.chrom, int(row.start), int(row.stop), row.strand
                bw_file = pyBigWig.open(rf'/analysis2/xuy/Alu/bw/SRP{srp}_EN_rep{rep}_sorted_cpm_{strand_file[strand]}.bw')
                counts = bw_file.values(chrom, start, stop)
                if counts != [0.0] * len(counts):
                    extended_counts = extend_to_uniform_length(counts, uniform_length)
                    idp.append(extended_counts)
    idp = np.array(idp)
    np.save(r'/analysis2/xuy/Alu/training_data/idp/nonsoloAlu.npy', idp)
    print(idp.shape)

        
if __name__ == '__main__':
	set_seed(42)
	# make_positive_dataset(3432)
	# make_negative_dataset(3432)
	make_independent_dataset(3432)