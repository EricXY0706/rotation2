import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import pyBigWig
# srp9_1, srp9_2 = np.load(r'/analysis2/xuy/Alu/training_data/pos/SRP_9_1.npy'), np.load(r'/analysis2/xuy/Alu/training_data/pos/SRP_9_2.npy')
# srp14_1, srp14_2 = np.load(r'/analysis2/xuy/Alu/training_data/pos/SRP_14_1.npy'), np.load(r'/analysis2/xuy/Alu/training_data/pos/SRP_14_2.npy')
# print(srp9_1.shape, srp9_2.shape)
# print(srp14_1.shape, srp14_2.shape)
# srp_9_14_1 = np.concatenate((srp9_1, srp14_1), axis=0)
# srp_9_14_2 = np.concatenate((srp9_2, srp14_2), axis=0)
# print(srp_9_14_1.shape, srp_9_14_2.shape)
# np.save(r'/analysis2/xuy/Alu/training_data/pos/SRP_9_14_1.npy', srp_9_14_1)
# np.save(r'/analysis2/xuy/Alu/training_data/pos/SRP_9_14_2.npy', srp_9_14_2)

# for srp in [9, 14]:
# 	Alu_SRP_total_file = pybedtools.BedTool(rf'/analysis2/xuy/Alu/ALU.SRP{srp}.bed')
# 	Alu_SRP_total = []
# 	for line in Alu_SRP_total_file:
# 		line = list(line)
# 		Alu_SRP_total.append((line[0], line[1], line[2], line[-1]))

# 	solo_alu_SRP_file = open(rf'/analysis2/xuy/Alu/soloAlu_SRP_9_14/SRP{srp}_copy.txt').readlines()[:-1]
# 	solo_alu_SRP = []
# 	for line in solo_alu_SRP_file:
# 		strand = line.split(',')[0][-1]
# 		line = line.split(',')[1]
# 		chrom = line.split(':')[0]
# 		start, stop = line.split(':')[1].split('-')[0], line.split(':')[1].split('-')[1][:-1]
# 		solo_alu_SRP.append((chrom, str(int(start)-1), stop, strand))
# 	solo_alu_SRP.sort()

# 	print(f'SRP{srp}\ntotal: {len(Alu_SRP_total)}, solo: {len(solo_alu_SRP)}')
# 	alu_SRP = list(set(Alu_SRP_total) - set(solo_alu_SRP))
# 	alu_SRP.sort()
# 	print(f'nonsolo: {len(alu_SRP)}\n')
# 	alu_SRP_ = [list(item) for item in alu_SRP]
# 	df = pd.DataFrame(alu_SRP_)
# 	df.columns = ['chrom', 'start', 'stop', 'strand']
# 	df.to_csv(rf'/analysis2/xuy/Alu/nonsoloAlu_SRP_9_14/Alu_nonsolo_SRP{srp}.csv', index=False)
# 	df.to_csv(rf'/analysis2/xuy/Alu/nonsoloAlu_SRP_9_14/Alu_nonsolo_SRP{srp}.bed', sep='\t', header=False, index=False)

# file = pd.read_csv(r'/analysis2/xuy/RNA-seq/E250005786/6-ExpressionData/RepeatMasker_Unique/soloAlu_prediction_1.csv')
# with open(r'/analysis2/xuy/RNA-seq/E250005786/6-ExpressionData/RepeatMasker_Unique/soloAlu_prediction_2.txt', 'w+') as fileobj:
#     for index, row in file.iterrows():
#         if row.ko == 'NELFCD':
#             print(f'{row.Alu}, {row.ko}, {row.chomosome}:{row.start}-{row.stop}, {row.strand}, {row.probability}')
        # fileobj.write(f'{row.Alu}, {row.ko}, {row.chomosome}:{row.start}-{row.stop}, {row.strand}, {row.probability}\n')
a = np.array([1, 2, 3])
print(a / np.sum(a))