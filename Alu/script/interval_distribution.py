import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from collections import Counter
fig, ax = plt.subplots()

intervals_9_1 = pd.read_csv(r'/analysis2/xuy/Alu/soloAlu_SRP_9_14_intervals/SRP9_rep1_intervals.csv')
file = intervals_9_1
fd_length = file[file['strand'] == '+']['total_length'].tolist()
rv_length = file[file['strand'] == '-']['total_length'].tolist()
v1 = ax.violinplot(fd_length, positions=[1], vert=True, showmeans=True, showmedians=True)
v2 = ax.violinplot(rv_length, positions=[2], vert=True, showmeans=True, showmedians=True)

intervals_9_2 = pd.read_csv(r'/analysis2/xuy/Alu/soloAlu_SRP_9_14_intervals/SRP9_rep2_intervals.csv')
file = intervals_9_2
fd_length = file[file['strand'] == '+']['total_length'].tolist()
rv_length = file[file['strand'] == '-']['total_length'].tolist()
v3 = ax.violinplot(fd_length, positions=[3], vert=True, showmeans=True, showmedians=True)
v4 = ax.violinplot(rv_length, positions=[4], vert=True, showmeans=True, showmedians=True)

intervals_14_1 = pd.read_csv(r'/analysis2/xuy/Alu/soloAlu_SRP_9_14_intervals/SRP14_rep1_intervals.csv')
file = intervals_14_1
fd_length = file[file['strand'] == '+']['total_length'].tolist()
rv_length = file[file['strand'] == '-']['total_length'].tolist()
v5 = ax.violinplot(fd_length, positions=[5], vert=True, showmeans=True, showmedians=True)
v6 = ax.violinplot(rv_length, positions=[6], vert=True, showmeans=True, showmedians=True)

intervals_14_2 = pd.read_csv(r'/analysis2/xuy/Alu/soloAlu_SRP_9_14_intervals/SRP14_rep2_intervals.csv')
file = intervals_14_2
fd_length = file[file['strand'] == '+']['total_length'].tolist()
rv_length = file[file['strand'] == '-']['total_length'].tolist()
v7 = ax.violinplot(fd_length, positions=[7], vert=True, showmeans=True, showmedians=True)
v8 = ax.violinplot(rv_length, positions=[8], vert=True, showmeans=True, showmedians=True)

labels = ['SRP9_rep1_forward', 'SRP9_rep1_reverse', 'SRP9_rep2_forward', 'SRP9_rep2_reverse', 'SRP14_rep1_forward', 'SRP14_rep1_reverse',
          'SRP14_rep2_forward', 'SRP14_rep2_reverse']
x = range(len(labels) + 1)
plt.xticks(x[1:], labels, fontsize=6, rotation=15)
plt.savefig(r'/analysis2/xuy/Alu/soloAlu_SRP_9_14_intervals/interval_distribution.pdf', dpi=300, bbox_inches='tight')