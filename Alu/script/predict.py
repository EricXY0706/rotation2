import pandas as pd
import numpy as np
import joblib
from utils import extend_to_uniform_length
import pyBigWig
from tqdm import *

model = joblib.load(r'/analysis2/xuy/Alu/models/svm_linear.joblib')
results_all, results_1, results_2,  results_3 = [], [], [], []
for ko in tqdm(['AP2S1', 'CMTR1', 'DLD', 'NELFA', 'NELFB', 'NELFCD']):
# for ko in tqdm(['CMTR1', 'NELFA', 'NELFB', 'NELFCD']):
    for type in ['increase', 'decrease']:
        file = pd.read_table(rf'/analysis2/xuy/RNA-seq/E250005786/6-ExpressionData/RepeatMasker_Unique/{ko}_Alu_{type}.txt')
        for index, row in file.iterrows():
            row = list(row)[0].split('_')
            alu = row[0]
            if alu == 'FLAM':
                alu += row[1]
            start, stop, strand = int(row[-3]), int(row[-2]), row[-1]
            chrom = ''.join(item for item in row if item[:3] == 'chr')
            result, probs = [], []
            for rep in [1, 2, 3]:
                strand_file = {'+': 'plus', '-': 'minus'}
                bw_file = pyBigWig.open(rf'/analysis2/xuy/RNA-seq/E250005786/3-deeptools_CPM_bw_byStrand/{ko}_rep{rep}.{strand_file[strand]}.CPM.bw')
                counts = bw_file.values(chrom, start, stop)
                extended_counts = np.array(extend_to_uniform_length(counts, 3432)).reshape(1, -1)
                y_pre, y_prob = model.predict(extended_counts), model.predict_proba(extended_counts)[:, 1][0]
                if y_pre == 1:
                    result.append([alu, ko, type, chrom, start, stop, strand, rep, round(y_prob, 2)])
                    results_all.append([alu, ko, type, chrom, start, stop, strand, rep, round(y_prob, 2)])
                    probs.append(round(y_prob, 2))
            if len(result) >= 1:
                avg_prob = np.mean(np.array(probs))
                globals()['results_'+str(len(result))].append([alu, ko, type, chrom, start, stop, strand, round(avg_prob, 2)])
                
df = pd.DataFrame(results_all)
df.columns = ['Alu', 'ko', 'increase/decrease', 'chomosome', 'start', 'stop', 'strand', 'rep', 'probability']
df.to_csv(rf'/analysis2/xuy/RNA-seq/E250005786/6-ExpressionData/RepeatMasker_Unique/soloAlu_prediction_all.csv', index=False)
for num in [1, 2, 3]:
	df = pd.DataFrame(globals()['results_'+str(num)])
	df.columns = ['Alu', 'ko', 'increase/decrease', 'chomosome', 'start', 'stop', 'strand', 'probability']
	df.to_csv(rf'/analysis2/xuy/RNA-seq/E250005786/6-ExpressionData/RepeatMasker_Unique/soloAlu_prediction_{num}.csv', index=False)