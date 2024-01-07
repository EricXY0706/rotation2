import pandas as pd
import numpy as np
import random
from matplotlib import pyplot as plt
from utils import set_seed, extend_to_uniform_length
import pyBigWig
from scipy import stats

set_seed(42)
pos_1, pos_2 = np.load(r'/analysis2/xuy/Alu/training_data/pos/SRP_9_14_1.npy'), np.load(r'/analysis2/xuy/Alu/training_data/pos/SRP_9_14_2.npy')
file_1, file_2 = pd.read_csv(r'/analysis2/xuy/Alu/soloAlu_SRP_9_14_intervals/SRP9_14_rep1_intervals.csv'), pd.read_csv(r'/analysis2/xuy/Alu/soloAlu_SRP_9_14_intervals/SRP9_14_rep2_intervals.csv')
train_test_ratio, SRP_9_size, updown_size = 0.8, 4282, 3000
pos_train_index = random.sample(range(0, pos_1.shape[0]), int(pos_1.shape[0] * train_test_ratio))
strand_file = {'+': 'fd', '-': 'rv'}
train_representations = np.zeros((updown_size * 2 + pos_1.shape[1]))
for index in pos_train_index:
	for rep in [1, 2]:
		chrom, strand = globals()['file_'+str(rep)].loc[index, 'chromosome'], globals()['file_'+str(rep)].loc[index, 'strand']
		start = globals()['file_'+str(rep)].loc[index, 'start_base'] + globals()['file_'+str(rep)].loc[index, 'start_interval']
		stop = globals()['file_'+str(rep)].loc[index, 'stop_base'] + globals()['file_'+str(rep)].loc[index, 'stop_interval']
		bw_file_route = ''
		if 0 < index < SRP_9_size:
			bw_file_route = bw_file_route.join(rf'/analysis2/xuy/Alu/bw/SRP9_EN_rep{rep}_sorted_cpm_{strand_file[strand]}.bw')
		else:
			bw_file_route = bw_file_route.join(rf'/analysis2/xuy/Alu/bw/SRP14_EN_rep{rep}_sorted_cpm_{strand_file[strand]}.bw')
		bw_file = pyBigWig.open(bw_file_route)
		upstream = bw_file.values(chrom, start-updown_size, start)
		downstream = bw_file.values(chrom, stop, stop+updown_size)
		center = globals()['pos_'+str(rep)][index, :].tolist()
		counts = upstream + center + downstream
	counts = np.array(counts) / 2
	train_representations += counts
train_representations /= len(pos_train_index)
train_representations /= np.sum(train_representations)

prediction_results = pd.read_csv(r'/analysis2/xuy/RNA-seq/E250005786/6-ExpressionData/RepeatMasker_Unique/soloAlu_prediction_all.csv')
baselines = [0.5, 0.6, 0.7, 0.8, 0.9]
strand_file = {'+': 'plus', '-': 'minus'}
for baseline in baselines:
    prediction_result = prediction_results[(baseline <= prediction_results['probability']) & (prediction_results['probability'] < (baseline + 0.1))]
    prediction_representations = np.zeros((updown_size * 2 + pos_1.shape[1]))
    for index, row in prediction_result.iterrows():
        ko, chrom, start, stop, strand, rep = row.ko, row.chomosome, row.start, row.stop, row.strand, row.rep
        bw_file = pyBigWig.open(rf'/analysis2/xuy/RNA-seq/E250005786/3-deeptools_CPM_bw_byStrand/{ko}_rep{rep}.{strand_file[strand]}.CPM.bw')
        extended_counts = extend_to_uniform_length(bw_file.values(chrom, start, stop), pos_1.shape[1])
        upstream = bw_file.values(chrom, start-updown_size, start)
        downstream = bw_file.values(chrom, stop, stop+updown_size)
        counts = upstream + extended_counts + downstream
        prediction_representations += np.array(counts)
    prediction_representations /= prediction_result.shape[0]
    prediction_representations /= np.sum(prediction_representations)
    KL = stats.entropy(train_representations, prediction_representations)
    plt.figure()
    plt.plot(list(range(len(train_representations))), train_representations, color='red', label='Training set')
    plt.plot(list(range(len(prediction_representations))), prediction_representations, color='blue', label='Prediction')
    plt.xticks([0, updown_size / 2, updown_size, updown_size + pos_1.shape[1] / 2, updown_size + pos_1.shape[1], updown_size * 1.5 + pos_1.shape[1], updown_size * 2 + pos_1.shape[1]],
			[0, 'upstream', updown_size, 'extended counts', updown_size + pos_1.shape[1], 'downstream', updown_size * 2 + pos_1.shape[1]], fontsize=8)
    plt.title(f'threshold: {round(baseline, 1)}~{round(baseline+0.1, 1)} (KLD={round(KL, 2)})')
    plt.legend()
    plt.savefig(rf'/analysis2/xuy/RNA-seq/E250005786/6-ExpressionData/RepeatMasker_Unique/soloAlu_prediction_all_{round(baseline, 1)}_{round(baseline+0.1, 1)}.pdf', dpi=300, bbox_inches='tight')