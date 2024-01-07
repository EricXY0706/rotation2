from sklearn.mixture import GaussianMixture
from sklearn.svm import SVC
from sklearn.metrics import roc_curve, auc, confusion_matrix
import numpy as np
from utils import set_seed
import random
import time
from tqdm import *
import joblib

set_seed(42)
# Load data
pos_1, pos_2 = np.load(r'/analysis2/xuy/Alu/training_data/pos/SRP_9_14_1.npy'), np.load(r'/analysis2/xuy/Alu/training_data/pos/SRP_9_14_2.npy')
neg = np.load(r'/analysis2/xuy/Alu/training_data/neg/neg.npy')
train_test_ratio = 0.8
pos_train_index = random.sample(range(0, pos_1.shape[0]), int(pos_1.shape[0] * train_test_ratio))
pos_test_index = list(set(list(range(0, pos_1.shape[0]))) - set(pos_train_index))
neg_train_index = random.sample(range(0, neg.shape[0]), int(neg.shape[0] * train_test_ratio))
neg_test_index = list(set(list(range(0, neg.shape[0]))) - set(neg_train_index))
x_train = np.concatenate((pos_1[pos_train_index, :], pos_2[pos_train_index, :], neg[neg_train_index, :]), axis=0)
x_test = np.concatenate((pos_1[pos_test_index, :], pos_2[pos_test_index, :], neg[neg_test_index, :]), axis=0)
y_train = np.array([1] * (len(pos_train_index) * 2) + [0] * len(neg_train_index))
y_test = np.array([1] * (len(pos_test_index) * 2) + [0] * len(neg_test_index))
# x_train: (27740, 3432) x_test: (6936, 3432)
# y_train: (27740,) y_test: (6936,)

# Model
start_time = time.time()
svm = SVC(kernel='linear', probability=True)
svm.fit(x_train, y_train)
y_pre, y_prob = svm.predict(x_test), svm.predict_proba(x_test)[:, 1]
fpr, tpr, _ = roc_curve(y_test, y_prob)
AUC = round(auc(fpr, tpr), 4)
confusions = confusion_matrix(y_test, y_pre)
TP, FP, TN, FN = confusions[1, 1], confusions[0, 1], confusions[0, 0], confusions[1, 0]
ACC = round((TP + TN) / (TP + FP + TN + FN), 4)
print(f'AUC: {AUC}, ACC: {ACC}')
print(f'Time elapsed: {round((time.time() - start_time) / 60, 2)} min.\n')
joblib.dump(svm, r'/analysis2/xuy/Alu/models/svm_linear.joblib')

# independent test
idp = np.load(r'/analysis2/xuy/Alu/training_data/idp/nonsoloAlu.npy')
y_idp = svm.predict(idp)
acc = round((np.sum(y_idp) / idp.shape[0]), 4)
print(f'independent test acc: {acc}')