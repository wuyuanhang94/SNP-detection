
# read labels
import numpy as np
file_read = open("phenotype.txt","r")
labels = [];
lines = file_read.readlines()
for line in lines:
    labels.append(line[0])
file_read.close()
output_y = np.array(labels)


# load GenData.mat, aa for 0, Aa for 1, and AA for 2
import scipy.io as sio # read .mat file
mat_contents = sio.loadmat('GenData.mat')
mat_GenData = mat_contents['GenData']
mat_GenData = mat_GenData + 1
[num_samples, num_features] = np.shape(mat_GenData)
output_y = np.array(labels)

# exclude MAF < 0.05
indices_MAF = []
for i in xrange(num_features):
    cnt1 = sum(mat_GenData[:,i] == 0)
    cnt2 = sum(mat_GenData[:,i] == 1)
    MAF = (cnt1 * 2 + cnt2) / (2.0 * num_samples)
    if MAF > 0.5: 
        MAF = 1 - MAF
    if MAF < 0.05:
        indices_MAF.append(i)
if len(indices_MAF) != 0:
    mat_GenData = np.delete(mat_GenData, indices_MAF, 1)
else:
    print('no MAF < 0.5')

# main body
from sklearn.ensemble import RandomForestClassifier
rfc = RandomForestClassifier(n_estimators = 1000, oob_score = True, \
                            verbose = 1000, n_jobs = -1, max_depth = 8)
rate_q = 0.1 # the larger rate_q, the smaller num_iter
rate_r = 0.9
SN_Rank = []
num_snp = num_features
num_iter = np.int(np.floor(np.log2(4.0/num_snp) / np.log2(1 - rate_q))) + 1
print(num_iter)
for i in xrange(num_iter):
    weight_snp = np.array(num_snp * [0.0])
    # calculate the distance
    max_distance_matrix = np.zeros(shape = (num_snp, num_snp), dtype = int)
    max_distance_sample = np.zeros(shape = (num_snp, 2), dtype = int)
    for row in np.arange(0, num_samples / 2): # 0 to 499
        for row_i in np.arange(0, num_samples / 2): # inner for 0 label
            dist = sum(mat_GenData[row] != mat_GenData[row_i])
            if max_distance_matrix[row, 0] < dist:
                max_distance_matrix[row, 0] = dist
                max_distance_sample[row, 0] = row_i
        for row_i in np.arange(num_samples / 2, num_samples): # outer for 0 label
            dist = sum(mat_GenData[row] != mat_GenData[row_i])
            if max_distance_matrix[row, 1] < dist:
                max_distance_matrix[row, 1] = dist
                max_distance_sample[row, 1] = row_i
    for row in np.arange(num_samples / 2, num_samples): # 500 to 999
        for row_i in np.arange(0, num_samples / 2): # outer for 1 label
            dist = sum(mat_GenData[row] != mat_GenData[row_i])
            if max_distance_matrix[row, 1] < dist:
                max_distance_matrix[row, 1] = dist
                max_distance_sample [row, 1]= row_i
        for row_i in np.arange(num_samples / 2, num_samples): # inner for 1 label
            dist = sum(mat_GenData[row] != mat_GenData[row_i])
            if max_distance_matrix[row, 0] < dist:
                max_distance_matrix[row, 0] = dist
                max_distance_sample[row, 0] = row_i
    # done the calculation of distance
    for j in xrange(num_samples):
        for k in xrange(num_snp):            
            weight_snp[k] = weight_snp[k] - (np.float(mat_GenData[j][k] \
            == mat_GenData[max_distance_sample[j][0]][k]) - np.float(mat_GenData[j][k] \
            == mat_GenData[max_distance_sample[j][1]][k])) / np.float(num_snp)
    # weight rank reorder
    indices_sort = np.argsort(weight_snp)[::-1] # inverse order
    num_RF = np.int(np.ceil(rate_r * num_snp))
    indices_sort = indices_sort[-num_RF:]
    
    # RF training
    rfc.fit(mat_GenData[:, indices_sort], output_y)
    print(rfc.oob_score_)
    importances = rfc.feature_importances_
    indices_importance = np.argsort(importances)[::-1]
    num_del = np.int(np.ceil(rate_q * num_snp))
    SN_Rank = SN_Rank + list(indices_importance[-num_del:][::-1])
    np.delete(mat_GenData, indices_importance[-num_del:], 1) # delete column
    num_snp = num_snp - num_del
# the algorithm is terminated
SN_Rank = np.array(SN_Rank[::-1]) + 1
# SN_Rank[:20]

# plot the feature importances for 9445 snps
import matplotlib.pyplot as plt
num_feature = len(SN_Rank)
plt.figure()
plt.title("Distribution of feature importance of Random Forest")
x_axis = np.arange(1,num_feature + 1, dtype = np.int)
plt.plot(x_axis, importances, 'b')
plt.xlabel("Order")
plt.ylabel("Feature importance")
plt.axis([0, 9950])
plt.grid(True)
plt.show()