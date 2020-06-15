# -*- coding: utf-8 -*-
# read labels
file_read = open("multi_phenos.txt","r")
labels = [];
lines = file_read.readlines()
for line in lines:
    labels.append(line[:-1].split(' '))
file_read.close()

# load GenData.mat, AA for 2, aa for 0, and Aa for 1
import scipy.io as sio
import numpy as np
mat_contents = sio.loadmat('GenData.mat')
mat_GenData = mat_contents['GenData']
mat_GenData = mat_GenData + 1
output_y = np.array(labels, dtype = np.int)

# training
import matplotlib.pyplot as plt
from sklearn.cross_validation import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import ElasticNet
alpha = 0.001
names = ["Random Forest", "Elastic Net"]
clf1 = RandomForestClassifier(n_estimators = 1000, oob_score = True, \
    max_depth = 8, n_jobs = -1)
clf2 = ElasticNet(alpha=alpha, l1_ratio=0.2)
X_train, X_test, y_train, y_test = train_test_split(mat_GenData, output_y, test_size=0)
num_features = np.shape(X_train)[1]
clf1.fit(X_train, y_train)
clf2.fit(X_train, y_train)
rf_importances = clf1.feature_importances_
elnt_weights = np.abs(clf2.coef_.sum(axis=0))
elnt_idx = elnt_weights.argsort()[::-1]
# plot signal-noise ratio
plt.figure()
plt.title("Normalized importances of RF VS. ElasticNet (Multi_phenos)")
x_axis = np.arange(1,num_features + 1, dtype = np.int)
rf_handle, = plt.plot(x_axis, rf_importances / np.sum(rf_importances), 'r*', label = 'Random Forest')
elnt_handle, = plt.plot(x_axis, elnt_weights / np.sum(elnt_weights), 'b+', label = 'Elastic Net')
plt.legend(handles = [rf_handle, elnt_handle])
plt.xlabel("Order")
plt.ylabel("Feature importance")
plt.axis([0, num_features, np.min(rf_importances / np.sum(rf_importances)), np.max(elnt_weights / np.sum(elnt_weights))])
plt.grid(True)
plt.show()

# plot most important weights for Elastic Net
plt.figure()
plt.title("Ordered weights of Elastic Net")
elnt_handle, = plt.plot(x_axis, elnt_weights[elnt_idx], 'b-', label = 'Elastic Net')
plt.legend(handles = [elnt_handle])
plt.xlabel("Order")
plt.ylabel("Weights")
plt.axis([0, num_features, np.min(elnt_weights), np.max(elnt_weights)])
plt.grid(True)
plt.show()

# re-training RF
num_RF = np.int(num_features * 0.1)
rf_idx_training = elnt_idx[:num_RF] 
clf1.fit(X_train[:,rf_idx_training], y_train)
rf_importances = clf1.feature_importances_
rf_idx = rf_importances.argsort()[::-1]
len_rf = np.int(0.1 * len(rf_importances))
print(rf_idx_training[:len_rf])

# plot importances
plt.figure()
rf_handle, = plt.plot(rf_idx_training[:len_rf], rf_importances[rf_idx[:len_rf]], 'r*', label = 'Random Forest + Elastic Net')
plt.legend(handles = [rf_handle])
plt.title("Chosen SNPs of RF + Elastic Net")
plt.xlabel("Order")
plt.ylabel("Feature Importances")
plt.axis([0, num_features, np.min(rf_importances[rf_idx[:len_rf]]), np.max(rf_importances[rf_idx[:len_rf]])])
plt.grid(True)
plt.show()
