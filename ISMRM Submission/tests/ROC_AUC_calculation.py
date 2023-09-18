# Python script to test ROC_AUC curves

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
import sys

indxs = np.arange(0,1000,1)#range(1000)

# Generate data (Weights of people)
Weights = np.random.uniform(50,150, 1000)


# Healthy
Healthy_score = 1.0*(Weights < 100)

# Add noise to Weights
Measured_Weights = Weights + np.random.uniform(-25,15,1000)



plt.figure()
plt.scatter(indxs[Healthy_score == 1], Measured_Weights[Healthy_score == 1], c = 'tab:green')
plt.scatter(indxs[Healthy_score != 1], Measured_Weights[Healthy_score != 1], c = 'tab:red')



print(1 - roc_auc_score(Healthy_score, Measured_Weights))



# Generate ROC curve
thresholds = np.arange(50,150,2)
YoudenIndxs = []

plt.figure()

# Set threshold
for threshold in thresholds:

    # Caclulate true positive rate
    TP_rate = np.sum( 1.0*(Measured_Weights < threshold)*(Healthy_score == 1) )/np.sum( 1.0*(Measured_Weights < threshold) )

    # Calculate false positive rate (1- true negative rate)
    FP_rate = 1 - np.sum( 1.0*(Measured_Weights > threshold)*(Healthy_score != 1) )/np.sum( 1.0*(Measured_Weights > threshold) )

    plt.scatter(FP_rate, TP_rate)

    # Calculate Youden Index
    J = TP_rate - FP_rate
    YoudenIndxs.append(J)
    
plt.ylim(0,1)
plt.xlim(0,1)
    
    
     
plt.figure()
plt.scatter(thresholds, YoudenIndxs)

plt.show()







