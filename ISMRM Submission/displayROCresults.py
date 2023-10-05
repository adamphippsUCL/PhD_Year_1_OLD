# Python script to display ROC results

# Import relevant libraries
import matplotlib.pyplot as plt
import numpy as np
import glob
import os


# Find model numbers which have been analysed
ModelNumbers = [ int(os.path.split(os.path.split(folder)[0])[1][-1]) for folder in glob.glob('ROC results/Model */') ]

# Initialise figure
fig = plt.figure()
ax = fig.add_subplot(111)

for ModelNum in ModelNumbers:
    
    # Read in ROC results
    fpr = np.load(f'ROC results/Model {ModelNum}/fpr.npy')
    tpr = np.load(f'ROC results/Model {ModelNum}/tpr.npy')
    
    ax.plot(fpr,tpr,'-*')
    
    
plt.show()


