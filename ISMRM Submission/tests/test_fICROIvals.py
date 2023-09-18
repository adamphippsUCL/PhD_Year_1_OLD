# Python script to experiment with extracting fIC values from ROIs

# Import relevant libraries
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.io import loadmat


# Specify patient number
PatNum = 'INN_129'

# === Load fIC, fEES, and fVASC results

# Original 
fIC_original = loadmat(f'VERDICT outputs/{PatNum}/Original/fIC.mat')['fIC']
fIC_original = np.moveaxis( fIC_original , -1, 0)

fEES_original = loadmat(f'VERDICT outputs/{PatNum}/Original/fEES.mat')['fEES']
fEES_original = np.moveaxis( fEES_original , -1, 0)

fVASC_original = loadmat(f'VERDICT outputs/{PatNum}/Original/fVASC.mat')['fVASC']
fVASC_original = np.moveaxis( fVASC_original , -1, 0)


# No VASC 
fIC_NoVASC = loadmat(f'VERDICT outputs/{PatNum}/No VASC/fIC.mat')['fIC']
fIC_NoVASC = np.moveaxis( fIC_NoVASC , -1, 0)

fEES_NoVASC = loadmat(f'VERDICT outputs/{PatNum}/No VASC/fEES.mat')['fEES']
fEES_NoVASC = np.moveaxis( fEES_NoVASC , -1, 0)

fVASC_NoVASC = loadmat(f'VERDICT outputs/{PatNum}/No VASC/fVASC.mat')['fVASC']
fVASC_NoVASC = np.moveaxis( fVASC_NoVASC , -1, 0)


# Reduced Rs
fIC_ReducedRs = loadmat(f'VERDICT outputs/{PatNum}/Reduced Rs/fIC.mat')['fIC']
fIC_ReducedRs = np.moveaxis( fIC_ReducedRs , -1, 0)

fEES_ReducedRs = loadmat(f'VERDICT outputs/{PatNum}/Reduced Rs/fEES.mat')['fEES']
fEES_ReducedRs = np.moveaxis( fEES_ReducedRs , -1, 0)

fVASC_ReducedRs = loadmat(f'VERDICT outputs/{PatNum}/Reduced Rs/fVASC.mat')['fVASC']
fVASC_ReducedRs = np.moveaxis( fVASC_ReducedRs , -1, 0)

# === Load ROI mask

LesionMask = np.load(f'ROIs/{PatNum}/lesion.npy')
LesionMask = (LesionMask != 0)


# == Extract fIC, fEES, and fVASC values
Lesion_fICs_original = fIC_original[LesionMask]
Lesion_fEESs_original = fEES_original[LesionMask]
Lesion_fVASCs_original = fVASC_original[LesionMask]

Lesion_fICs_NoVASC = fIC_NoVASC[LesionMask]
Lesion_fEESs_NoVASC = fEES_NoVASC[LesionMask]
Lesion_fVASCs_NoVASC = fVASC_NoVASC[LesionMask]

Lesion_fICs_ReducedRs = fIC_ReducedRs[LesionMask]
Lesion_fEESs_ReducedRs = fEES_ReducedRs[LesionMask]
Lesion_fVASCs_ReducedRs = fVASC_ReducedRs[LesionMask]

Lesion_original_ratio = Lesion_fICs_original/(Lesion_fICs_original+Lesion_fEESs_original)


plt.figure()
plt.plot(Lesion_fICs_original)
plt.plot(Lesion_fICs_ReducedRs)
plt.plot(Lesion_fICs_NoVASC)
plt.show()