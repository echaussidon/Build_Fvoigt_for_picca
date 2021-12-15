#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Edmond CHAUSSIDON

import glob
import fitsio
import numpy as np
import matplotlib.pyplot as plt

"""

Need to use the delta files from picca of the considered simulation or data release

The output is used in build_Fvoigt.py to compute a weighted average of f(n,r)

We compute the weight for each lambda --> it is the same weight used in the lyman-alpha autocorrelation

"""


files = glob.glob("delta*.gz")

lamb = np.arange(3600, 5501, 1)
weight = np.zeros(lamb.size)

i = 0
nb_files = len(files)
print(" ")
print(i, " /" + str(nb_files))

for file in files:
	data = fitsio.FITS(file)
	for hdu in data[1:]:
		weight += np.interp(lamb, 10**hdu["LOGLAM"][:], hdu["WEIGHT"][:], left=0, right=0)
	i += 1
	print(i, " /" +str(nb_files))

weight_norm = weight/np.sum(weight)

data = np.transpose(np.concatenate((lamb, weight, weight_norm)).reshape(3,lamb.size))
np.savetxt('weight_lambda.txt', data)

print("    ** Calcul terminee :) ")

#plt.plot(lamb, weight)
#plt.show()
