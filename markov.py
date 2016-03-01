# -*- coding: utf-8 -*-
"""
Created on Mon Feb 29 21:16:29 2016
Changing the frequency of a markov chain
@author: michaud
"""

import numpy as np

# input matrix
Ptwo = np.ones((2,2))
Ptwo[0,0] = 0.8
Ptwo[0,1] = 0.2
Ptwo[1,0] = 0.0
Ptwo[1,1] = 1.0
print "two-year transition = \n",Ptwo

# get eigen values and vectors
w,v = np.linalg.eig(Ptwo)

# take root of eigen values
sw = np.diag(np.sqrt(w))

# reconstruct 
Pone = np.dot(np.dot(v, sw) , np.linalg.inv(v))

print "one-year transition = \n",Pone

# check 

Ptwo_check = np.dot(Pone,Pone)

print "check for two transitions \n", Ptwo_check