import sys
import pickle
import numpy as np

#********************************************************************#
# load model
loaded_model = pickle.load(open(str('OM3_training_clean.model'), 'rb'))

# load features file
csvfile = open(sys.argv[1],'r')

# Classifying protein
for key in csvfile:
    
    X = np.asarray(key.split(',')[1::])
    p = loaded_model.predict_proba(X)
    pred = loaded_model.predict(X)
    print (X[0] + '\t' + pred + '\t' + p + '\n')