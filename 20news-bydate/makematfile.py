from scipy.sparse import coo_matrix
from scipy.io     import savemat
import numpy

value = []
row = []
col = []

with open('./matlab/train.data') as f:
    for line in f:
        triple = line.split();
        value.append(float(triple[2]))
        row.append(int(triple[0]) - 1)
        col.append(int(triple[1]) - 1)

trainMat = coo_matrix((value,(row,col)));
a = trainMat.tocsr()

value = []
row = []
col = []
with open('./matlab/test.data') as f:
    for line in f:
       triple = line.split();
       value.append(float(triple[2]))
       row.append(int(triple[0]) - 1)  
       col.append(int(triple[1]) - 1)


testMat = coo_matrix((value,(row,col))); 
labelsList = []

with open('./matlab/train.label') as f:
    for line in f:
        labelsList.append(float(line))

train_labels = numpy.array(labelsList)


labelsList = []
with open('./matlab/test.label') as f:
    for line in f:
        labelsList.append(float(line))

test_labels = numpy.array(labelsList)

sparseMat = coo_matrix((value,(row,col)))
savemat('../matlab/20newsgroup.mat', {'Problem' : {'name': '20newsgroup', 'train_data' : trainMat,  'train_labels': train_labels, 'test_data' : testMat, 'test_labels':  test_labels}})
