from scipy.sparse import coo_matrix
from scipy.io     import savemat
import os
graphFiles = [x for x in os.listdir('./matlab/snap_mats/') if x[-4:] == '.txt']
for filename in graphFiles:

    data = []
    row = []
    col = []
    with open('./matlab/snap_mats/'+filename) as f:
        print(filename)
        for line in f:
            if line[0] == '#':
                continue
            pair = line.split();
            data.append(1.0);
            row.append(int(pair[0]));
            col.append(int(pair[1]));

    nodes = row+col
    nodes = list(dict.fromkeys(nodes))
    nod2idx = {nodes[i] : i for i in range(len(nodes))}
    for i in range(len(row)):
        row[i] = nod2idx[row[i]]
        col[i] = nod2idx[col[i]]

    
    sparseMat = coo_matrix((data,(row,col)),(len(nodes),len(nodes)))
    savemat('./matlab/snap_mats/'+filename[:-3] + 'mat', {'Problem' : {'A' : sparseMat, 'name' : filename[:-4]}});
