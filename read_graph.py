from scipy.sparse import coo_matrix
from scipy.io     import savemat


filename = "cit-HepPh.txt"

data = []
row = []
col = []
with open(filename) as f:
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
savemat('test.mat', {'A' : sparseMat});
