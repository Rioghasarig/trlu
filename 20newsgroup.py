from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import TfidfVectorizer
from scipy.sparse import csr_matrix
from scipy.io     import savemat
import numpy

#cats = ['alt.atheism', 'comp.graphics']
newsgroups_train = fetch_20newsgroups(subset='train', remove=('headers','footers','quotes'))
newsgroups_test = fetch_20newsgroups(subset='test', remove=('headers','footers','quotes'))
vectorizer = TfidfVectorizer(stop_words='english', max_df=.5, min_df=2,max_features=10000)

train_vectors = vectorizer.fit_transform(newsgroups_train.data)
test_vectors = vectorizer.transform(newsgroups_test.data)
vocab = vectorizer.get_feature_names()

trainMat = csr_matrix(train_vectors,dtype=float)
train_labels = newsgroups_train.target

testMat = csr_matrix(test_vectors,dtype=float)
test_labels = newsgroups_test.target

savemat('./matlab/20newsgroup-py.mat', {'Problem' : {'name': '20newsgroup', 'train_data' : trainMat, 'train_labels': train_labels, 'test_data': testMat, 'test_labels': test_labels}})
with open('./matlab/vocab-py.txt', 'w') as vocab_file:
    for word in vocab:
        vocab_file.write(word+'\n')


