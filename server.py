import os
import flask
from flask import Flask
from flask_cors import CORS

import scanpy as sc
import numpy as np
import pandas as pd
#import dropkick as dk
import scipy.stats

# create Flask app
app = Flask(__name__)
CORS(app)

# make-em-up data
my_data = ['Hi Client. I am Server response '+str(idx+1)+'. Nice to meet you.' for idx in range(100)]
adata = np.load("AS1_4258-raw.npy")
dkscore = np.load("AS1_4258-dkscore.npy")
sorted_ind = np.argsort(dkscore)
top_sort = sorted_ind[-50:]
top_adata = adata[top_sort,:]
print(top_adata.shape)
print(adata.shape)
print(dkscore.shape)
sample = ['Hi Victoria, this is the Pearson Correlation Coefficient: ']
#droplet = adata[0,:]

@app.route('/grab_data', methods=['GET','POST'])
def get_signals():
    # our way of grabbing data from the client
    client_data = flask.request.json
    new_item = client_data['new_item']
    new_data = my_data[new_item]

    # this is returned to the client - note, it will take data_obj, and convert it to a Javascript object
    data_obj = {'val':new_data}
    return flask.jsonify(data_obj)
#
@app.route('/get_genes', methods=['GET','POST'])
def get_genes():
    # get genes from the client
    gene_data = adata[0]
    gene_data = gene_data.tolist()

    # returned to the client
    return flask.jsonify(gene_data)

@app.route('/get_dkscore', methods=['GET','POST'])
def get_dkscore():
    # get dropkick score from the client
    score = dkscore.tolist()

    return flask.jsonify(score)

@app.route('/get_pearson_coef', methods=['GET','POST'])
def get_pearson_coef():
    # get pearson coef and convert for client
    coef_data = flask.request.json
    new_item = coef_data['new_item']
    new_data = sample[new_item]

    # returned to client
    coef_obj = {'val': new_data}
    return flask.jsonify(coef_obj)

@app.route('/calc_pearson_coef', methods=['GET','POST'])
def calc_pearson_coef():
    # get minScore and maxScore data for calculating pearson coef
    scoreRange = flask.request.json
    minScore = scoreRange[0]
    maxScore = scoreRange[1]
    
    # calculate pearson correlation coefficient between two droplets
    score_ind = np.where((dkscore >= minScore) & (dkscore <= maxScore))
    score_adata = adata[score_ind,:]
    coef_list = []
    for i in range(score_ind[0].shape[0] - 1):
        for j in range(score_ind[0].shape[0] - 1):
            coef_num = scipy.stats.pearsonr(score_adata[:,i][0], score_adata[:,j][0])[0]
            coef = {'droplet1': int(score_ind[0][i]), 'droplet2': int(score_ind[0][j]), 'coef_num': coef_num}
            coef_list.append(coef)

    return flask.jsonify(coef_list)
    

# execute the application (by default, it should be hosted at localhost:5000, which you will see in the output)
if __name__ == '__main__':
    app.run()
