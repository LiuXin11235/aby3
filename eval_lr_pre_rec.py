import os
import argparse
import subprocess
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import Counter as cnt
import logging
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, f1_score, average_precision_score, precision_score, recall_score
import sys
import pickle
import yaml
from multiprocessing import Pool


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--party_n", type=str, help="number of data parties", default="3")
    args = parser.parse_args()  
    
    # read the prediction output
    prediction_output = pd.read_csv('prediction_output.csv')
    y_true = prediction_output['true_label'].values
    y_pred = prediction_output['pred_label'].values
    # read the pruned labels
    for i in range(int(args.party_n)):
        try:
            pruned_label = pd.read_csv('lr_train/pruned_test_label_party_%d.csv'%(i), header=None)
            y_true = np.concatenate([y_true, pruned_label.values.flatten()])
            y_pred = np.concatenate([y_pred, np.zeros(pruned_label.shape[0])])
        except FileNotFoundError as e:
            continue

    pre_score = precision_score(y_true=y_true, y_pred=y_pred)
    rec_score = recall_score(y_true=y_true, y_pred=y_pred)
    f1_score = f1_score(y_true=y_true, y_pred=y_pred)
    print('Precision: %.4f, Recall: %.4f, F1: %.4f'%(pre_score, rec_score, f1_score))
    print('See the prediction output in prediction_output.csv')  
