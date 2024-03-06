#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 22/02/2021
# Author: Sadettin Y. Ugurlu & David McDonald

import os
import pickle
import random
import sys
import warnings
import xlsxwriter
import numpy as np
from autogluon.tabular import TabularDataset, TabularPredictor
import pandas as pd
from sklearn.metrics import average_precision_score
import numpy as np
from sklearn.metrics import f1_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import recall_score
from sklearn.metrics import precision_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import Normalizer,StandardScaler,RobustScaler,MinMaxScaler
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
import sklearn
from sklearn.linear_model import RidgeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import GradientBoostingClassifier
import random
import numpy as np
import mxnet as mx
from autogluon.tabular import TabularPredictor


# -Training proteins:
training=['1AO0', '1RX2', '1DD7', '3IAD', '2CLH', '3AO1', '3IYD', '1V4S', '1I7S', '2YC3', '3MK6', '3IDB', '1Z8D', '3I0R', '2D5Z', '3F3V', '3EPS', '1EFA', '3MKS', '2EWN', '1XJE', '2BU2', '1COZ', '1HAK', '3GR4', '3RZ3', '1EGY', '2ZMF', '1PJ3', '3PTZ', '2XO8', '1SHJ', '1DB1', '1CE8', '1S9J', '1QTI', '2Q5O', '2OI2', '1ESM', '2POC', '2X1L', '1XTU', '2BND', '2I80', '3GCP', '2AL4', '1X88', '3O2M', '3CQD', '3FIG', '3HO8', '1LTH', '1FAP', '3HV8', '3GVU', '3PJG', '3H30', '1T49', '1RD4', '2V92', '2C2B', '3FZY', '3NJQ', '3UO9', '1W96', '2GS7', '3IJG', '1ZDS', '3F6G', '2PUC', '2R1R', '2VGI', '1KP8', '3OS8', '1W25', '3PEE', '3QEL', '1LDN', '1XLS', '1PFK', '3IRH', '1FTA', '2QF7', '3BEO', '3ZLK', '4AVB', '1QW7', '4B9Q', '1TUG', '1PEQ']

# -Filtered protein based on TM-scores:
# There is no protein pairs having higher than 0.5 TM-scores accross bencmarks and training set
test_1=['3QOP', '11BG', '2XJC', '1H9G', '3QH0', '4BO2', '1UXV', '4I1R', '4AW0', '2Q8M', '1NJJ', '3F9N', '2HIM', '1DKU', '1W0F', '1OF6', '3GCD', '2I7N', '2BE9', '3N1V', '3LAJ', '2HVW', '4ETZ', '4HSG', '3O96', '4OO9', '3HNC', '4NBN', '1JLR', '1FX2', '3E5U', '4EBW', '3HO6', '1FIY', '2JFN', '3PYY', '4BZB', '4MBS', '4B1F', '3KGF', '4C7B', '2VPR', '1HKB', '2A69', '4BQH', '2Y0P', '3LW0', '3LU6', '3KF0', '3PXF', '1M8P', '2RD5', '4JAF', '4EO6', '2YLO', '3KCG']
test_2=['3QOP', '11BG', '2XJC', '1H9G', '3QH0', '4BO2', '1UXV', '4I1R', '4AW0', '2Q8M', '1NJJ', '3F9N', '2HIM', '1DKU', '1W0F', '1OF6', '3GCD', '2I7N', '2BE9', '3N1V', '3LAJ', '2HVW', '4ETZ', '4HSG', '3O96', '4OO9', '3HNC', '4NBN', '1JLR', '1FX2', '3E5U', '4EBW', '3HO6', '1FIY', '2JFN', '3PYY', '4BZB', '4MBS', '4B1F', '3KGF', '4C7B', '2VPR', '1HKB', '2A69', '4BQH', '2Y0P', '3LW0', '3LU6', '3KF0', '3PXF', '1M8P', '2RD5', '4JAF', '4EO6', '2YLO', '3KCG']
test_3=['1J07', '1I72', '1YP2', '3BRK', '1ECB', '4CFH', '4EAG', '3KH5', '3L76', '1WQW', '4OP0', '1UWH', '1CKK', '3J41', '3I54', '2OZ6', '2FSZ', '1CSM', '4UUU', '4DQW', '1KMP', '1NV7', '2VD4', '1NE7', '4LZ5', '4U5B', '4PKN', '3E2A', '1VEA', '3TUV', '1L5G', '3BLW', '4MQT', '1O0S', '1S9I', '3D2P', '1FCJ', '1TBF', '2K31', '2PA3', '2H06', '2QMX', '3OF1', '2VK1', '1A3W', '4IP7', '1XMS', '3CMU', '1HK8', '3RSR', '2ONB', '2NW8', '1I6K', '4NES', '2JJX', '3NWY', '1PZO', '4OR2', '4RQZ', '1Q3E', '2PTM', '2VVT', '4Q0A', '4B6E', '3PMA', '3HWS', '1UM8', '2BTY', '1AZX', '1XXA', '3KJN', '1MC0', '2C18', '4LEG', '4TPW', '4DLR', '4QSK', '3UVV', '3HL8', '3LMH', '1OJ9', '3RHW', '3P4W', '2Q6H', '4JKT', '2QXL', '3QKU', '3AUX', '3AV0', '3THO', '4GQQ', '4M0Y', '1T4G', '2Y39', '3DBA', '3K8S', '3ATH', '4H39', '1BM7', '3FUD', '3JPY', '3F1O', '3ZM9', '1BJ4', '4TPT', '3CEV', '3ZFZ', '4FXY', '4M1P', '4NIL', '4OHF', '4CLL', '4PPU', '4Q9M', '3QAK', '4PHU', '3UT3', '2O8B', '2RDE', '4BXC', '4QPL', '3PNA']


#  PASSER gave number for each protein in training set. We used their number to shuffle and select 50 different training set.

index_dictionary = {'0': '1AO0', '1': '1RX2', '2': '1DD7', '3': '3IAD', '4': '2CLH', '5': '3AO1', '6': '3IYD',
                        '7': '1V4S', '8': '1I7S', '9': '2YC3', '10': '3MK6', '11': '3IDB', '12': '1Z8D', '13': '3I0R',
                        '14': '2D5Z', '15': '3F3V', '16': '3EPS', '17': '1EFA', '18': '3MKS', '19': '2EWN',
                        '20': '1XJE', '21': '2BU2', '22': '1COZ', '23': '1HAK', '24': '3GR4', '25': '3RZ3',
                        '26': '1EGY', '27': '2ZMF', '28': '1PJ3', '29': '3PTZ', '30': '2XO8', '31': '1SHJ',
                        '32': '1DB1', '33': '1CE8', '34': '1S9J', '35': '1QTI', '36': '2Q5O', '37': '2OI2',
                        '38': '1ESM', '39': '2POC', '40': '2X1L', '41': '1XTU', '42': '2BND', '43': '2I80',
                        '44': '3GCP', '45': '2AL4', '46': '1X88', '47': '3O2M', '48': '3CQD', '49': '3HO8',
                        '50': '1FAP', '51': '3HV8', '52': '3GVU', '53': '3PJG', '54': '3H30', '55': '1T49',
                        '56': '1RD4', '57': '2V92', '58': '2C2B', '59': '3FZY', '60': '3NJQ', '61': '3UO9',
                        '62': '1W96', '63': '2GS7', '64': '3IJG', '65': '1ZDS', '66': '3F6G', '67': '2PUC',
                        '68': '2R1R', '69': '2VGI', '70': '1KP8', '71': '3OS8', '72': '1W25', '73': '3PEE',
                        '74': '3QEL', '75': '1LDN', '76': '1XLS', '77': '1PFK', '78': '3IRH', '79': '1FTA',
                        '80': '3BEO', '81': '3ZLK', '82': '4AVB', '83': '4B9Q', '84': '1TUG', '85': '1PEQ'}

seed = 42
random.seed(seed)
np.random.seed(seed)
mx.random.seed(seed)
TabularPredictor.random_seed = seed

# turn off warnings
warnings.filterwarnings('ignore')

current_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(current_dir, "data")

# To undersampled based on protein requires splited data.
def split_training_set(filename):
    training_protein_files = f"{current_dir}/data/training_one_by_one"
    if not os.path.exists(training_protein_files):
        os.makedirs(training_protein_files)
    os.chdir(current_dir)
    data = pd.read_csv(f"{filename}.csv",low_memory=False)
    # Split the 'protein_pocket' column by '_' and expand it into two separate columns
    data[['protein_id', 'pocket_number']] = data['protein_pocket'].str.split('_', expand=True)

    os.chdir(training_protein_files)
    protein_list = []
    for key, value in data.iterrows():
        protein_list.append(value["protein_id"])
    protein_list = set(protein_list)

    for target in protein_list:
        split_data = data[data["protein_id"] == str(target)]
        #shuffled_data = split_data.sample(frac=1)
        split_data.to_csv('%s.csv' % target, header=True)

# In research, we split training set to use different training set and compare performance of models on three set.
# In order to represent the pipeline, we will use WHOLE TRAINING SET, and test peroformance of models ONLY ONCE.
# Please find different splits in supporting directory. Once the code is executed, it will produce 51 MODELS USED IN STUDY.
def train_a_model(different_split_dictionary):
    os.chdir(current_dir)
    # protein to train model to make undersampling
    # undersampling should be protein by protein
    for key_split,value_split in different_split_dictionary.items():
        trainIndex=value_split[0]
        testIndex=value_split[1]
        frames=[]
        for protein_index in trainIndex:
            protein = index_dictionary.get(str(protein_index))
            data_temp_train = pd.read_csv(f"{data_dir}/training_one_by_one/{protein}.csv",index_col=1,low_memory=False)
            ##############################################################################################################
            # PASSER provides delete negative based on fpocket order or random. Here, we used random undersampling
            # Random Undersapling:
            # The number of pockets ARE NOT IN ORDER.
            # calculate maximum negative number (it should be 5 times higher than posivite at most)
            positive_number = sum(data_temp_train["Label"]) * 5
            selected_positive = positive_number
            total_pocket_number = len(data_temp_train["Label"])
            drop_number = total_pocket_number - positive_number
            droop_list = []
            for key_1, value in data_temp_train.iterrows():
                if total_pocket_number > selected_positive:
                    if value["Label"] == 0:
                        droop_list.append(key_1)
            # Pockets are not orders, so delete first part or last part is totally in random.
            droop_list = (droop_list[:int(drop_number)])
            # Drop rows
            data_temp_train = data_temp_train.drop(droop_list)
            # take frames in a list
            frames.append(data_temp_train)
        # Concat frames to create training set
        training_data = pd.concat(frames)
        ####################################################################################################################
        # PASSER2.0 only uses Fpocket features
        fpocket_features = ["Label", 'Score', 'Druggability Score', 'Number of Alpha Spheres', "Total SASA", 'Polar SASA', 'Apolar SASA',
                   'Volume', 'Mean local hydrophobic density', 'Mean alpha sphere radius', 'Mean alp. sph. solvent access',
                   'Apolar alpha sphere proportion', 'Hydrophobicity score', 'Volume score', 'Polarity score',
                   'Charge score',
                   "Proportion of polar atoms",
                   'Alpha sphere density', 'Cent. of mass - Alpha Sphere max dist', 'Flexibility',]

        train_data_fpocket = training_data.loc[:, fpocket_features]

        train_data_fpocket = TabularDataset(train_data_fpocket)
        # Once auto_stacking is True, Autogluon automatically select validation data from training set.
        TabularPredictor(label="Label", eval_metric='roc_auc').fit(
            train_data_fpocket, auto_stack=True,
        )
        ####################################################################################################################
        ####################################################################################################################
        # MEF-Allosite uses four models. In order train four models:
        # 1- Gradient Boosting + Boruta features to train the first model
        gradient_boruta = ['Label', 'Druggability Score', 'Number of Alpha Spheres', 'Apolar SASA',
                           'biopython.csvcharge_at_pH', '_SecondaryStrD1100']

        train_data_gra_bor = training_data.loc[:, gradient_boruta]

        train_data_gra_bor = TabularDataset(train_data_gra_bor)
        # Once auto_stacking is True, Autogluon automatically select validation data from training set.
        TabularPredictor(label="Label", eval_metric='roc_auc').fit(
            train_data_gra_bor, auto_stack=True,
        )
        ###################################################################################################################
        # 2- Random Forest + Boruta features to train the second model
        random_boruta = ["Label", 'Score', 'Druggability Score', 'Number of Alpha Spheres', 'Total SASA',
                         'Polar SASA', 'Apolar SASA', 'Volume', 'Mean local hydrophobic density',
                         'Apolar alpha sphere proportion', 'Hydrophobicity score', 'Charge score',
                         'Proportion of polar atoms', 'Alpha sphere density', 'Cent. of mass - Alpha Sphere max dist',
                         'biopython.csvmw', 'biopython.csvcharge_at_pH', 'QSOSW35', 'QSOSW37', 'QSOgrant37',
                         '_SecondaryStrC2', '_PolarityC1']
        train_data_ran_bor = training_data.loc[:, random_boruta]

        train_data_ran_bor = TabularDataset(train_data_ran_bor)
        # Once auto_stacking is True, Autogluon automatically select validation data from training set.
        TabularPredictor(label="Label", eval_metric='roc_auc').fit(
            train_data_ran_bor, auto_stack=True,
        )
        ####################################################################################################################
        # 3- Adaboosting model based feature selection to train the third model
        adaboosting = ["Label",'Score', 'Number of Alpha Spheres', 'Mean alpha sphere radius',
                       'Druggability Score', 'biopython.csvaromaticty',
                       'biopython.csvisoelectric_point']
        train_data_ada = training_data.loc[:, adaboosting]

        train_data_ada = TabularDataset(train_data_ada)
        # Once auto_stacking is True, Autogluon automatically select validation data from training set.
        TabularPredictor(label="Label", eval_metric='roc_auc').fit(
            train_data_ada, auto_stack=True,
        )
        ####################################################################################################################
        # 3- Gradient boosting model based feature selection to train the third model
        gradient_boosting = ["Label",'Number of Alpha Spheres', 'biopython.csvmw', 'Druggability Score',
                             'QSOgrant34', 'Score', 'QSOgrant23', 'Mean local hydrophobic density',
                             'QSOgrant35', 'QSOSW17', 'QSOSW23', 'QSOgrant17', 'Volume', 'QSOSW35',
                             'GG', '_SecondaryStrD1075', 'QR', 'biopython.csvcharge_at_pH',
                             'Flexibility', 'biopython.csvaromaticty', 'Mean alpha sphere radius',
                             'biopython.csvinstability']
        train_data_gra_boo = training_data.loc[:, gradient_boosting]

        train_data_gra_boo = TabularDataset(train_data_gra_boo)
        # Once auto_stacking is True, Autogluon automatically select validation data from training set.
        TabularPredictor(label="Label", eval_metric='roc_auc').fit(
            train_data_gra_boo, auto_stack=True,
        )


def rename_models(model_path):
    model_list = os.listdir(
        model_path)
    os.chdir(model_path)
    # During the training, we trained passer first then our four models. That's why we renamed them
    repeat_number = 0
    # Her tekrar için 5 modeli birlikte işleyin
    for repeat, _ in enumerate(range(0, len(model_list), 5)):
        # Model isimlerini daha anlamlı bir şekilde yeniden adlandırın
        os.rename(sorted(model_list)[repeat * 5], f"{repeat_number}_passer")
        os.rename(sorted(model_list)[repeat * 5 + 1], f"{repeat_number}_gradient_boruta")
        os.rename(sorted(model_list)[repeat * 5 + 2], f"{repeat_number}_random_boruta")
        os.rename(sorted(model_list)[repeat * 5 + 3], f"{repeat_number}_adaboosting")
        os.rename(sorted(model_list)[repeat * 5 + 4], f"{repeat_number}_gradient_boosting")
        # Bir sonraki tekrara geçin
        repeat_number += 1




########################################################################################################################
# split_training_set("training_ready")
different_split_dictionary = {'first' : ([0, 3, 4, 6, 7, 8, 10, 15, 16, 19, 21, 24, 25, 26, 28, 29, 30, 31, 33, 34, 37, 38, 41, 42, 44, 46, 47, 48, 50, 53, 54, 56, 57, 59, 64, 65, 66, 69, 70, 71, 72, 73, 74, 75, 76, 78, 79, 82, 84, 85,43, 68, 58, 1, 80, 67, 22, 23, 5, 45, 39, 9, 52, 36, 63, 13, 12, 62],
                               [81, 60, 35, 18, 55, 51, 17, 20, 11, 32, 2, 49, 27, 14, 77, 83, 40, 61]),
                    'second' : ([0, 1, 2, 3, 5, 6, 7, 12, 14, 15, 17, 21, 23, 27, 28, 29, 31, 32, 33, 35, 38, 40, 42, 43, 46, 48, 50, 51, 52, 53, 54, 55, 58, 60, 61, 62, 63, 64, 65, 66, 71, 74, 75, 76, 78, 81, 82, 83, 84, 85, 47, 4, 20, 18, 77, 26, 19, 80, 39, 68, 30, 16, 73, 36, 9, 8, 45, 56],
                               [49, 13, 79, 41, 37, 24, 67, 72, 69, 44, 57, 34, 59, 10, 25, 11, 70, 22]),
                    'third' : ([1, 3, 4, 5, 6, 7, 10, 13, 14, 17, 18, 20, 23, 25, 26, 28, 30, 33, 36, 40, 42, 43, 44, 45, 47, 49, 50, 52, 54, 56, 57, 58, 60, 63, 65, 67, 68, 69, 70, 71, 72, 74, 75, 76, 77, 78, 79, 81, 82, 85, 19, 59, 21, 64, 2, 84, 24, 27, 8, 51, 48, 53, 0, 41, 34, 38, 37, 83],
                               [9, 73, 61, 31, 66, 35, 11, 32, 29, 12, 15, 39, 22, 62, 46, 80, 16, 55])}
train_a_model(different_split_dictionary)
rename_models(model_path=f"{current_dir}/AutogluonModels")