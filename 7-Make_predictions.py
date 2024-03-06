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
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, VotingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
import sklearn
from sklearn.linear_model import RidgeClassifier, SGDClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.svm import SVC,LinearSVC
from sklearn.dummy import DummyClassifier
from sklearn.metrics import roc_auc_score
from sklearn.metrics import recall_score, precision_score,accuracy_score
from tabulate import tabulate
import statistics

# turn off warnings
warnings.filterwarnings('ignore')


# -Filtered protein based on TM-scores:
# There is no protein pairs having higher than 0.5 TM-scores accross bencmarks and training set
test_1=['3QOP', '11BG', '2XJC', '1H9G', '3QH0', '4BO2', '1UXV', '4I1R', '4AW0', '2Q8M', '1NJJ', '3F9N', '2HIM', '1DKU', '1W0F', '1OF6', '3GCD', '2I7N', '2BE9', '3N1V', '3LAJ', '2HVW', '4ETZ', '4HSG', '3O96', '4OO9', '3HNC', '4NBN', '1JLR', '1FX2', '3E5U', '4EBW', '3HO6', '1FIY', '2JFN', '3PYY', '4BZB', '4MBS', '4B1F', '3KGF', '4C7B', '2VPR', '1HKB', '2A69', '4BQH', '2Y0P', '3LW0', '3LU6', '3KF0', '3PXF', '1M8P', '2RD5', '4JAF', '4EO6', '2YLO', '3KCG']
test_2=['3QOP', '11BG', '2XJC', '1H9G', '3QH0', '4BO2', '1UXV', '4I1R', '4AW0', '2Q8M', '1NJJ', '3F9N', '2HIM', '1DKU', '1W0F', '1OF6', '3GCD', '2I7N', '2BE9', '3N1V', '3LAJ', '2HVW', '4ETZ', '4HSG', '3O96', '4OO9', '3HNC', '4NBN', '1JLR', '1FX2', '3E5U', '4EBW', '3HO6', '1FIY', '2JFN', '3PYY', '4BZB', '4MBS', '4B1F', '3KGF', '4C7B', '2VPR', '1HKB', '2A69', '4BQH', '2Y0P', '3LW0', '3LU6', '3KF0', '3PXF', '1M8P', '2RD5', '4JAF', '4EO6', '2YLO', '3KCG']
test_3=['1J07', '1I72', '1YP2', '3BRK', '1ECB', '4CFH', '4EAG', '3KH5', '3L76', '1WQW', '4OP0', '1UWH', '1CKK', '3J41', '3I54', '2OZ6', '2FSZ', '1CSM', '4UUU', '4DQW', '1KMP', '1NV7', '2VD4', '1NE7', '4LZ5', '4U5B', '4PKN', '3E2A', '1VEA', '3TUV', '1L5G', '3BLW', '4MQT', '1O0S', '1S9I', '3D2P', '1FCJ', '1TBF', '2K31', '2PA3', '2H06', '2QMX', '3OF1', '2VK1', '1A3W', '4IP7', '1XMS', '3CMU', '1HK8', '3RSR', '2ONB', '2NW8', '1I6K', '4NES', '2JJX', '3NWY', '1PZO', '4OR2', '4RQZ', '1Q3E', '2PTM', '2VVT', '4Q0A', '4B6E', '3PMA', '3HWS', '1UM8', '2BTY', '1AZX', '1XXA', '3KJN', '1MC0', '2C18', '4LEG', '4TPW', '4DLR', '4QSK', '3UVV', '3HL8', '3LMH', '1OJ9', '3RHW', '3P4W', '2Q6H', '4JKT', '2QXL', '3QKU', '3AUX', '3AV0', '3THO', '4GQQ', '4M0Y', '1T4G', '2Y39', '3DBA', '3K8S', '3ATH', '4H39', '1BM7', '3FUD', '3JPY', '3F1O', '3ZM9', '1BJ4', '4TPT', '3CEV', '3ZFZ', '4FXY', '4M1P', '4NIL', '4OHF', '4CLL', '4PPU', '4Q9M', '3QAK', '4PHU', '3UT3', '2O8B', '2RDE', '4BXC', '4QPL', '3PNA']

current_dir = os.path.dirname(os.path.abspath(__file__))

def model_fn(model_dir):
    """loads model from previously saved artifact"""
    model = TabularPredictor.load(model_dir)  # or model = MultiModalPredictor.load(model_dir) for example
    model.persist_models()  # This line only works for TabularPredictor
    return model

def split_training_set(filename,test_number):
    training_protein_files = f"{current_dir}/data/test_one_by_one/test_{test_number}"
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


def validation_models(model_path, test_list, test_number, repeat_number=3, model_selection="MEF-AlloSite"):
    fpocket_features = ["Label", 'Score', 'Druggability Score', 'Number of Alpha Spheres', "Total SASA",
                        'Polar SASA', 'Apolar SASA',
                        'Volume', 'Mean local hydrophobic density', 'Mean alpha sphere radius',
                        'Mean alp. sph. solvent access',
                        'Apolar alpha sphere proportion', 'Hydrophobicity score', 'Volume score', 'Polarity score',
                        'Charge score',
                        "Proportion of polar atoms",
                        'Alpha sphere density', 'Cent. of mass - Alpha Sphere max dist', 'Flexibility', ]
    gradient_boruta = ["Label", 'Druggability Score', 'Number of Alpha Spheres', 'Apolar SASA',
                       'biopython.csvcharge_at_pH', '_SecondaryStrD1100']
    random_boruta = ["Label", 'Score', 'Druggability Score', 'Number of Alpha Spheres', 'Total SASA',
                     'Polar SASA', 'Apolar SASA', 'Volume', 'Mean local hydrophobic density',
                     'Apolar alpha sphere proportion', 'Hydrophobicity score', 'Charge score',
                     'Proportion of polar atoms', 'Alpha sphere density', 'Cent. of mass - Alpha Sphere max dist',
                     'biopython.csvmw', 'biopython.csvcharge_at_pH', 'QSOSW35', 'QSOSW37', 'QSOgrant37',
                     '_SecondaryStrC2', '_PolarityC1']
    adaboosting = ["Label",'Score', 'Number of Alpha Spheres', 'Mean alpha sphere radius',
                   'Druggability Score', 'biopython.csvaromaticty',
                   'biopython.csvisoelectric_point']
    gradient_boosting = ["Label",'Number of Alpha Spheres', 'biopython.csvmw', 'Druggability Score',
                         'QSOgrant34', 'Score', 'QSOgrant23', 'Mean local hydrophobic density',
                         'QSOgrant35', 'QSOSW17', 'QSOSW23', 'QSOgrant17', 'Volume', 'QSOSW35',
                         'GG', '_SecondaryStrD1075', 'QR', 'biopython.csvcharge_at_pH',
                         'Flexibility', 'biopython.csvaromaticty', 'Mean alpha sphere radius',
                         'biopython.csvinstability']

    average_f1_1 = []
    average_ave_pre = []
    average_roc_score = []


    for repeat in range(0,repeat_number):
        passer = model_fn(f"{model_path}/{repeat}_passer")
        gradient_boruta_model = model_fn(f"{model_path}/{repeat}_gradient_boruta")
        random_boruta_model = model_fn(f"{model_path}/{repeat}_random_boruta")
        adaboosting_model = model_fn(f"{model_path}/{repeat}_adaboosting")
        gradient_boosting_model = model_fn(f"{model_path}/{repeat}_gradient_boosting")

        average_pre_list = []
        roc_score_list = []
        f1_1_list = []

        for protein in test_list:
            filename = f"{current_dir}/data/test_one_by_one/test_{test_number}/{protein}.csv"
            data_test = pd.read_csv(filename, low_memory=False)
            passer_test = data_test.loc[:, fpocket_features]
            gradient_boruta_test = data_test.loc[:, gradient_boruta]
            random_boruta_test = data_test.loc[:, random_boruta]
            adaboosting_test = data_test.loc[:, adaboosting]
            gradient_boosting_test = data_test.loc[:, gradient_boosting]

            # find the location of 1 on y_true
            y_true = data_test["Label"]

            # make prediction for unseen data
            y_proba_0 = passer.predict_proba(passer_test)
            y_proba_1 = gradient_boruta_model.predict_proba(gradient_boruta_test)
            y_proba_2 = random_boruta_model.predict_proba(random_boruta_test)
            y_proba_3 = adaboosting_model.predict_proba(adaboosting_test)
            y_proba_4 = gradient_boosting_model.predict_proba(gradient_boosting_test)
            if model_selection == "MEF-AlloSite":
                mef_allo_site_proba = np.mean(
                    [y_proba_1.iloc[:, 1], y_proba_2.iloc[:, 1], y_proba_3.iloc[:, 1], y_proba_4.iloc[:, 1]], axis=0)
                selected_model_proba=mef_allo_site_proba
            elif model_selection == "PaSSer2.0":
                selected_model_proba=y_proba_0.iloc[:,1]
            highest_1_indices = np.argsort(selected_model_proba)[-1:]

            # Create an array of zeros
            y_pred_1 = np.zeros(len(selected_model_proba))

            # Set the highest 3 probabilities to 1
            y_pred_1[highest_1_indices] = 1

            # Here, only three metrics have been used to improve clearity
            aver_ = average_precision_score(y_true, selected_model_proba)
            roc_auc_score_ = roc_auc_score(y_true, selected_model_proba)
            f1_score_1 = f1_score(y_true=y_true, y_pred=y_pred_1)

            average_pre_list.append(aver_)
            roc_score_list.append(roc_auc_score_)
            f1_1_list.append(f1_score_1)

        average_ave_pre.append(statistics.mean(average_pre_list))
        average_roc_score.append(statistics.mean(roc_score_list))
        average_f1_1.append(statistics.mean(f1_1_list))

    return statistics.mean(average_ave_pre),statistics.mean(average_roc_score),statistics.mean(average_f1_1)

def validate_on_a_test(test_list,filename,test_number):
    #split_training_set(filename=filename,
    #                   test_number=test_number)

    #
    mef_av_pre,mef_av_roc,mef_av_f1=(validation_models(model_path=f"{current_dir}/AutogluonModels",
                            test_list=test_list,
                            test_number=test_number,
                            repeat_number=3,
                            model_selection="MEF-AlloSite")
          )
    par_av_pre,par_av_roc,par_av_f1=(validation_models(model_path=f"{current_dir}/AutogluonModels",
                            test_list=test_list,
                            test_number=test_number,
                            repeat_number=3,
                            model_selection="PaSSer2.0")
          )

    print("In our experiments, we repeated 51 times, then average scores")
    print(f"The average metric scores for 3 RUN (Test {test_number}):")
    results = [
        ("MEF-AlloSite", mef_av_pre, mef_av_roc, mef_av_f1),
        ("PaSSer2.0", par_av_pre, par_av_roc, par_av_f1)
    ]

    # Table headers
    headers = ["Model", "Average Pre", "ROC AUC", "F1 Score"]

    # Print table
    print(tabulate(results, headers=headers, tablefmt="pretty"))



validate_on_a_test(test_list=test_1,
                   filename="test_1_ready",
                   test_number=1)
validate_on_a_test(test_list=test_2,
                   filename="test_2_ready",
                   test_number=2)
validate_on_a_test(test_list=test_3,
                   filename="test_3_ready",
                   test_number=3)