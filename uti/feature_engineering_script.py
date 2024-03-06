import os
import sys
import glob
import os
from script import docking_based_cavity_detection
from script import psvina
from script import protein_uti
from script import vina_uti
from script import smina_feature
import pandas as pd
import numpy as np
import csv
import subprocess
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from script import protein_similarity_search
from script import protein_uti
from sklearn import model_selection
from sklearn.ensemble import HistGradientBoostingClassifier

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

# two-method: 1) remove similar features and then select features
# 2) select features then remove similar features


# Generate example data
"""filename="/home/yavuz/yavuz_proje/allosteric_binding_site_3_may/data/training_data/train_7.csv"
#filename="/home/yavuz/yavuz_proje/allosteric_binding_site_3_may/data/training_data/1XJE.csv"
data = pd.read_csv("%s"%filename, index_col=0)
column_means = data.mean()
data.fillna(column_means, inplace=True)

selection_list = ['sequence_len', 'CYS', 'ASP', 'SER', 'GLN', 'LYS', 'ILE', 'PRO', 'THR', 'PHE', 'ASN', 'GLY', 'HIS', 'LEU', 'ARG', 'TRP', 'ALA', 'VAL', 'GLU', 'TYR', 'MET', 'mw', 'aromaticty', 'instability', 'isoelectric_point', 'helix', 'turn', 'sheet', 'gravy', 'charge_at_pH', 'molar_extinction_coefficient_0', 'molar_extinction_coefficient_1', 'protein_volume', 'total_atom', 'atom_density', 'total_area',
                  'pocket_number', 'Score', 'Druggability Score', 'Number of Alpha Spheres', 'Total SASA', 'Polar SASA', 'Apolar SASA', 'Volume', 'Mean local hydrophobic density', 'Mean alpha sphere radius', 'Mean alp. sph. solvent access', 'Apolar alpha sphere proportion', 'Hydrophobicity score', 'Volume score', 'Polarity score', 'Charge score', 'Proportion of polar atoms', 'Alpha sphere density', 'Cent. of mass - Alpha Sphere max dist', 'Flexibility',
                  '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
                  'Total_area', 'Polar', 'Apolar', 'Nitrogen', 'Carbon', 'Oxygen', 'Sulfur', '_PolarizabilityC1', '_PolarizabilityC2', '_PolarizabilityC3', '_SolventAccessibilityC1', '_SolventAccessibilityC2', '_SolventAccessibilityC3', '_SecondaryStrC1', '_SecondaryStrC2', '_SecondaryStrC3', '_ChargeC1', '_ChargeC2', '_ChargeC3', '_PolarityC1', '_PolarityC2', '_PolarityC3', '_NormalizedVDWVC1', '_NormalizedVDWVC2', '_NormalizedVDWVC3', '_HydrophobicityC1', '_HydrophobicityC2', '_HydrophobicityC3', '_PolarizabilityT12', '_PolarizabilityT13', '_PolarizabilityT23', '_SolventAccessibilityT12', '_SolventAccessibilityT13', '_SolventAccessibilityT23', '_SecondaryStrT12', '_SecondaryStrT13', '_SecondaryStrT23', '_ChargeT12', '_ChargeT13', '_ChargeT23', '_PolarityT12', '_PolarityT13', '_PolarityT23', '_NormalizedVDWVT12', '_NormalizedVDWVT13', '_NormalizedVDWVT23', '_HydrophobicityT12', '_HydrophobicityT13', '_HydrophobicityT23', '_PolarizabilityD1001', '_PolarizabilityD1025', '_PolarizabilityD1050', '_PolarizabilityD1075', '_PolarizabilityD1100', '_PolarizabilityD2001', '_PolarizabilityD2025', '_PolarizabilityD2050', '_PolarizabilityD2075', '_PolarizabilityD2100', '_PolarizabilityD3001', '_PolarizabilityD3025', '_PolarizabilityD3050', '_PolarizabilityD3075', '_PolarizabilityD3100', '_SolventAccessibilityD1001', '_SolventAccessibilityD1025', '_SolventAccessibilityD1050', '_SolventAccessibilityD1075', '_SolventAccessibilityD1100', '_SolventAccessibilityD2001', '_SolventAccessibilityD2025', '_SolventAccessibilityD2050', '_SolventAccessibilityD2075', '_SolventAccessibilityD2100', '_SolventAccessibilityD3001', '_SolventAccessibilityD3025', '_SolventAccessibilityD3050', '_SolventAccessibilityD3075', '_SolventAccessibilityD3100', '_SecondaryStrD1001', '_SecondaryStrD1025', '_SecondaryStrD1050', '_SecondaryStrD1075', '_SecondaryStrD1100', '_SecondaryStrD2001', '_SecondaryStrD2025', '_SecondaryStrD2050', '_SecondaryStrD2075', '_SecondaryStrD2100', '_SecondaryStrD3001', '_SecondaryStrD3025', '_SecondaryStrD3050', '_SecondaryStrD3075', '_SecondaryStrD3100', '_ChargeD1001', '_ChargeD1025', '_ChargeD1050', '_ChargeD1075', '_ChargeD1100', '_ChargeD2001', '_ChargeD2025', '_ChargeD2050', '_ChargeD2075', '_ChargeD2100', '_ChargeD3001', '_ChargeD3025', '_ChargeD3050', '_ChargeD3075', '_ChargeD3100', '_PolarityD1001', '_PolarityD1025', '_PolarityD1050', '_PolarityD1075', '_PolarityD1100', '_PolarityD2001', '_PolarityD2025', '_PolarityD2050', '_PolarityD2075', '_PolarityD2100', '_PolarityD3001', '_PolarityD3025', '_PolarityD3050', '_PolarityD3075', '_PolarityD3100', '_NormalizedVDWVD1001', '_NormalizedVDWVD1025', '_NormalizedVDWVD1050', '_NormalizedVDWVD1075', '_NormalizedVDWVD1100', '_NormalizedVDWVD2001', '_NormalizedVDWVD2025', '_NormalizedVDWVD2050', '_NormalizedVDWVD2075', '_NormalizedVDWVD2100', '_NormalizedVDWVD3001', '_NormalizedVDWVD3025', '_NormalizedVDWVD3050', '_NormalizedVDWVD3075', '_NormalizedVDWVD3100', '_HydrophobicityD1001', '_HydrophobicityD1025', '_HydrophobicityD1050', '_HydrophobicityD1075', '_HydrophobicityD1100', '_HydrophobicityD2001', '_HydrophobicityD2025', '_HydrophobicityD2050', '_HydrophobicityD2075', '_HydrophobicityD2100', '_HydrophobicityD3001', '_HydrophobicityD3025', '_HydrophobicityD3050', '_HydrophobicityD3075', '_HydrophobicityD3100', 'GAFF TOTAL BOND STRETCHING ENERGY', 'GAFF TOTAL ANGLE BENDING ENERGY', 'GAFF TOTAL TORSIONAL ENERGY', 'GAFF TOTAL IMPROPER-TORSIONAL ENERGY', 'GAFF TOTAL VAN DER WAALS ENERGY', 'GAFF TOTAL ELECTROSTATIC ENERGY', 'GAFF TOTAL ENERGY', 'UFF TOTAL BOND STRETCHING ENERGY', 'UFF TOTAL ANGLE BENDING ENERGY', 'UFF TOTAL TORSIONAL ENERGY', 'UFF TOTAL VAN DER WAALS ENERGY', 'UFF TOTAL ENERGY', 'Ghemical TOTAL BOND STRETCHING ENERGY', 'Ghemical TOTAL ANGLE BENDING ENERGY', 'Ghemical TOTAL TORSIONAL ENERGY', 'Ghemical TOTAL VAN DER WAALS ENERGY', 'Ghemical TOTAL ELECTROSTATIC ENERGY', 'Ghemical TOTAL ENERGY', 'MMFF94 TOTAL BOND STRETCHING ENERGY', 'MMFF94 TOTAL ANGLE BENDING ENERGY', 'MMFF94 TOTAL STRETCH BENDING ENERGY', 'MMFF94 TOTAL TORSIONAL ENERGY', 'MMFF94 TOTAL OUT-OF-PLANE BENDING ENERGY', 'MMFF94 TOTAL VAN DER WAALS ENERGY', 'MMFF94 TOTAL ELECTROSTATIC ENERGY', 'MMFF94 TOTAL ENERGY']
selection_list = ['pocket_number', 'Score', 'Druggability Score', 'Number of Alpha Spheres', 'Total SASA', 'Polar SASA',
                  'Apolar SASA', 'Volume', 'Mean local hydrophobic density', 'Mean alpha sphere radius',
                  'Mean alp. sph. solvent access', 'Apolar alpha sphere proportion', 'Hydrophobicity score',
                  'Volume score', 'Polarity score', 'Charge score', 'Proportion of polar atoms', 'Alpha sphere density',
                  'Cent. of mass - Alpha Sphere max dist', 'Flexibility',

'sequence_len', 'CYS', 'ASP', 'SER', 'GLN', 'LYS', 'ILE', 'PRO', 'THR', 'PHE', 'ASN', 'GLY', 'HIS', 'LEU', 'ARG', 'TRP', 'ALA', 'VAL', 'GLU', 'TYR', 'MET', 'mw', 'aromaticty', 'instability', 'isoelectric_point', 'helix', 'turn', 'sheet', 'gravy', 'charge_at_pH', 'molar_extinction_coefficient_0', 'molar_extinction_coefficient_1', 'protein_volume', 'total_atom', 'atom_density', 'total_area',
                  '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
                  'Total_area', 'Polar', 'Apolar', 'Nitrogen', 'Carbon', 'Oxygen', 'Sulfur', '_PolarizabilityC1', '_PolarizabilityC2', '_PolarizabilityC3', '_SolventAccessibilityC1', '_SolventAccessibilityC2', '_SolventAccessibilityC3', '_SecondaryStrC1', '_SecondaryStrC2', '_SecondaryStrC3', '_ChargeC1', '_ChargeC2', '_ChargeC3', '_PolarityC1', '_PolarityC2', '_PolarityC3', '_NormalizedVDWVC1', '_NormalizedVDWVC2', '_NormalizedVDWVC3', '_HydrophobicityC1', '_HydrophobicityC2', '_HydrophobicityC3', '_PolarizabilityT12', '_PolarizabilityT13', '_PolarizabilityT23', '_SolventAccessibilityT12', '_SolventAccessibilityT13', '_SolventAccessibilityT23', '_SecondaryStrT12', '_SecondaryStrT13', '_SecondaryStrT23', '_ChargeT12', '_ChargeT13', '_ChargeT23', '_PolarityT12', '_PolarityT13', '_PolarityT23', '_NormalizedVDWVT12', '_NormalizedVDWVT13', '_NormalizedVDWVT23', '_HydrophobicityT12', '_HydrophobicityT13', '_HydrophobicityT23', '_PolarizabilityD1001', '_PolarizabilityD1025', '_PolarizabilityD1050', '_PolarizabilityD1075', '_PolarizabilityD1100', '_PolarizabilityD2001', '_PolarizabilityD2025', '_PolarizabilityD2050', '_PolarizabilityD2075', '_PolarizabilityD2100', '_PolarizabilityD3001', '_PolarizabilityD3025', '_PolarizabilityD3050', '_PolarizabilityD3075', '_PolarizabilityD3100', '_SolventAccessibilityD1001', '_SolventAccessibilityD1025', '_SolventAccessibilityD1050', '_SolventAccessibilityD1075', '_SolventAccessibilityD1100', '_SolventAccessibilityD2001', '_SolventAccessibilityD2025', '_SolventAccessibilityD2050', '_SolventAccessibilityD2075', '_SolventAccessibilityD2100', '_SolventAccessibilityD3001', '_SolventAccessibilityD3025', '_SolventAccessibilityD3050', '_SolventAccessibilityD3075', '_SolventAccessibilityD3100', '_SecondaryStrD1001', '_SecondaryStrD1025', '_SecondaryStrD1050', '_SecondaryStrD1075', '_SecondaryStrD1100', '_SecondaryStrD2001', '_SecondaryStrD2025', '_SecondaryStrD2050', '_SecondaryStrD2075', '_SecondaryStrD2100', '_SecondaryStrD3001', '_SecondaryStrD3025', '_SecondaryStrD3050', '_SecondaryStrD3075', '_SecondaryStrD3100', '_ChargeD1001', '_ChargeD1025', '_ChargeD1050', '_ChargeD1075', '_ChargeD1100', '_ChargeD2001', '_ChargeD2025', '_ChargeD2050', '_ChargeD2075', '_ChargeD2100', '_ChargeD3001', '_ChargeD3025', '_ChargeD3050', '_ChargeD3075', '_ChargeD3100', '_PolarityD1001', '_PolarityD1025', '_PolarityD1050', '_PolarityD1075', '_PolarityD1100', '_PolarityD2001', '_PolarityD2025', '_PolarityD2050', '_PolarityD2075', '_PolarityD2100', '_PolarityD3001', '_PolarityD3025', '_PolarityD3050', '_PolarityD3075', '_PolarityD3100', '_NormalizedVDWVD1001', '_NormalizedVDWVD1025', '_NormalizedVDWVD1050', '_NormalizedVDWVD1075', '_NormalizedVDWVD1100', '_NormalizedVDWVD2001', '_NormalizedVDWVD2025', '_NormalizedVDWVD2050', '_NormalizedVDWVD2075', '_NormalizedVDWVD2100', '_NormalizedVDWVD3001', '_NormalizedVDWVD3025', '_NormalizedVDWVD3050', '_NormalizedVDWVD3075', '_NormalizedVDWVD3100', '_HydrophobicityD1001', '_HydrophobicityD1025', '_HydrophobicityD1050', '_HydrophobicityD1075', '_HydrophobicityD1100', '_HydrophobicityD2001', '_HydrophobicityD2025', '_HydrophobicityD2050', '_HydrophobicityD2075', '_HydrophobicityD2100', '_HydrophobicityD3001', '_HydrophobicityD3025', '_HydrophobicityD3050', '_HydrophobicityD3075', '_HydrophobicityD3100', 'GAFF TOTAL BOND STRETCHING ENERGY', 'GAFF TOTAL ANGLE BENDING ENERGY', 'GAFF TOTAL TORSIONAL ENERGY', 'GAFF TOTAL IMPROPER-TORSIONAL ENERGY', 'GAFF TOTAL VAN DER WAALS ENERGY', 'GAFF TOTAL ELECTROSTATIC ENERGY', 'GAFF TOTAL ENERGY', 'UFF TOTAL BOND STRETCHING ENERGY', 'UFF TOTAL ANGLE BENDING ENERGY', 'UFF TOTAL TORSIONAL ENERGY', 'UFF TOTAL VAN DER WAALS ENERGY', 'UFF TOTAL ENERGY', 'Ghemical TOTAL BOND STRETCHING ENERGY', 'Ghemical TOTAL ANGLE BENDING ENERGY', 'Ghemical TOTAL TORSIONAL ENERGY', 'Ghemical TOTAL VAN DER WAALS ENERGY', 'Ghemical TOTAL ELECTROSTATIC ENERGY', 'Ghemical TOTAL ENERGY', 'MMFF94 TOTAL BOND STRETCHING ENERGY', 'MMFF94 TOTAL ANGLE BENDING ENERGY', 'MMFF94 TOTAL STRETCH BENDING ENERGY', 'MMFF94 TOTAL TORSIONAL ENERGY', 'MMFF94 TOTAL OUT-OF-PLANE BENDING ENERGY', 'MMFF94 TOTAL VAN DER WAALS ENERGY', 'MMFF94 TOTAL ELECTROSTATIC ENERGY', 'MMFF94 TOTAL ENERGY'
                  ]
print(len(selection_list))
X = data.loc[:, selection_list]
y = data["label"]
data=X"""
#######################################################################################################################################
#######################################################################################################################################
# ensemble remove similar features
def remove_similar_features(data,
                            threshold = 0.95):
    import pandas as pd
    import numpy as np
    from sklearn.feature_selection import SelectKBest, f_classif


    # Step 2: Compute pairwise correlations
    correlation_matrix = data.corr().abs()

    # Step 3: Remove highly correlated features
     # Set the correlation threshold
    corr_features = set()  # Set to store correlated feature pairs

    # Find correlated features and add them to the set
    for i in range(len(correlation_matrix.columns)):
        for j in range(i):
            if correlation_matrix.iloc[i, j] >= threshold:
                corr_features.add(correlation_matrix.columns[i])
                corr_features.add(correlation_matrix.columns[j])

    # Remove the correlated features from the DataFrame
    data_filtered = data.drop(columns=corr_features)
    selected_features=(data_filtered.columns.tolist())
    print(len(selected_features))
    print(selected_features)
    return selected_features

# p-value
def remove_similar_features_p_value(data, threshold=0.95):
    import pandas as pd
    from sklearn.feature_selection import VarianceThreshold

    # Compute pairwise correlations
    correlation_matrix = data.corr().abs()

    # Identify highly correlated features
    corr_features = set()
    for i in range(len(correlation_matrix.columns)):
        for j in range(i):
            if correlation_matrix.iloc[i, j] >= threshold:
                corr_features.add(correlation_matrix.columns[i])
                corr_features.add(correlation_matrix.columns[j])

    # Drop highly correlated features
    data_filtered = data.drop(columns=corr_features)

    # Apply VarianceThreshold to remove low-variance features
    selector = VarianceThreshold()
    data_filtered = selector.fit_transform(data_filtered)

    # Get the selected feature indices
    selected_indices = selector.get_support(indices=True)

    # Get the selected feature names
    selected_features = data.columns[selected_indices].tolist()

    return selected_features

# Recursive Feature Elimination (RFE)
def remove_similar_features_rfe(data, num_features=100):
    import pandas as pd
    from sklearn.feature_selection import RFE
    from sklearn.ensemble import RandomForestRegressor

    # Step 1: Separate features and target variable
    y = data['label']
    X = data.drop(columns=['label'])

    # Step 2: Perform Recursive Feature Elimination (RFE)
    estimator = RandomForestRegressor()  # You can change the estimator to a different model if desired
    rfe = RFE(estimator, n_features_to_select=num_features)
    X_selected = rfe.fit_transform(X, y)

    # Step 3: Get the selected feature indices
    selected_indices = rfe.get_support(indices=True)

    # Step 4: Get the selected feature names
    selected_features = X.columns[selected_indices].tolist()

    print(len(selected_features))
    print(selected_features)

    return selected_features

# select k
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_regression

import pandas as pd
from sklearn.feature_selection import SelectKBest, f_regression

def remove_similar_features_selectkbest(data, num_features=100):
    # Step 1: Separate features and target variable
    y = data['label']
    X = data.drop(columns=['label'])

    # Step 2: Check if num_features is within a valid range
    if num_features == 'all':
        k_value = X.shape[1]  # Select all features
    else:
        k_value = min(num_features, X.shape[1])

    # Step 3: Perform feature selection using SelectKBest
    selector = SelectKBest(score_func=f_regression, k=k_value)
    X_selected = selector.fit_transform(X, y)

    # Step 4: Get the selected feature indices
    selected_indices = selector.get_support(indices=True)

    # Step 5: Get the selected feature names
    selected_features = X.columns[selected_indices].tolist()

    print(len(selected_features))
    print(selected_features)

    # Step 6: Transform the dataset using the selected features
    data_filtered = data.iloc[:, selected_indices]

    return selected_features

# Example usage:
#selected_features = remove_similar_features(data, threshold=0.95, num_features=10)


# show correlation
def heat_map(data):
    # Pairwise correlation
    correlation_matrix = data.corr()
    print("\nCorrelation Matrix:")
    print(correlation_matrix)

    import seaborn as sns
    import matplotlib.pyplot as plt

    # Assuming you have a correlation matrix called 'correlation_matrix'
    # You can calculate it using pandas DataFrame's 'corr' method

    # Create a custom colormap
    cmap = sns.diverging_palette(250, 10, as_cmap=True)

    # Create a correlation matrix plot
    plt.figure(figsize=(100, 80))

    # Increase figure size and rotate x-axis labels
    plt.figure(figsize=(12, 10))
    plt.xticks(rotation=90)

    # Decrease font size
    sns.set(font_scale=0.8)
    sns.heatmap(correlation_matrix, cmap=cmap, annot=True, fmt=".2f", vmin=-1, vmax=1, center=0)
    plt.title('Correlation Matrix')
    plt.show()
    plt.savefig('correlation_matrix.png')

def feature_set_similarty_detection(feature_set_selections):
    import numpy as np


    # Create an empty dictionary to store the Q statistics
    q_statistics = {}

    # Iterate through each pair of features
    for i in range(len(feature_set_selections)):
        for j in range(i + 1, len(feature_set_selections)):
            selection1 = set(feature_set_selections[i])
            selection2 = set(feature_set_selections[j])

            # Calculate the intersection and union of the two selections
            intersection = len(selection1.intersection(selection2))
            union = len(selection1.union(selection2))

            # Calculate the Q statistic for the pair
            q_stat = intersection / union

            # Store the Q statistic in the dictionary
            q_statistics[(i, j)] = q_stat

    # Print the pair-wise Q statistics
    for pair, q_stat in q_statistics.items():
        print(f"Pair {pair}: Q-statistic = {q_stat}")
    return q_statistics

#######################################################################################################################################
#######################################################################################################################################
# ensemble feature selection
# 1 Ensemble Method: Boruta
def boruta_selection_orginal(X,y):
    import random
    random.seed(42)  # Set a random seed for reproducibility
    from boruta import BorutaPy
    from sklearn.ensemble import RandomForestClassifier

    # Assuming you have features 'X' and target variable 'y'

    # Create a RandomForestClassifier object
    clf = RandomForestClassifier()

    # Create a Boruta object
    boruta_selector = BorutaPy(estimator=clf, n_estimators='auto', max_iter=100)

    # Fit the selector to your data
    boruta_selector.fit(X.values, y.values)

    # Get the indices of the selected features
    selected_features_indices = boruta_selector.support_
    selected_features = X.columns[selected_features_indices].tolist()
    print(selected_features)
    print(len(selected_features))
    return selected_features
def boruta_selection(X,y):
    import random
    random.seed(42)  # Set a random seed for reproducibility
    from boruta import BorutaPy
    from sklearn.ensemble import RandomForestClassifier

    # Assuming you have features 'X' and target variable 'y'

    # Create a RandomForestClassifier object
    clf = RandomForestClassifier()

    # Create a Boruta object
    boruta_selector = BorutaPy(estimator=clf, n_estimators="auto")

    # Fit the selector to your data
    boruta_selector.fit(X.values, y.values)

    # Get the indices of the selected features
    selected_features_indices = boruta_selector.support_
    selected_features = X.columns[selected_features_indices].tolist()
    print(selected_features)
    print(len(selected_features))
    return selected_features
def boruta_selection_default(X,y):
    import random
    random.seed(42)  # Set a random seed for reproducibility
    from boruta import BorutaPy
    from sklearn.ensemble import RandomForestClassifier

    # Assuming you have features 'X' and target variable 'y'

    # Create a RandomForestClassifier object
    clf = RandomForestClassifier()

    # Create a Boruta object
    boruta_selector = BorutaPy(estimator=clf,two_step=True,random_state=42)

    # Fit the selector to your data
    boruta_selector.fit(X.values, y.values)

    # Get the indices of the selected features
    selected_features_indices = boruta_selector.support_
    selected_features = X.columns[selected_features_indices].tolist()
    print(selected_features)
    print(len(selected_features))
    return selected_features
def boruta_selection_more_feature(X,y):
    import random
    random.seed(42)  # Set a random seed for reproducibility
    from boruta import BorutaPy
    from sklearn.ensemble import RandomForestClassifier

    # Assuming you have features 'X' and target variable 'y'

    # Create a RandomForestClassifier object
    clf = RandomForestClassifier(n_estimators=5)

    # Create a Boruta object
    boruta_selector = BorutaPy(estimator=clf, random_state=42, alpha=0.001,n_estimators=100)

    # Fit the selector to your data
    boruta_selector.fit(X.values, y.values)

    # Get the indices of the selected features
    selected_features_indices = boruta_selector.support_

    # check ranking of features
    boruta_ranking=boruta_selector.ranking_
    selected_features = X.columns[selected_features_indices].tolist()
    boruta_ranking = X.columns[boruta_ranking].tolist()
    print(selected_features)
    print(len(selected_features))

    print("boruta_ranking")
    print(boruta_ranking)
    print(len(boruta_ranking))
    return selected_features

# 2 Embedded Method: Lasso (L1 Regularization)
def lasso_feature_selection(X,y):
    from sklearn.linear_model import Lasso

    lasso = Lasso(alpha=0.05)
    lasso.fit(X, y)

    #coef_dictionary=
    # Select features with non-zero coefficients
    selected_features = X.columns[lasso.coef_ != 0].tolist()
    print("selected_features")
    print(selected_features)
    return selected_features

def lasso_feature_selection_with_number(X, y, num_features=100):
    from sklearn.linear_model import Lasso

    lasso = Lasso(alpha=0.05)
    lasso.fit(X, y)

    # Select features with non-zero coefficients
    selected_features_indices = (lasso.coef_ != 0)
    selected_features_coefficients = lasso.coef_[selected_features_indices]

    # Sort coefficients in descending order and select the top 100 features
    sorted_indices = selected_features_coefficients.argsort()[::-1]
    top_100_indices = sorted_indices[:num_features]

    selected_features = X.columns[selected_features_indices][top_100_indices].tolist()
    print("selected_features")
    print(selected_features)
    return selected_features


# 3 Embeded all
def embeded_all(X,y):
    from sklearn.datasets import load_boston
    from sklearn.linear_model import Lasso, Ridge, ElasticNet
    from sklearn.ensemble import RandomForestRegressor


    # Lasso feature selection
    lasso = Lasso(alpha=0.1)
    lasso.fit(X, y)
    lasso_selected_features = X.columns[lasso.coef_ != 0].tolist()

    # Ridge Regression feature selection
    ridge = Ridge(alpha=0.1)
    ridge.fit(X, y)
    ridge_selected_features = X.columns[ridge.coef_ != 0].tolist()

    # Elastic Net feature selection
    elastic_net = ElasticNet(alpha=0.1, l1_ratio=0.5)
    elastic_net.fit(X, y)
    elastic_net_selected_features = X.columns[elastic_net.coef_ != 0].tolist()

    # Random Forest feature importance
    rf = RandomForestRegressor()
    rf.fit(X, y)
    rf_feature_importances = rf.feature_importances_
    mean_importance = np.mean(rf_feature_importances)
    rf_selected_features = X.columns[rf_feature_importances > mean_importance].tolist()

    # Print selected features
    print("Lasso Selected Features:", lasso_selected_features)
    print("Ridge Regression Selected Features:", ridge_selected_features)
    print("Elastic Net Selected Features:", elastic_net_selected_features)
    print("Random Forest Selected Features:", rf_selected_features)
    return lasso_selected_features,ridge_selected_features,elastic_net_selected_features,rf_selected_features
# 4 p-value
def p_value_selection(X,y):
    from sklearn.feature_selection import SelectFpr, chi2,SelectFwe,SelectFpr
    from sklearn.feature_selection import GenericUnivariateSelect, chi2
    # Assuming you have features 'X' and target variable 'y'

    # Set the desired alpha level for false positive rate
    alpha = 0.00005  # Significance level of 0.05

    # Create a SelectFwe object with the desired score function
    selector = SelectFwe(alpha=alpha)

    # Fit the selector to your data
    selector.fit(X, y)

    # Get the indices of the selected features
    selected_features_indices = selector.get_support(indices=True)

    # Select the corresponding feature names
    selected_features = X.columns[selected_features_indices].tolist()
    print(len(selected_features))
    print(selected_features)
    return  selected_features
from sklearn.feature_selection import SelectFwe

def p_value_selection_with_number(X, y, n=100):
    # Assuming you have features 'X' and target variable 'y'

    # Set the desired alpha level for false positive rate
    alpha = 0.05  # Significance level of 0.05

    # Create a SelectFwe object with the desired score function and alpha
    selector = SelectFwe(alpha=alpha)

    # Fit the selector to your data
    selector.fit(X, y)

    # Get the p-values of the features
    p_values = selector.pvalues_

    # Sort p-values in ascending order and get the top 'n' feature indices
    sorted_indices = p_values.argsort()
    top_n_indices = sorted_indices[:n]

    # Select the corresponding feature names
    selected_features = X.columns[top_n_indices].tolist()
    print("Selected Features based on p-values:")
    print(selected_features)
    return selected_features

from sklearn.datasets import load_boston
from sklearn.linear_model import Lasso, Ridge, ElasticNet
from sklearn.ensemble import RandomForestRegressor
# 5 ridge_feature_selection
def Ridge_feature_selection(X, y):
    # Ridge Regression feature selection
    ridge = Ridge()
    ridge.fit(X, y)

    # Get the absolute magnitudes of the coefficients
    coef_magnitudes = abs(ridge.coef_)

    # Sort the feature magnitudes in descending order
    sorted_indices = np.argsort(coef_magnitudes)[::-1]

    # Select the top 100 features
    top_indices = sorted_indices[:100]

    # Get the selected feature names
    ridge_selected_features = X.columns[top_indices].tolist()

    print("Ridge Regression Selected Features:", ridge_selected_features)
    return ridge_selected_features
# 6 elastic_net
def Elastic_feature_selection(X,y):
    from sklearn.datasets import load_boston
    from sklearn.linear_model import Lasso, Ridge, ElasticNet
    from sklearn.ensemble import RandomForestRegressor

    # Elastic Net feature selection
    elastic_net = ElasticNet(alpha=0.1, l1_ratio=0.5)
    elastic_net.fit(X, y)
    elastic_net_selected_features = X.columns[elastic_net.coef_ != 0].tolist()

    print("Elastic Net Selected Features:", elastic_net_selected_features)
    return elastic_net_selected_features

def elastic_feature_selection_with_number(X, y, num_features=100):
    from sklearn.linear_model import ElasticNet

    # Elastic Net feature selection
    elastic_net = ElasticNet(alpha=0.1, l1_ratio=0.5)
    elastic_net.fit(X, y)

    # Select features with non-zero coefficients
    selected_features_indices = (elastic_net.coef_ != 0)
    selected_features_coefficients = elastic_net.coef_[selected_features_indices]

    # Sort coefficients in descending order and select the top 100 features
    sorted_indices = selected_features_coefficients.argsort()[::-1]
    top_100_indices = sorted_indices[:num_features]

    selected_features = X.columns[selected_features_indices][top_100_indices].tolist()
    print("Elastic Net Selected Features:", selected_features)
    return selected_features


from sklearn.ensemble import RandomForestRegressor
import numpy as np
# RandomForestRegressor
def RandomForestRegressor_feature_selection(X, y,number=100):
    # Random Forest feature importance
    rf = RandomForestRegressor()
    rf.fit(X, y)
    rf_feature_importances = rf.feature_importances_

    # Sort feature importances in descending order
    sorted_indices = np.argsort(rf_feature_importances)[::-1]

    # Select the top 100 features
    top_indices = sorted_indices[:number]

    # Get the selected feature names
    rf_selected_features = X.columns[top_indices].tolist()

    print("Random Forest Selected Features:", rf_selected_features)
    return rf_selected_features
#######################################################################################################################################
#######################################################################################################################################
# aggregation
"""feature_set_selections = [
        ['feature1', 'feature2'],
        ['feature2', 'feature4'],
        ['feature1', 'feature3', 'feature5'],
        # ... more subsets
    ]"""
# 1) select_half_of_them
def select_half_of_them(feature_set_selections):
    from collections import defaultdict

    # Count the number of times each feature is selected
    feature_counts = defaultdict(int)
    for selection in feature_set_selections:
        for feature in selection:
            feature_counts[feature] += 1

    # Set a threshold or voting rule to select features
    threshold = len(feature_set_selections) * 0.5  # Select features chosen by at least 50% of the methods
    selected_features = [feature for feature, count in feature_counts.items() if count >= threshold]

    # Print the final selected features
    print(selected_features)

# 2) p-value
def p_value(feature_set_selections):
    from collections import Counter
    from scipy.stats import chi2_contingency

    def aggregate_feature_sets(feature_sets, p_value_threshold):
        # Flatten all feature sets into a single list
        all_features = [feature for feature_set in feature_sets for feature in feature_set]

        # Count the occurrences of each feature
        feature_counts = Counter(all_features)

        # Calculate p-value for each feature
        p_values = {}
        for feature, count in feature_counts.items():
            observed = [[count, len(feature_sets) - count], [len(feature_sets) - count, count]]
            chi2, p, _, _ = chi2_contingency(observed)
            p_values[feature] = p

        # Select features with p-values below the threshold
        aggregated_features = [feature for feature, p_value in p_values.items() if p_value < p_value_threshold]

        return aggregated_features

    # Set the p-value threshold
    p_value_threshold = 0.05

    # Aggregate the feature sets using the p-value threshold
    aggregated_features = aggregate_feature_sets(feature_set_selections, p_value_threshold)

    print(aggregated_features)
    return  aggregated_features

# 3) mean and median
def borda_cound_mean_median(feature_set_selections):
    from collections import defaultdict

    def borda_count(feature_set_selections):
        # Count the number of times each feature is selected
        feature_counts = defaultdict(int)
        for selection in feature_set_selections:
            for i, feature in enumerate(selection):
                feature_counts[feature] += len(selection) - i

        # Calculate the Borda scores
        borda_scores = {feature: score for feature, score in feature_counts.items()}

        return borda_scores

    def remove_features_by_threshold(borda_scores, threshold):
        # Calculate the mean and median of the Borda scores
        borda_values = list(borda_scores.values())
        mean = sum(borda_values) / len(borda_values)
        median = sorted(borda_values)[len(borda_values) // 2]

        # Remove features based on the specified thresholds
        selected_features = [feature for feature, score in borda_scores.items()
                             if score >= threshold ]

        return selected_features

    # Calculate Borda scores
    borda_scores = borda_count(feature_set_selections)

    # Set the thresholds
    mean_threshold = sum(borda_scores.values()) / len(borda_scores)
    median_threshold = sorted(borda_scores.values())[len(borda_scores) // 2]
    #print("mean_threshold",mean_threshold)
    #print("median_threshold",median_threshold)
    # Remove features based on the thresholds
    selected_features_mean = remove_features_by_threshold(borda_scores, mean_threshold)
    print("Selected features based on mean threshold:", selected_features_mean)

    selected_features_median = remove_features_by_threshold(borda_scores, median_threshold)
    print("Selected features based on median threshold:", selected_features_median)

    print("borda_scores",borda_scores)
    return  selected_features_mean, selected_features_median,borda_scores

