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
# keep similar features
def keep_similar_features_p_value(data, similar_features,threshold=0.95):
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

    # Filter correlated features based on similarity to any element in the given list
    filtered_corr_features = set()
    for feature in corr_features:
        for similar_feature in similar_features:
            if similar_feature in feature:
                filtered_corr_features.add(feature)

    # Keep similar features and remove the rest
    selected_features = [feature for feature in data.columns if feature in filtered_corr_features]

    return selected_features


