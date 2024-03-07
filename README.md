<h1> MEF-AlloSite: An accurate and robust Multimodel Ensemble Feature selection for the Allosteric Site identification model </h1>

This readme file documents all of the required steps to run MEF-AlloSite.

Note that the code was implemented and tested on a Linux operating system only.

<h2>How to set up the environment</h2>

**Before installation**

(i)-MathFeature

-Please put MathFeature files under the bin directory.

MathFeature:

https://github.com/Bonidia/MathFeature

**1- Conda Installation**:
-We have provided an Anaconda environment file for easy set-up. If you do not have Anaconda installed, you can get Miniconda from [HERE](https://docs.anaconda.com/free/miniconda/).

Then, install ```bash mamba```:
```bash
conda install -c conda-forge mamba
```

If an individual encounters difficulties in implementing mamba in their primary setting, proceed to establish a fresh environment for mamba by following these steps:
```bash
conda create -n mamba-env -c conda-forge mamba
conda activate mamba-env
```

**2- Conda Env Formation**:
- To establish a conda environment, please execute the following code:
```bash
mamba env create -n mef_allosite -f environment.yml
```
The conda environment already encompasses Fpocket [1], Biopython [2], PyBioMed [3], and MathFeature [4].

**3- Activate Conda Env**:
- Before running MEF-AlloSite, enable the conda environment.
```bash
mamba activate mef_allosite
```
It should be noted that to activate the mef_allosite environment, it may be necessary to close and re-open the terminal.

Next, proceed to install supplementary packages using the pip command as follows:
```bash
pip install -r requirements.txt
```


**4 - Make predictions**
- CSV files for Test 1, 2, and 3 exist in the repo as well as trained models. Therefore, you can make predictions on our test cases by executing:
```bash
python3 7-Make_predictions.py
```

**To train new models/or make predictions on novel proteins**
**Run Python Scripts one by one**:
- The script contains a list of proteins for training, testing, and testing. Consequently, the execution of the script results in the retrieval of all proteins utilised in the study.
```bash
python3 1-Download_PDBs.py
```

- The objective of the cavity detection-based technique is to identify and locate cavities on target proteins. Fpocket [1], similar to PaSSer2.0 [5], is utilised by MEF-AlloSite. The script does protein cleaning to subsequently find cavities.
```bash
python3 2-Find_pockets.py
```

- Please execute the script below to gather the Fpocket [1] feature. The script retrieves data from .info output files. (PaSSer2.0 [5] labelling strategy has been used).
```bash
python3 3-Prepare_fpocket_features.py
```

- The script checks whether an atom is missing structure on pockets.
- To possibly increase the performance of amino acid-based features, the number of residues side by pocket can be increased using the script. 

NOTE: FPOCKET results have been used in the research to prevent any bias in comparison analysis.
```bash
python3 4-Fix_pockets.py
```

- The script takes amino acid from a FIXED pocket and prepares Biopython [2], PyBioMed [3] and MathFeature [4] features. Then, merge with Fpocket [1] features.
```bash
python3 5-Prepare_Other_Features.py
```

- The prepared CSV file has been used to train models using AutoGluon three times using different splits of the training set. In the research, 51 different splits have been used to train models.
```bash
python3 6-train_a_model.py
```

- The script uses models to generate forecasts on three separate occasions and subsequently calculates the average of these predictions to provide the final performance metrics, such as F1 score, Average Precision and ROC AUC score.
```bash
python3 7-Make_predictions.py
```

**Simple output**

- Please find the example_three_run file to have an idea about performance output.

**Acknowledgements**

- We express our heartfelt gratitude to the creators of all the software components that constitute the PaSSer [5] pipeline.
- The authors of PaSSer [5] are acknowledged for their generous contribution of data.


**References**


[1] Le Guilloux, Vincent, Peter Schmidtke, and Pierre Tuffery. "Fpocket: an open source platform for ligand pocket detection." BMC bioinformatics 10.1 (2009): 1-11.

[2] Chapman, Brad, and Jeffrey Chang. "Biopython: Python tools for computational biology." ACM Sigbio Newsletter 20.2 (2000): 15-19.

[3] Dong, Jie, et al. "PyBioMed: a python library for various molecular representations of chemicals, proteins and DNAs and their interactions." Journal of cheminformatics 10 (2018): 1-11.

[4] Bonidia, Robson P., et al. "MathFeature: feature extraction package for DNA, RNA and protein sequences based on mathematical descriptors." Briefings in bioinformatics 23.1 (2022): bbab434.

[5] Xiao, Sian, Hao Tian, and Peng Tao. "PASSer2. 0: accurate prediction of protein allosteric sites through automated machine learning." Frontiers in Molecular Biosciences 9 (2022): 879251.

