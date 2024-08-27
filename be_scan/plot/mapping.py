"""
Author: Calvin XiaoYang Hu
Adapted from: Kevin Ngan - v2.2
Date: 240515

{Description: }
"""

import os, re, sys, time, warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import scipy.stats as stats
import scipy.spatial.distance as dist
import scipy.interpolate as interp
import scipy.optimize as sp_opt
import scipy.cluster as sp_cl


# import statsmodels.api as sm
import statsmodels.stats.multitest as smm
from statsmodels.nonparametric.smoothers_lowess import lowess
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from Bio.Data import IUPACData

import warnings
warnings.filterwarnings("ignore")

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

from be_scan.plot.clustering_helper import *

# goal is to input a dataframe of guide, mutations, logFC
# get the pdb id and have a 3D coordinates of the residues

def explode_standard(data, mutation_col): 
    # load the dataframe
    data[mutation_col] = data[mutation_col].str.match(r'^[A-Za-z]{1,3}\d{1,4}[A-Za-z.]{1,3}$')
    # explode the dataframe so that each row is now a mutation instead of a guide
    new_data = data.explode(mutation_col)
    # save the dataframe
    return new_data
1
def mapping(df_filepath, mutation_col, score_col, 
            AlphaFoldID, threshold=1): 
    
    # load the dataframe
    df_filepath = Path(df_filepath)
    df_data = pd.read_csv(df_filepath)
    df_data = explode_standard(df_data, mutation_col)

    # trim the dataframe based on threshold
    df_trimmed = df_data[df_data[score_col] >= threshold] 
    
    # load AlphaFold structure


    # align the trimmed list of mutations onto AlphaFold structure

    # produce a structure with residues highlighted (pymol?)
