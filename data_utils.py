from numpy import isnan
from numpy import nan
from numpy import log2
from statistics import median
from scipy.stats import ttest_ind
from statistics import mean
import pandas as pd

def check_present(row):
    bool_row = []
    for i in row:
        bool_row.append(not isnan(i))
    return sum(bool_row)

def check_three_of_each_type(row, cell_types=["1_B_", "1_T_"]):
    present_in_types = []
    for i in cell_types:
        cells_of_type = list(s for s in row.index if i in s)
        data_by_type = row.loc[cells_of_type]
        in_type = check_present(data_by_type)
        three_in_type = bool(in_type >= 3)
        present_in_types.append(three_in_type)
    if sum(present_in_types) == len(cell_types):
        return True
    else: return False
    
def check_presence_absence(row, cell_types=["1_B_", "1_T_"], min_reps=3):
    present_in_types = {}
    for i in cell_types:
        cells_of_type = list(s for s in row.index if i in s)
        data_by_type = row.loc[cells_of_type]
        in_type = check_present(data_by_type)
        present_in_types[i] = in_type
    if 0 in list(present_in_types.values()):#absent in one type
        if present_in_types[cell_types[0]] >= min_reps:
            return "B_cell"
        elif present_in_types[cell_types[1]]>= min_reps:
            return "T_cell"
        
def ttest_wrapper(row, cell_types = ["1_B_", "1_T_"]):
    split_row = []
    for i in cell_types:
        cells_of_type = list(s for s in row.index if i in s)
        split_row.append(row.loc[cells_of_type])
        
    tstat = ttest_ind(split_row[0],split_row[1])
    tstat = pd.Series(dict(statistic=tstat[0], pvalue=tstat[1]))
    return tstat
        
def get_fold_changes(row, cell_types=["1_B_", "1_T_"]):
    means = {}
    for i in cell_types:
        cells_of_type = list(s for s in row.index if i in s)
        data_by_type = row.loc[cells_of_type]
        means[i] = mean(data_by_type)
    return (means[cell_types[0]]-means[cell_types[1]])


#get higher-in-B and higher-in-T proteins
def is_altered(tscore, pvalue=.01,change_factor=2):
    log2_fold_change=log2(change_factor)
    if tscore['pvalue'] < pvalue:
        if tscore['log2(B)-log2(T)'] > log2_fold_change:
            #first type is statistically bigger
            return 'B cells'
        elif tscore['log2(B)-log2(T)'] < -log2_fold_change:
            #second type is statistically bigger
            return "T cells"
        
        
