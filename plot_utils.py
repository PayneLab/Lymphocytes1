from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt
import math
from numpy import isnan

import numpy as np
import pandas as pd

def PCA_by_type(data, cell_types=["B_cells","T_cells"]):
    pca = PCA(n_components=5)

    alist=data.values.flatten()
    alist= [a for a in alist if not isnan(a)]
    nan_appoximate = float(alist[math.ceil(float(len(alist))*.01)])
    pca_result = pca.fit_transform(np.nan_to_num(data.transpose(), nan=nan_appoximate))

    samples=np.array(data.columns.values)

    sns.set_style("white")
    figure = plt.figure()
    for cell_type in cell_types:
        cells_of_type = list(i for i,s in enumerate(samples) if cell_type in s)
        plt.scatter(pca_result[cells_of_type,0],pca_result[cells_of_type,1])

    plt.legend(cell_types, loc='upper left', bbox_to_anchor=(1.05, 1))
    return figure