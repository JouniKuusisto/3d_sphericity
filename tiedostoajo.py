# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 13:43:31 2021

@author: Jouni Kuusisto
"""

from sphericity import *
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt

import pandas as pd

mypath = "source_data/ES/analyysiin/"
files = [f for f in listdir(mypath) if isfile(join(mypath, f))]

df = pd.DataFrame(index = files, columns=("Case", "Volume", "LASP", "3DS"))
df["Case"] = [no + 1 for no in range(len(files))]

for file in files:
    print(file)
    case_mesh = loadMesh(join(mypath, file))
    volume, bisbal, kuusisto = printResults(case_mesh) # return volume, bisbal, kuusisto
    print("\n \n")
    plotTriangleCentersAndCenterOfMass(case_mesh, " Mesh volume (ml):" + " {:.2f}".format(volume/1000) + ", LASP: " + "{:.1f}".format(bisbal*100) + "%" + ", 3DS: " + "{:.1f}".format(kuusisto*100) + "%", save = True, save_name= "output2/" + str(df.at[file, "Case"]) + ".png")
    df.at[file, "Volume"] = volume/1000
    df.at[file, "LASP"] = bisbal
    df.at[file, "3DS"] = kuusisto
    
def blandAltman(files):
    dfBA = pd.DataFrame(index = files, columns = ("name", "Volume", "LASP", "3DS"))
    dfBA["3DS"] = df["3DS"]
    dfBA["LASP"] = df["LASP"]
    dfBA["3DS-LASP"] = df["3DS"] - df["LASP"]
    dfBA["mean"] = (df["3DS"] + df["LASP"]) / 2
    ba_bias = dfBA["3DS-LASP"].mean(axis=0)*100
    min_range = dfBA["mean"].min()*100
    max_range = dfBA["mean"].max()*100
    std = dfBA["3DS-LASP"].std()*100
    loa_upper = ba_bias + 2*std
    loa_lower = ba_bias - 2*std
    plt.plot([min_range, max_range], [ba_bias, ba_bias], c='r') # bias
    plt.plot([min_range, max_range], [loa_upper, loa_upper], c='y') # 2SD
    plt.plot([min_range, max_range], [loa_lower, loa_lower], c='y') # -2SD
    plt.scatter(dfBA["mean"]*100, dfBA["3DS-LASP"]*100)
    plt.ylim((-10,10))
    plt.yticks([tick for tick in range(-10, 12, 2)])
    plt.xlabel("(3DS+LASP) / 2 (%)")
    plt.ylabel("3DS-LASP (% points)")
    plt.title("Bland-Altman Plot of 3DS and LASP ")
    plt.savefig("BA-sphericity.png", dpi=300)
    plt.show()
    
    return ba_bias, loa_upper, loa_lower
    
bias, loa_upper, loa_lower = blandAltman(files)
    