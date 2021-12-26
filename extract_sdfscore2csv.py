#!/usr/bin/env python
from rdkit.Chem import PandasTools
from rdkit import RDConfig
import pandas as pd
import numpy as np
from sklearn import metrics
import os,sys,string


if len(sys.argv) != 5:
    print("")
    print("usage:  ")
    print(sys.argv[0]," <input sdf file>  <output csv file> <Score header> <agg method for score: min or max>")
    print("For example:")
    print(sys.argv[0],"active_fred_dock.sdf active_fred_score.csv min")
    print("Any question, Please feel free to contact me. info@molcalx.com")
    sys.exit()


sdf_file = sys.argv[1]
if not os.path.exists(sdf_file):
   #message = "Sorry, I cannot find the "%s" file."
   print("Sorry, Cannot find the %s file" % sdf_file)
   sys.exit()

df = PandasTools.LoadSDF(sdf_file,smilesName='SMILES',molColName='Molecule',
           includeFingerprints=True)

ofile = sys.argv[2]

score_head = sys.argv[3]

agg_method = sys.argv[4]

df = PandasTools.LoadSDF(sdf_file,smilesName='SMILES',molColName='Molecule', includeFingerprints=True)

df[score_head] = df[score_head].astype(float)
df_grouped = df.groupby('Title')
vs_score = df_grouped.agg({
    'Title':'first',
    'label':'first',
    score_head : agg_method,
})

vs_score.to_csv(ofile,index=False)
