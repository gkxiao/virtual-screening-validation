#!/usr/bin/env python
# coding: utf-8

import rdkit
from rdkit import Chem
import os,sys,string,argparse
from optparse import OptionParser

parser = argparse.ArgumentParser(description="Add Tile and label to Ledock sdf result.\n")
parser.add_argument('input',metavar='<input>',help="input sdf file")
parser.add_argument('output',metavar='<output>',help="output sdf file")
parser.add_argument('class_label',metavar='<label>',help="label should be active or decoy")
args = parser.parse_args()
ifile = args.input
ofile = args.output
class_label = args.class_label

suppl = Chem.SDMolSupplier(ifile, removeHs=False)
ofile = Chem.SDWriter(ofile)


for mol in suppl:
   if mol is not None:
      title = mol.GetProp('_Name')
      label = class_label
      mol.SetProp('Title',str(title))
      mol.SetProp('label',str(label))
      ofile.write(mol)

