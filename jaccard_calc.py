import numpy as np
import pandas as pd
import os
from sklearn.metrics import jaccard_score
import operator

########## read in GMT file and split into list of names and gene_lists

raw_library_data = []

# put path to GMT file here 
with open('/Users/jessiecheng/Documents/maayanlab/mygenesets_down.gmt.txt') as f:
    for line in f.readlines():
        raw_library_data.append(line.split("\t\t"))

name = []
gene_list = []

for i in range(len(raw_library_data)):
    raw_row = raw_library_data[i][0].split("\t")
    j = 0 
    while j < len(raw_row):
        if (raw_row[j] == "\n"):
            raw_row.pop(j)
        else:
            a = raw_row[j].strip()
            if (a == ""):
                raw_row.pop(j)
            else:
                raw_row[j] = a
                j = j+1
    name += [raw_row[0]]
    gene_list += [raw_row[2:]]

# optional - print to see what name and gene_list look like
print(name)
print(gene_list)
print(len(name))
print(len(gene_list))

########## calculate jaccard distances and store matrix in a new csv file 

def jaccard_calc(setA, setB):
    num = 0
    denom = len(setA) + len(setB)
    for s in setA:
        if (s in setB or s.lower() in setB):
            num = num+1
            denom = denom-1
    return num/denom

cols = []
for i in range(len(name)):
    col = [0] * (len(name)+1)
    fields = name[i].split(",")
    col[0] = fields[0] + ", " + fields[3]  ## you can choose which fields to include in set names
    for j in range(len(name)):
        col[j+1] = jaccard_calc(gene_list[i], gene_list[j])
    cols.append(col)

name = ['Name'] + name
df = pd.DataFrame(np.array(cols), columns=name)

# optional - print to see what the matrix looks like
print(df)

# exports jaccard matrix as a csv file, rename file as you wish
df.to_csv('mygenesets_down_jaccard.csv')







