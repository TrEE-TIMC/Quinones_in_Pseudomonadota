import os
import glob
import sys
import argparse
import pandas as pd
import numpy as np
from functools import reduce
# from utils.fonctions import fonctions_utiles as fu
import quinones

sys.path.append("..")
sys.path.append("/home/choberts/QUINEVOL/scripts/")
# pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser()
parser.add_argument('-wd', '--working_directory', required=True,
                    help="directory containing the results file")
parser.add_argument('-t', '--tax_level', required=False,
                    help="precise on which taxonomy level to filter")
parser.add_argument('-v', '--value', required=False,
                    help="precise which clade, if tax_level = phylum,\
                    precise Proteobacteria or Cyanobacteria")
args = parser.parse_args()


def remove_trail_slash(path):
    if path.endswith('/'):
        path = path[:-1]
    return path


def read_table(tab):
    if "proteins" in tab:
        dataframe = pd.read_csv(tab, sep="\t")
    else:
        dataframe = pd.read_csv(tab, sep="\t", index_col=0)
    return dataframe


def dataframe_differences(macsyres, macsyprofileres):
    '''
    Some genomes don't have a system (number of genes doesn't
    reach the threshold). So, they are missing in the macsyres
    results but are present in macsyprofile results.
    '''
    macsyprof_only = macsyprofileres[~macsyprofileres['genome'].isin(
        macsyres['genome'])]
    return macsyprof_only


def merged_annotation_table(macsyres, macsyprofileres, macsyprof_only):
    '''
    pour pouvoir obtenir la table complète d'annotation :
        - reporter les résultats pour les hits hors systèmes
        dans la table avec les systèmes
        - méthode nécessite une correspondance parfaite des
        tables -> je dois supprimer les génomes qui n'ont pas
        de systèmes pour les remettre (concat) dans un second temps

    => Avec la méthode ci-dessous > la table contient tous les hits
    Si je veux sélectionner les hits (par ex MenE) par rapport à
    leur place dans un système, il faudra faire autrement.
    '''
    macsyprof = macsyprofileres[~macsyprofileres['genome'].isin(
        macsyprof_only['genome'])]
    macsyprof = macsyprof.reset_index(drop=True)
    mask = macsyprof[quinones.quinones_genes] >= 1
    merged_dataframe = macsyres.where(~mask, macsyprof)
    merged_dataframe = pd.concat([merged_dataframe, macsyprof_only])
    return merged_dataframe


def count_gene_pw(df_pw, pw_name):
    """
    insert in the given dataframe column with the count of genes
    for the given pathway
    """
    for index, row in df_pw.iterrows():
        c = 0
        for values in quinones.quinones_pw[pw_name]:
            if row[values] > 0:
                c += 1
            df_pw.at[index, pw_name] = int(c)
    return df_pw


def check_pw(df, pw_name):
    """
    Gives the information if the condition is fulfilled to
    validate a pathway (pw_name provided) if the pw is validated,
    in the column pw_name_found the value is set to 1
    (ex : MK_pw_found : 1) if not, the value is set to 0
    """
    for index, row in df.iterrows():
        if pw_name == "UQ-O2dep":
            if row["UQ-com_pw"] >= 5 and row["UQ-O2dep_pw"] >= 1:
                df.at[index, pw_name+'_found'] = 1
            else:
                df.at[index, pw_name+'_found'] = 0
        if pw_name == "UQ-O2indep":
            if row["UQ-com_pw"] >= 5 and row["UQ-O2indep_pw"] >= 3:
                df.at[index, pw_name+'_found'] = 1
            else:
                df.at[index, pw_name+'_found'] = 0
        if pw_name == "MK-classical":
            if row["MK-classical_pw"] >= 6:
                df.at[index, pw_name+'_found'] = 1
            else:
                df.at[index, pw_name+'_found'] = 0
        if pw_name == "MK-futalosine":
            if row["MK-futalosine_pw"] >= 6:
                df.at[index, pw_name+'_found'] = 1
            else:
                df.at[index, pw_name+'_found'] = 0
        if pw_name == "RQ":
            if row["RQ_pw"] >= 1 and (row["UQ-O2indep_pw"] == 3
                                      or row["UQ-O2dep_pw"] >= 1):
                df.at[index, pw_name+'_found'] = 1
            else:
                df.at[index, pw_name+'_found'] = 0
        if pw_name == "PQ":
            if row["PQ_pw"] >= 6:
                df.at[index, pw_name+'_found'] = 1
            else:
                df.at[index, pw_name+'_found'] = 0
    return df


def q_pathways(dataframe):
    for pw in quinones.quinones_pw.keys():
        count_gene_pw(dataframe, pw)
    check_pw(dataframe, "UQ-O2dep")
    check_pw(dataframe, "UQ-O2indep")
    check_pw(dataframe, "MK-classical")
    check_pw(dataframe, "MK-futalosine")
    check_pw(dataframe, "RQ")
    check_pw(dataframe, "PQ")


directory = remove_trail_slash(args.working_directory)
print(directory)
macsy_tab = glob.glob(directory+'/macsyres*tax*')[0]
macsyprofile_tab = glob.glob(directory+'/macsyprofileres*tax*')[0]
prot_pos_tab = glob.glob(directory+'/*proteins_pos*')[0]


macsyres = read_table(macsy_tab)
macsyprofile = read_table(macsyprofile_tab)
macsyprof_only = dataframe_differences(macsyres, macsyprofile)
merged_df = merged_annotation_table(macsyres, macsyprofile, macsyprof_only)


if args.value:
    if args.value == 'Proteobacteria':
        merged_df = merged_df[merged_df['class'] != 'Epsilonproteobacteria']

if args.value:
    if args.value == 'Proteobacteria':
        merged_df = merged_df[merged_df['phylum'] != 'Proteobacteria']


q_pathways(merged_df)

prot_pos = read_table(prot_pos_tab)

dataset = directory.strip("data/")

# os.mkdir('results/'+dataset)

merged_df.to_csv('results/'+dataset+'/annotation_results_tax.tsv', sep='\t')
