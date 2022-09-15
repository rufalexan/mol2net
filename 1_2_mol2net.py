#python3













'''
================================
1_2_mol2net -- MOLECULAR NETWORK
================================

AUTHOR: Alexander Ruf (rufalexan@gmail.com)

If you use this code for your work, please cite the corresponding paper:

10.1021/acs.analchem.2c01271

, this github repository and/or the corresponding Zenodo DOI:

@software{rufalexan_2022_7025094,
  author       = {rufalexan},
  title        = {rufalexan/mol2net: v0.1.0},
  month        = aug,
  year         = 2022,
  publisher    = {Zenodo},
  version      = {v0.1.0},
  doi          = {10.5281/zenodo.7025094},
  url          = {https://doi.org/10.5281/zenodo.7025094}
}

'''




















from numpy import *
import pandas as pd
import string
from collections import namedtuple
from scipy import stats
from itertools import combinations 
import matplotlib as mpl
#from matplotlib import cm
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import datetime

begin_time = datetime.datetime.now()



################################
##                          DATA
################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##name = 'pubchem_database_CHNO_1-10000'
name = 'TEST_pubchem_database_CHNO_10001-20000'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
matrix = pd.read_csv(name + '.tsv', sep='\t', skiprows=0).to_numpy()
##name = 'PubChem 1-10000'              #################################
name = 'TEST_PubChem 10001-20000'              #################################

id = around(matrix[:,0].reshape(len(matrix),1).astype(float),6)
mz = around(matrix[:,2].reshape(len(matrix),1).astype(float),6)
mz2 = mz
intensity = around(matrix[:,3].reshape(len(matrix),1).astype(float),6)
data = concatenate([mz, id, intensity], axis=1)
data2 = data



################################
##                    list match
################################
name_reactions = 'TRAFOS'
reactions_input = pd.read_csv(name_reactions + '.tsv', sep='\t', skiprows=0).to_numpy()
h = reactions_input[:,0]
c = reactions_input[:,1]
n = reactions_input[:,2]
o = reactions_input[:,3]
s = reactions_input[:,4]
m_h = 1.007825
m_c = 12.
m_o = 15.994915
m_n = 14.003074
m_s = 31.972071
masses = (m_h, m_c, m_n, m_o, m_s)
elements = transpose(array([h, c, n, o, s]))
mass_differences = around(array(abs(dot(elements, masses)).reshape(len(reactions_input),1), dtype=float),6)

def prepend(list, str): 
    str  += '{0}'
    list = [str.format(i) for i in list] 
    return(list)

hydrogen = 'H'
carbon = 'C'
nitrogen = 'N'
oxygen = 'O'
sulfur = 'S'

hydrogens = prepend(h, hydrogen) 
carbons = prepend(c, carbon)
nitrogens = prepend(n, nitrogen)
oxygens = prepend(o, oxygen)
sulfurs = prepend(s, sulfur)

formulas = [str(carbons[i]) + str(hydrogens[i]) + str(nitrogens[i]) + str(oxygens[i]) + str(sulfurs[i]) + str(' ') for i in range(len(carbons))]

formulas = [sub.replace('C1H', 'CH') for sub in formulas] 
formulas = [sub.replace('H1N', 'HN') for sub in formulas] 
formulas = [sub.replace('N1O', 'NO') for sub in formulas] 
formulas = [sub.replace('C0', '') for sub in formulas] 
formulas = [sub.replace('H0', '') for sub in formulas] 
formulas = [sub.replace('N0', '') for sub in formulas] 
formulas = [sub.replace('O0', '') for sub in formulas] 
formulas = [sub.replace('O1', 'O') for sub in formulas]
formulas = [sub.replace('S0', '') for sub in formulas]
formulas = [sub.replace('S1', 'S') for sub in formulas]
formulas = array(formulas).reshape(len(hydrogens),1)

reactions = concatenate([elements, mass_differences, formulas], axis=1)

header_reactions = (array(['H', 'C', 'N', 'O', 'S', 'mass_difference', 'molecular_formula'])).reshape(len(reactions.T),1).T
reactions2 = concatenate((header_reactions, reactions))

df_reactions = pd.DataFrame({'mass_difference': reactions[:, 5], 'molecular_formula': reactions[:, 6], 'H': reactions[:, 0], 'C': reactions[:, 1], 'N': reactions[:, 2], 'O': reactions[:, 3], 'S': reactions[:, 4]}).reindex(columns=(['molecular_formula', 'mass_difference', 'H', 'C', 'N', 'O', 'S']))









################################
##                    md matches
################################
md_matches = namedtuple('md_matches', 'md_matches hits')
new_list = []
for md in mass_differences:
	for element in mz:
		if element+md in mz2:
			tmp_new_list = md_matches(element, element+md)
			new_list.append(tmp_new_list)

matches = array(new_list).reshape(len(new_list),2)
matches = matches[argsort(matches[:,0])]
df_sources = pd.DataFrame({'mz_sources': data[:, 0], 'id_sources': data[:, 1]})
df_targets = pd.DataFrame({'mz_targets': data[:, 0], 'id_targets': data[:, 1]})
df_matches = pd.DataFrame({'mz_sources': matches[:, 0], 'mz_targets': matches[:, 1]})
df_source_match = df_matches.merge(df_sources, how='left', on=['mz_sources'])
df_target_match = df_source_match.merge(df_targets, how='left', on=['mz_targets'])
mass_diff = around(df_target_match['mz_targets'] - df_target_match['mz_sources'],6)
type = full((len(mass_diff),1), 'Undirected')
type = pd.DataFrame({'type': type[:, 0]})
Label = full((len(mass_diff),1), 'x')
Label = pd.DataFrame({'Label': Label[:, 0]})
Weight = full((len(mass_diff),1), 1)
Weight = pd.DataFrame({'Weight': Weight[:, 0]})

df_full_edge_matrix = pd.concat([df_target_match, type, Label, Weight, mass_diff.rename('mass_difference')], axis=1)
edge_list = ( df_full_edge_matrix.reindex(columns=['id_sources', 'id_targets', 'type', 'Weight', 'mz_sources', 'mz_targets', 'Label', 'mass_difference']) ).rename(columns={'id_sources': 'Source', 'id_targets': 'Target'})
df_mz_md = (pd.DataFrame({'mass_difference': reactions[:, 5], 'molecular_formula': reactions[:, 6]})).reindex(columns=['molecular_formula', 'mass_difference'])
df_mz_md['mass_difference'] = df_mz_md['mass_difference'].astype(str)
edge_list['mass_difference'] = edge_list['mass_difference'].astype(str)
df_edge_list = edge_list.merge(df_mz_md, how='left', on=['mass_difference'])

df_edge_list['Source'] = df_edge_list['Source'].astype(int)
df_edge_list['Target'] = df_edge_list['Target'].astype(int)







################################
##      freq of mass differences
################################
freq_of_md = stats.itemfreq(mass_diff)
df_freq_of_md = pd.DataFrame({'mass_difference': freq_of_md[:, 0], 'freq_of_md': freq_of_md[:, 1]})
df_freq_of_md = df_freq_of_md.reindex(columns=['mass_difference', 'freq_of_md'])

df_freq_of_md['mass_difference'] = df_freq_of_md['mass_difference'].astype(str)
df_reactions['mass_difference'] = df_reactions['mass_difference'].astype(str)
df_reactions_freq = df_reactions.merge(df_freq_of_md, how='left', on=['mass_difference'])
#df_reactions_freq = pd.concat([df_reactions, df_freq_of_md], axis=1, sort=False)

df_reactions_freq.plot.bar(x='molecular_formula', y='freq_of_md')







################################
##                        export
################################
df_edge_list_gephi = df_edge_list

df_edge_list_gephi2 = df_edge_list_gephi[['Source', 'Target', 'type', 'mass_difference', 'mz_sources', 'mz_targets', 'Label', 'molecular_formula']].reindex(columns=['Source', 'Target', 'type', 'mass_difference', 'mz_sources', 'mz_targets', 'Label', 'molecular_formula']).rename(columns={'mass_difference':'Weight', 'mz_sources': 'mz_source', 'mz_targets':'mz_target'})       # Weights = MDs

formulas_list = formulas.tolist()

df_edge_list_gephi3 = df_edge_list_gephi.merge(df_freq_of_md, how='left', on=['mass_difference'])
df_edge_list_gephi3.to_csv('EDGES_' + name + str(formulas_list) + '.tsv', sep='\t', index=False)    # with frequencies of MDs






print(name)
print(df_reactions_freq)













print('Time code = '+str(datetime.datetime.now() - begin_time)+' [h:min:sec]')







