#python3











'''
==============================
3_netAna  --  NETWORK ANALYSIS
==============================

AUTHOR: Alexander Ruf (rufalexan@gmail.com)

If you use this code for your work, please cite the corresponding paper, this github repository and/or the corresponding Zenodo DOI.
'''







import datetime
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import cm
from networkx.drawing.nx_agraph import graphviz_layout
import pydot
import matplotlib as mpl
import matplotlib.ticker as ticker
begin_time = datetime.datetime.now()
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#name    = 'PubChem 1-10000[H2 ], [O ], [CO ], [H3N ]'
#name2   = 'pubchem_database_CHNO_1_10000'
#name    = 'PubChem 10001-20000[H2 ], [O ], [CO ], [H3N ]'
#name2   = 'pubchem_database_CHNO_10001_20000'
#name    = 'Ice residue[H2 ], [O ], [CO ], [H3N ]'
name    = 'EDGES_[H2 ], [O ], [CO ], [H3N ]_min_set'
#name     = 'EDGES_[H2 ], [O ], [CO ], [H3N ], [CHN ], [CH3N ], [H2O ], [CO2 ], [CH2O ], [CHNO ]'
#name    = 'Ice residue[O ], [CO ], [H3N ]'             # min_set - H2
#name    = 'Ice residue[H2 ], [CO ], [H3N ]'            # min_set - O
#name    = 'Ice residue[H2 ], [O ], [H3N ]'             # min_set - CO
#name    = 'Ice residue[H2 ], [O ], [CO ]'              # min_set - NH3
name2   = 'NODES_aurelien_matrix'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_edges = pd.read_csv(name + '.tsv', sep='\t', header=0)
df_nodes = pd.read_csv(name2 + '.tsv', sep='\t', header=0)
G = nx.from_pandas_edgelist(df_edges, 'Source', 'Target', create_using=nx.Graph())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#name = 'PubChem subset 1-10000'
#name = 'PubChem subset 10001-20000'
#name = 'Ice residue'
name2 = 'H$_2$O:CH$_3$OH:NH$_3$'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####################
#         COMPONENTS
####################
Gcc = list(nx.connected_components(G))
d = {name:k for k,comp in enumerate(Gcc) for name in comp}               #dict(enumerate(Gcc))
df_comp = pd.DataFrame.from_dict(d, orient='index', columns=['component']).rename_axis('id').reset_index()
df_nodes2 = df_nodes.merge(df_comp, how='left', on='id')
#nans = df_nodes2[df_nodes2.isna().any(axis=1)]
df_nodes_cc = df_nodes2[~df_nodes2['component'].isnull()]
df_source = df_edges.rename(columns={'Source':'id'})
df_edges2 = df_source.merge(df_comp, how='left', on='id')
df_edges2 = df_edges2.rename(columns={'component':'component source'})
df_edges2 = df_edges2.rename(columns={'id':'Source'})
df_target = df_edges2.rename(columns={'Target':'id'})
df_edges3 = df_target.merge(df_comp, how='left', on='id')
df_edges3 = df_edges3.rename(columns={'component':'component target'})
df_edges3 = df_edges3.rename(columns={'id':'Target'})
df_edges3 = df_edges3.rename(columns={'mass_difference':'mass difference'})
df_edges3 = df_edges3.rename(columns={'molecular_formula':'molecular formula'})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FILTER_COMPONENTS  =   [0,1]  # [3,12,35,1,17]         #     np.inf         list(range(0,10+1))          PubChem 1-10000 [24,2]      PubChem 10001 [24,3]       min_set - H2   [7,17,34,38,6]           min_set - O    [13,5,8,25,15]            min_set - CO   [18,6,37,50,34]            min_set - NH3  [3,12,35,1,17]

#df_edges_cc = df_edges3[(df_edges3['component source'] <= FILTER_COMPONENTS)]     
#df_nodes_cc = df_nodes_cc[(df_nodes_cc['component'] <= FILTER_COMPONENTS)]
#FILTER_COMPONENTS = '0-'+ str(FILTER_COMPONENTS)

#df_edges_cc = df_edges3[(df_edges3['component source'] == FILTER_COMPONENTS)]   # das brauche ich wegen list/isin nicht mehr?
#df_nodes_cc = df_nodes_cc[(df_nodes_cc['component'] == FILTER_COMPONENTS)]

df_edges_cc = df_edges3[(df_edges3['component source'].isin(FILTER_COMPONENTS))]
df_nodes_cc = df_nodes_cc[(df_nodes_cc['component'].isin(FILTER_COMPONENTS))]

#df_nodes_cc['lg_abundance_int'] = np.log10(df_nodes_cc['abundance_int'])           # because range of intensities is too large for colorbar legend
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

G_cc = nx.from_pandas_edgelist(df_edges_cc, 'Source', 'Target', create_using=nx.Graph())

####################
#             COLORS
####################
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                 PUBCHEM
#color_type = 'Compound class'      #    H  C  N  O         H/C  N/C  O/C        exact_mass  DBE      Compound class    Compound class_aromatics      X_C  complexity  logP   component         ### "/" filenames cannot be saved !

#########sample_colors = pd.DataFrame({color_type: ['Other',     'CH_aromatics', 'Acids',  'Alcohols', 'Ketones',  'Esters', 'Amides', 'Nitriles', 'Amines', 'N_cyclic_aliphatics', 'N_aromatics'],
#########                        'sample_color':   ['whitesmoke', 'black',       'lightseagreen','blue',     'skyblue',  'cyan',   'gold',   'lime',     'green',  'magenta',              'red'] })         

#sample_colors = pd.DataFrame({color_type: ['CH_aromatics', 'Acids',  'Alcohols', 'Ketones',  'Esters', 'Amides', 'Nitriles', 'Amines', 'N_cyclic_aliphatics', 'N_aromatics'],
#                        'sample_color':   ['black',       'lightseagreen','blue',     'skyblue',  'cyan',   'gold',   'lime',     'green',  'magenta',              'red'] })         
##########                                   CH            |----             CHO               ---|  |- CHNO-||---      CHN     ---||---         N heterocycles        ---|

##sample_colors = pd.DataFrame({color_type: ['Other',      'CH_aromatics', 'N_aromatics', 'CH_aromatics_N_aromatics'],
##                        'sample_color':   ['whitesmoke', 'black',        'red',         'blue'] })         
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                     ICE_RESIDUE
color_type = 'Ice composition H2O_CH3OH_NH3'                 # Ice composition H2O_CH3OH_NH3    H  C  N  O  H/C  N/C  O/C  mz_theo_neutral  abundance_int - lg_abundance_int  mu  DBE  N/O  O/H         component         ### "/" filenames cannot be saved !

sample_colors = pd.DataFrame({color_type: ['3_1_0.2', '3_1_1',   '3_1_1 16h', '3_1_5', '10_1_1',    '3_1_1 overirradiation'], 
                        'sample_color':   ['blue',    '#15B01A', 'green',     'red',   'lightgray', 'orange']})
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

carac = df_nodes_cc[['id', color_type]]
carac = carac.set_index('id')
carac = carac.reindex(G_cc.nodes())
carac[color_type]=pd.Categorical(carac[color_type])
carac[color_type].cat.codes
##  DISCRETE, MANUAL COLORS (e.g. colored by samples)
sample_names = df_nodes_cc[color_type].unique()                       ###### no CH_aromatics dabei, und dann nan dabei
carac2 = carac.merge(sample_colors, how='left', on=color_type)            
node_colors = carac2['sample_color'].tolist()























############
###   LABELS
############
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##node_label = 'molecular_formula'
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##node_labels = df_nodes_cc[['id', node_label]]
##df_source = df_edges_cc.rename(columns={'Source':'id'})
##df_edges2 = df_source.merge(node_labels, how='left', on='id')
##df_edges2 = df_edges2.rename(columns={'id':'Source'})
##df_edges2 = df_edges2.rename(columns={node_label:node_label+'_source'})
##df_target = df_edges2.rename(columns={'Target':'id'})
##df_edges2 = df_target.merge(node_labels, how='left', on='id')
##df_edges2 = df_edges2.rename(columns={'id':'Target'})
##df_edges2 = df_edges2.rename(columns={node_label:node_label+'_target'})
###df_edges2 = df_edges2.rename(columns={'molecular_formula_x':node_label})
###df_edges2[node_label+'_source'].astype(str)
###df_edges2[node_label+'_target'].astype(str)

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##edge_label = 'molecular formula'
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##df_edges2 = df_edges2.rename(columns={edge_label:'Edge_label'})
##G_cc = nx.from_pandas_edgelist(df_edges2, node_label+'_source', node_label+'_target', edge_attr='Edge_label', create_using=nx.Graph())
##edge_labels = nx.get_edge_attributes(G_cc, "Edge_label")

##label_size=30         #12
##node_sizes=2000       #400

##plt.figure()
##nx.draw_networkx_edge_labels(G_cc, pos = nx.nx_pydot.graphviz_layout(G_cc), edge_labels=edge_labels, font_size=label_size)
##nx.draw(G_cc, pos = nx.nx_pydot.graphviz_layout(G_cc), with_labels=True, node_color=node_colors, node_size=node_sizes, width=1, font_size=label_size, alpha=.5)
###nx.draw(G_cc, pos = nx.nx_pydot.graphviz_layout(G_cc), with_labels=True, node_color=carac[color_type].cat.codes, cmap=plt.cm.jet, node_size=node_size, width=.3, font_size=label_size, alpha=.5)     

##plt.show()

























################################
###          DEGREE DISTRIBUTION
################################
##m=0
##degree_freq = nx.degree_histogram(G_cc)
##degrees = range(len(degree_freq))

####plt.figure(figsize=(6, 5)) 
#####plt.loglog(degrees[m:], degree_freq[m:],'ko-') 
####plt.plot(degrees[m:], degree_freq[m:],'ko-', label='') 
####plt.axvline(x=np.mean(degrees[m:]), color='k', linestyle='--', alpha=.3)
####plt.axhline(y=np.mean(degree_freq[m:]), color='k', linestyle='--', alpha=.3)
####plt.xlim(-.5,26.5)
####plt.xlabel('Degree')
####plt.ylabel('Frequency')
#####plt.legend()
####plt.title(name+', Component '+str(FILTER_COMPONENTS), wrap=True)
######plt.savefig(name+', Comp'+str(FILTER_COMPONENTS)+'    degree distribution'+'.png')


##G = G_cc

##degrees = [val for (node, val) in G.degree()]

##df_nodes_cc['degrees'] = degrees


###df_nodes_cc.to_csv('__________________.txt', header="'\t'.join(list(df_export.columns))", index=None, sep='\t')
###df_edges_cc.to_csv('__________________edges.txt', header="'\t'.join(list(df_export.columns))", index=None, sep='\t')

















#df_nodes_cc[(df_nodes_cc['id'] == 122)]


#df_nodes_cc.plot.scatter('mz_theo_neutral', 'lg_abundance_int')
####df_nodes_cc['abundance_int'].value_counts().plot()
#####df_nodes_cc['abundance_int'].cumsum().plot()
#plt.show()






#######df = pd.DataFrame(df_edges_cc['Source'].unique())
#######coverage_nodes = int(round(df.shape[0] / df_nodes_cc.shape[0] * 100, 0))     # maybe not correct as some nodes can be also uniquely as targets
#coverage_nodes = int(round(df_nodes_cc.shape[0] / df_nodes.shape[0] * 100, 0))
#print(name2+', '+name+'  --  Component '+str(FILTER_COMPONENTS))
#print('Coverage_nodes = '+str(coverage_nodes)+' %')










############################
#####         AVERAGE VALUES
############################
#C_min = df_nodes_cc['C'].min().astype(int)
#C_max = df_nodes_cc['C'].max().astype(int)
#H_min = df_nodes_cc['H'].min().astype(int)
#H_max = df_nodes_cc['H'].max().astype(int)
#N_min = df_nodes_cc['N'].min().astype(int)
#N_max = df_nodes_cc['N'].max().astype(int)
#O_min = df_nodes_cc['O'].min().astype(int)
#O_max = df_nodes_cc['O'].max().astype(int)
#mz_min = df_nodes_cc['mz_theo_neutral'].min().astype(int)
#mz_max = df_nodes_cc['mz_theo_neutral'].max().astype(int)
#DBE_min = df_nodes_cc['DBE'].min().astype(int)
#DBE_max = df_nodes_cc['DBE'].max().astype(int)

#C_median = df_nodes_cc['C'].median().astype(int)
#H_median = df_nodes_cc['H'].median().astype(int)
#N_median = df_nodes_cc['N'].median().astype(int)
#O_median = df_nodes_cc['O'].median().astype(int)
#mz_median = df_nodes_cc['mz_theo_neutral'].median().astype(int)
#DBE_median = df_nodes_cc['DBE'].median().astype(int)

##print('C_min = '+str(df_nodes_cc['C'].min()))
##print('C_max = '+str(df_nodes_cc['C'].max()))
##print('H_min = '+str(df_nodes_cc['H'].min()))
##print('H_max = '+str(df_nodes_cc['H'].max()))
##print('N_min = '+str(df_nodes_cc['N'].min()))
##print('N_max = '+str(df_nodes_cc['N'].max()))
##print('O_min = '+str(df_nodes_cc['O'].min()))
##print('O_max = '+str(df_nodes_cc['O'].max()))
##print('mz_min = '+str(df_nodes_cc['mz_theo_neutral'].min()))
##print('mz_max = '+str(df_nodes_cc['mz_theo_neutral'].max()))
##print('DBE_min = '+str(df_nodes_cc['DBE'].min()))
##print('DBE_max = '+str(df_nodes_cc['DBE'].max()))

##PERCENT = int(df_nodes_cc.shape[0] / df_nodes.shape[0] * 100)
#PERCENT = np.round(df_nodes_cc.shape[0] / df_nodes.shape[0] * 100, 1)

#print('%'+name+'  --  Component '+str(FILTER_COMPONENTS))
##print('C$_{'+str(C_min)+'-'+str(C_max)+'}$H$_{'+str(H_min)+'-'+str(H_max)+'}$N$_{'+str(N_min)+'-'+str(N_max)+'}$O$_{'+str(O_min)+'-'+str(O_max)+'}$, '+'m/z = '+str(mz_min)+'-'+str(mz_max)+' amu, '+'DBE = '+str(DBE_min)+'-'+str(DBE_max))
##print('C$_{'+str(C_median)+'}$H$_{'+str(H_median)+'}$N$_{'+str(N_median)+'}$O$_{'+str(O_median)+'}$, \\\ '+'m/z = '+str(mz_median)+' amu, '+'DBE = '+str(DBE_median)+' \\\ (median values) \\\ '+str(df_nodes_cc.shape[0])+' molecular formulas ('+str(PERCENT)+' % of all formulas) \\\ ')
##print('C$_{'+str(C_median)+'}$H$_{'+str(H_median)+'}$N$_{'+str(N_median)+'}$O$_{'+str(O_median)+'}$, \\\ '+'m/z = '+str(mz_median)+' amu, '+'DBE = '+str(DBE_median)+' \\\ (median values) \\\ '+str(df_nodes_cc.shape[0])+' molecular formulas \\\ ')
#print('C$_{'+str(C_median)+'}$H$_{'+str(H_median)+'}$N$_{'+str(N_median)+'}$O$_{'+str(O_median)+'}$, \\\ '+'DBE = '+str(DBE_median)+' \\vspace{.5cm} ')













################################
###          ELEMENT HISTOGRAMS
################################
#####~~~~~~~~~~~~~~
#ELEMENT = 'H'
#####~~~~~~~~~~~~~~
#df_element_freq = pd.DataFrame(df_nodes_cc[ELEMENT].value_counts())

#plt.figure()
#pd.DataFrame(df_nodes_cc[ELEMENT].value_counts(normalize=False)).sort_index().plot(kind='bar', rot=0, ylabel='', legend=False, color='black', width=.5)
#plt.title('Component '+str(FILTER_COMPONENTS))          
#plt.xticks(())
#plt.yticks(())
##plt.yscale('log')
#plt.xlabel('#'+str(ELEMENT), fontsize=25)
##plt.ylabel('Frequency (log scale, normalized)')
#plt.ylabel('Frequency', fontsize=25)
#plt.savefig(name+', Comp'+str(FILTER_COMPONENTS)+ '    element_components_'+str(ELEMENT)+'.png')













#################################
####          COMPONENT HISTOGRAM
#################################
#df_comp_freq = pd.DataFrame(df_nodes_cc['component'].value_counts()).rename(columns={'component': 'freq_component'})

##plt.figure()
##pd.DataFrame(df_nodes_cc['component'].value_counts(normalize=True)).sort_index().plot(kind='bar', rot=0, ylabel='', legend=False, color='black', width=1)
###pd.DataFrame(df_nodes_cc['component'].value_counts(normalize=False)).sort_index().plot(kind='bar', rot=0, ylabel='', legend=False, color='black', width=1)
##plt.title('Component '+str(FILTER_COMPONENTS))          
##plt.xticks(np.arange(0, len(df_comp_freq), 10))
##plt.yticks()
##plt.yscale('log')
##plt.xlabel('Component number')
##plt.ylabel('Frequency (log scale, normalized)')
###plt.ylabel('Frequency (log scale)')
##plt.savefig(name+', Comp'+str(FILTER_COMPONENTS)+ '    hist_components'+'.png')

#print(name)
#print(df_comp_freq)





















####################
##     NORMAL LAYOUT
####################
#plt.figure()
#plt.title(name, wrap=True, fontsize=20)
##nx.draw(G_cc, pos = nx.nx_pydot.graphviz_layout(G_cc), with_labels=False, node_color=node_colors, node_size=20, width=.03, font_size=5)
#nx.draw(G_cc, pos = nx.nx_pydot.graphviz_layout(G_cc), with_labels=False, node_color=node_colors, node_size=20, width=.1)
#plt.plot([], [], sample_colors.values[0][1], marker='o', markersize=10, label=sample_colors.values[0][0])     # Legend
#plt.plot([], [], sample_colors.values[1][1], marker='o', markersize=10, label=sample_colors.values[1][0])
#plt.plot([], [], sample_colors.values[2][1], marker='o', markersize=10, label=sample_colors.values[2][0])
#plt.plot([], [], sample_colors.values[3][1], marker='o', markersize=10, label=sample_colors.values[3][0])
#plt.plot([], [], sample_colors.values[4][1], marker='o', markersize=10, label=sample_colors.values[4][0])
#plt.plot([], [], sample_colors.values[5][1], marker='o', markersize=10, label=sample_colors.values[5][0])
##plt.plot([], [], sample_colors.values[6][1], marker='o', markersize=10, label=sample_colors.values[6][0])
##plt.plot([], [], sample_colors.values[7][1], marker='o', markersize=10, label=sample_colors.values[7][0])
##plt.plot([], [], sample_colors.values[8][1], marker='o', markersize=10, label=sample_colors.values[8][0])
##plt.plot([], [], sample_colors.values[9][1], marker='o', markersize=10, label=sample_colors.values[9][0])
##plt.legend()                      #loc='upper left'
##plt.savefig(name+', Comp'+str(FILTER_COMPONENTS)+'    '+color_type+'    network'+'.png')
#plt.show()














####################################################
####    CONTINOUS COLOR GRADIENT (e.g. element maps)
####################################################
#nodes = G_cc.nodes()

#plt.figure()
#plt.title(color_type, wrap=True, fontsize=20)
#nx.draw(G_cc, pos = nx.nx_pydot.graphviz_layout(G_cc), with_labels=False, node_color=carac[color_type].cat.codes, cmap=plt.cm.jet, node_size=20, width=.1, font_size=10)
##pos = nx.spring_layout(G_cc)
##pos = nx.fruchterman_reingold_layout(G_cc)
#pos = nx.nx_pydot.graphviz_layout(G_cc)
#ec = nx.draw_networkx_edges(G_cc, pos, alpha=0.2)
#nc = nx.draw_networkx_nodes(G_cc, pos, nodelist=nodes, node_color=carac[color_type], with_labels=False, node_size=5, cmap=plt.cm.jet) 
#nc = nx.draw_networkx_nodes(G_cc, pos, nodelist=nodes, node_color=carac[color_type], node_size=5, cmap=plt.cm.jet) 
#plt.colorbar(nc)
#plt.axis('off')
##plt.savefig(name+', Comp'+str(FILTER_COMPONENTS)+'    '+color_type+'    network'+'.png')
#plt.show()

































###############################
##     COMPOUND CLASS HISTOGRAM   --  PubChem
###############################
#df_comp_class_freq = pd.DataFrame(df_nodes_cc['Compound class'].value_counts()).rename(columns={'Compound class': 'freq_component'})

#plt.figure()
#pd.DataFrame(df_nodes_cc['Compound class'].value_counts(normalize=True)).sort_index().plot(kind='bar', rot=90, ylabel='', legend=False, color='black', width=1)
##pd.DataFrame(df_nodes_cc['Compound class'].value_counts(normalize=False)).sort_index().plot(kind='bar', rot=90, ylabel='', legend=False, color='black', width=1)
#plt.title(str(name)+', Component '+str(FILTER_COMPONENTS))          
#plt.yticks()
#plt.xlabel('Compound class')
#plt.ylabel('Frequency (normalized)')
##plt.ylabel('Frequency')
#plt.savefig(name+', Comp'+str(FILTER_COMPONENTS)+ '    hist_comp_class'+'.svg')










###############################
##     COMPOUND CLASS PIE CHART   --  PubChem
###############################

#df_nodes_cc['comp_class_counts'] = df_nodes_cc['Compound class'].map(df_nodes_cc['Compound class'].value_counts())
#df = df_nodes_cc[['Compound class', 'comp_class_counts']]
#df = df.drop_duplicates(subset = ['comp_class_counts'])
#sample_colors = pd.DataFrame({color_type: ['CH_aromatics', 'Acids',  'Alcohols', 'Ketones',  'Esters', 'Amides', 'Nitriles', 'Amines', 'N_cyclic_aliphatics', 'N_aromatics'],
#                        'sample_color':   ['black',       'lightseagreen','blue',     'skyblue',  'cyan',   'gold',   'lime',     'green',  'magenta',              'red'] })         
#sample_colors = sample_colors.rename(columns={'DBE':'Compound class'})
#df = df.merge(sample_colors, how='left', on='Compound class')


#plt.figure()
#df.plot(kind='pie', y = 'comp_class_counts', autopct='%1.1f%%', 
# startangle=90, shadow=False, labels=df['Compound class'], legend = False, fontsize=14, colors=df['sample_color'])
#plt.savefig(name+', Comp'+str(FILTER_COMPONENTS)+ '    pie_comp_class'+'.svg')
##plt.show()






























#######################
##  FULL-SCREEN NETWORK
#######################
#plt.figure()
#nx.draw(G_cc, pos = nx.nx_pydot.graphviz_layout(G_cc), with_labels=False, node_color=node_colors, node_size=15, width=.1)
##nx.draw(G_cc, pos = nx.nx_pydot.graphviz_layout(G_cc), with_labels=True, node_color=carac[color_type].cat.codes, cmap=plt.cm.jet, node_size=7, width=.1, font_size=10)
#plt.show()

































##########################
##     CLUSTER TRANSITIONS                              fig 2: min_set     directions between nh3_poor to medium, etc.. 
##         WITH DIRECTIONS
##########################
##~~~~~~~~~~~~~~~~~~~~~~~
#FROM_TO = 'Ice composition H2O_CH3OH_NH3'
##~~~~~~~~~~~~~~~~~~~~~~~
#df_samples = df_nodes_cc[['id', FROM_TO]]
#df_edges4 = df_edges_cc.rename(columns={'Source':'id'})
#df_edges4 = df_edges4.merge(df_samples, how='left', on='id')
#df_edges4 = df_edges4.rename(columns={FROM_TO: FROM_TO+str('_Source'), 'id':'Source'})
#df_edges5 = df_edges4.rename(columns={'Target':'id'})
#df_edges5 = df_edges5.merge(df_samples, how='left', on='id')
#df_edges5 = df_edges5.rename(columns={FROM_TO: FROM_TO+str('_Target'), 'id':'Target'})

##~~~~~~~~~~~~~~~~~~~~~~~
#SAMPLE_START = '3_1_0.2'
#SAMPLE_END   = '3_1_1'

##SAMPLE_START = '3_1_1'
##SAMPLE_END   = '10_1_1'

##SAMPLE_START = '10_1_1'
##SAMPLE_END   = '3_1_5'
##~~~~~~~~~~~~~~~~~~~~~~~
#df_edges6 = df_edges5[(df_edges5[FROM_TO+str('_Source')] == SAMPLE_START) & (df_edges5[FROM_TO+str('_Target')] == SAMPLE_END)]
#df_edges6_backward= df_edges5[(df_edges5[FROM_TO+str('_Source')] == SAMPLE_END) & (df_edges5[FROM_TO+str('_Target')] == SAMPLE_START)]
#products6  = df_edges6['mz_targets'] - df_edges6['mz_sources']
#products6_backward  = df_edges6_backward['mz_targets'] - df_edges6_backward['mz_sources']

##trafo_names = ['H2', 'CO', 'O', 'H3N']           ### automatize later
##trafo_names = df_edges['molecular_formula'].unique().tolist()       ## doesn't work.. maybe because of order below
#trafo_names = ['H2', 'O', 'CO', 'H3N', 'CHN', 'CH3N', 'H2O', 'CO2', 'CH2O', 'CHNO']
#for x in trafo_names:
#    TRAFO = x
#    globals()['df_edges6_%s' % x] = df_edges6[(df_edges6['molecular formula'].str.strip() == TRAFO)]
#    globals()['products6_%s' % x]  = globals()['df_edges6_%s' % x]['mz_targets'] - globals()['df_edges6_%s' % x]['mz_sources']   
#    globals()['df_edges6_backward_%s' % x] = df_edges6_backward[(df_edges6_backward['molecular formula'].str.strip() == TRAFO)]
#    globals()['products6_backward_%s' % x] = globals()['df_edges6_backward_%s' % x]['mz_targets'] - globals()['df_edges6_backward_%s' % x]['mz_sources']

##forward = [products6_H2.shape[0], products6_O.shape[0], products6_CO.shape[0], products6_H3N.shape[0]]
##backward = [products6_backward_H2.shape[0], products6_backward_O.shape[0], products6_backward_CO.shape[0], products6_backward_H3N.shape[0]]
##df = pd.DataFrame({'Forward+Backward': np.add(forward, backward), 'Forward': forward,  'Backward': backward}, index=l)
#forward = [products6_H2.shape[0], products6_O.shape[0], products6_CO.shape[0], products6_H3N.shape[0], products6_CHN.shape[0], products6_CH3N.shape[0], products6_H2O.shape[0], products6_CO2.shape[0], products6_CH2O.shape[0], products6_CHNO.shape[0]]
#backward = [products6_backward_H2.shape[0], products6_backward_O.shape[0], products6_backward_CO.shape[0], products6_backward_H3N.shape[0], products6_backward_CHN.shape[0], products6_backward_CH3N.shape[0], products6_backward_H2O.shape[0], products6_backward_CO2.shape[0], products6_backward_CH2O.shape[0], products6_backward_CHNO.shape[0]]
#df = pd.DataFrame({'Forward': forward,  'Backward': backward}, index=trafo_names)

#fontsize = 12
##plt.figure()
##df.plot.bar(rot=0, stacked=True)
##df.plot(kind='bar', rot=0,  width=.8, color={"Forward+Backward": "lightgray", "Forward": "black", 'Backward': 'w'},  edgecolor='black')
#df.plot(kind='bar', rot=0,  width=.8, color={"Forward": "black", 'Backward': 'w'},  edgecolor='black', stacked=True, fontsize=fontsize)
#plt.title(SAMPLE_START+r'$\rightarrow$'+SAMPLE_END+'    ('+str(name2)+', Component ' +str(FILTER_COMPONENTS)+')', fontsize=fontsize+2)
##plt.xticks()
##plt.yticks()
#plt.xlabel('Transformations', fontsize=fontsize+2)
#plt.ylabel('Frequency', fontsize=fontsize+2)
#plt.legend(fontsize=fontsize)
#plt.savefig(name+'_Component'+str(FILTER_COMPONENTS)+'_'+SAMPLE_START+SAMPLE_END+'    trafo_bar_chart'+'.png')

##plt.show()






###   -----------------------          old long version     -----------------------------------------

##########################
##     CLUSTER TRANSITIONS                              fig 2: min_set     directions between nh3_poor to medium, etc.. 
##         WITH DIRECTIONS
##########################
#df_samples = df_nodes[['id', 'Ice composition H2O_CH3OH_NH3']]
#df_edges4 = df_edges_cc.rename(columns={'Source':'id'})
#df_edges4 = df_edges4.merge(df_samples, how='left', on='id')
#df_edges4 = df_edges4.rename(columns={'Ice composition H2O_CH3OH_NH3': 'H2O_CH3OH_NH3_Source', 'id':'Source'})
#df_edges5 = df_edges4.rename(columns={'Target':'id'})
#df_edges5 = df_edges5.merge(df_samples, how='left', on='id')
#df_edges5 = df_edges5.rename(columns={'Ice composition H2O_CH3OH_NH3': 'H2O_CH3OH_NH3_Target', 'id':'Target'})

##~~~~~~~~~~~~~~~~~~~~~~~
#SAMPLE_START = '3_1_0.2'
#SAMPLE_END   = '3_1_1'
##~~~~~~~~~~~~~~~~~~~~~~~
#df_edges6 = df_edges5[(df_edges5['H2O_CH3OH_NH3_Source'] == SAMPLE_START) & (df_edges5['H2O_CH3OH_NH3_Target'] == SAMPLE_END)]
##~~~~~~~~~~~~~~~~~~~~~~~
#SAMPLE_START = '3_1_1'
#SAMPLE_END   = '3_1_0.2'
##~~~~~~~~~~~~~~~~~~~~~~~
#df_edges6_backward = df_edges5[(df_edges5['H2O_CH3OH_NH3_Source'] == SAMPLE_START) & (df_edges5['H2O_CH3OH_NH3_Target'] == SAMPLE_END)]

##~~~~~~~~~~~~~~~~~~~~~~~
#SAMPLE_START = '3_1_1'
#SAMPLE_END   = '10_1_1'
##~~~~~~~~~~~~~~~~~~~~~~~
#df_edges7 = df_edges5[(df_edges5['H2O_CH3OH_NH3_Source'] == SAMPLE_START) & (df_edges5['H2O_CH3OH_NH3_Target'] == SAMPLE_END)]
##~~~~~~~~~~~~~~~~~~~~~~~
#SAMPLE_START = '10_1_1'
#SAMPLE_END   = '3_1_1'
##~~~~~~~~~~~~~~~~~~~~~~~
#df_edges7_backward = df_edges5[(df_edges5['H2O_CH3OH_NH3_Source'] == SAMPLE_START) & (df_edges5['H2O_CH3OH_NH3_Target'] == SAMPLE_END)]

##~~~~~~~~~~~~~~~~~~~~~~~
#SAMPLE_START = '10_1_1'
#SAMPLE_END   = '3_1_5'
##~~~~~~~~~~~~~~~~~~~~~~~
#df_edges8 = df_edges5[(df_edges5['H2O_CH3OH_NH3_Source'] == SAMPLE_START) & (df_edges5['H2O_CH3OH_NH3_Target'] == SAMPLE_END)]
##~~~~~~~~~~~~~~~~~~~~~~~
#SAMPLE_START = '3_1_5'
#SAMPLE_END   = '10_1_1'
##~~~~~~~~~~~~~~~~~~~~~~~
#df_edges8_backward = df_edges5[(df_edges5['H2O_CH3OH_NH3_Source'] == SAMPLE_START) & (df_edges5['H2O_CH3OH_NH3_Target'] == SAMPLE_END)]

#products6  = df_edges6['mz_targets'] - df_edges6['mz_sources']
#products7  = df_edges7['mz_targets'] - df_edges7['mz_sources']
#products8  = df_edges8['mz_targets'] - df_edges8['mz_sources']

#df_edges6_H2 = df_edges6[(df_edges6['molecular formula'].str.strip() == 'H2')]
#products6_H2  = df_edges6_H2['mz_targets'] - df_edges6_H2['mz_sources']
#df_edges6_backward_H2 = df_edges6_backward[(df_edges6_backward['molecular formula'].str.strip() == 'H2')]
#products6_backward_H2  = df_edges6_backward_H2['mz_targets'] - df_edges6_backward_H2['mz_sources']
#df_edges7_H2 = df_edges7[(df_edges7['molecular formula'].str.strip() == 'H2')]
#products7_H2  = df_edges7_H2['mz_targets'] - df_edges7_H2['mz_sources']
#df_edges7_backward_H2 = df_edges7_backward[(df_edges7_backward['molecular formula'].str.strip() == 'H2')]
#products7_backward_H2  = df_edges7_backward_H2['mz_targets'] - df_edges7_backward_H2['mz_sources']
#df_edges8_H2 = df_edges8[(df_edges8['molecular formula'].str.strip() == 'H2')]
#products8_H2  = df_edges8_H2['mz_targets'] - df_edges8_H2['mz_sources']
#df_edges8_backward_H2 = df_edges8_backward[(df_edges8_backward['molecular formula'].str.strip() == 'H2')]
#products8_backward_H2  = df_edges8_backward_H2['mz_targets'] - df_edges8_backward_H2['mz_sources']

#df_edges6_O = df_edges6[(df_edges6['molecular formula'].str.strip() == 'O')]
#products6_O  = df_edges6_O['mz_targets'] - df_edges6_O['mz_sources']
#df_edges6_backward_O = df_edges6_backward[(df_edges6_backward['molecular formula'].str.strip() == 'O')]
#products6_backward_O  = df_edges6_backward_O['mz_targets'] - df_edges6_backward_O['mz_sources']
#df_edges7_O = df_edges7[(df_edges7['molecular formula'].str.strip() == 'O')]
#products7_O  = df_edges7_O['mz_targets'] - df_edges7_O['mz_sources']
#df_edges7_backward_O = df_edges7_backward[(df_edges7_backward['molecular formula'].str.strip() == 'O')]
#products7_backward_O  = df_edges7_backward_O['mz_targets'] - df_edges7_backward_O['mz_sources']
#df_edges8_O = df_edges8[(df_edges8['molecular formula'].str.strip() == 'O')]
#products8_O  = df_edges8_O['mz_targets'] - df_edges8_O['mz_sources']
#df_edges8_backward_O = df_edges8_backward[(df_edges8_backward['molecular formula'].str.strip() == 'O')]
#products8_backward_O  = df_edges8_backward_O['mz_targets'] - df_edges8_backward_O['mz_sources']

#df_edges6_CO = df_edges6[(df_edges6['molecular formula'].str.strip() == 'CO')]
#products6_CO  = df_edges6_CO['mz_targets'] - df_edges6_CO['mz_sources']
#df_edges6_backward_CO = df_edges6_backward[(df_edges6_backward['molecular formula'].str.strip() == 'CO')]
#products6_backward_CO  = df_edges6_backward_CO['mz_targets'] - df_edges6_backward_CO['mz_sources']
#df_edges7_CO = df_edges7[(df_edges7['molecular formula'].str.strip() == 'CO')]
#products7_CO  = df_edges7_CO['mz_targets'] - df_edges7_CO['mz_sources']
#df_edges7_backward_CO = df_edges7_backward[(df_edges7_backward['molecular formula'].str.strip() == 'CO')]
#products7_backward_CO  = df_edges7_backward_CO['mz_targets'] - df_edges7_backward_CO['mz_sources']
#df_edges8_CO = df_edges8[(df_edges8['molecular formula'].str.strip() == 'CO')]
#products8_CO  = df_edges8_CO['mz_targets'] - df_edges8_CO['mz_sources']
#df_edges8_backward_CO = df_edges8_backward[(df_edges8_backward['molecular formula'].str.strip() == 'CO')]
#products8_backward_CO  = df_edges8_backward_CO['mz_targets'] - df_edges8_backward_CO['mz_sources']

#df_edges6_H3N = df_edges6[(df_edges6['molecular formula'].str.strip() == 'H3N')]
#products6_H3N  = df_edges6_H3N['mz_targets'] - df_edges6_H3N['mz_sources']
#df_edges6_backward_H3N = df_edges6_backward[(df_edges6_backward['molecular formula'].str.strip() == 'H3N')]
#products6_backward_H3N  = df_edges6_backward_H3N['mz_targets'] - df_edges6_backward_H3N['mz_sources']
#df_edges7_H3N = df_edges7[(df_edges7['molecular formula'].str.strip() == 'H3N')]
#products7_H3N  = df_edges7_H3N['mz_targets'] - df_edges7_H3N['mz_sources']
#df_edges7_backward_H3N = df_edges7_backward[(df_edges7_backward['molecular formula'].str.strip() == 'H3N')]
#products7_backward_H3N  = df_edges7_backward_H3N['mz_targets'] - df_edges7_backward_H3N['mz_sources']
#df_edges8_H3N = df_edges8[(df_edges8['molecular formula'].str.strip() == 'H3N')]
#products8_H3N  = df_edges8_H3N['mz_targets'] - df_edges8_H3N['mz_sources']
#df_edges8_backward_H3N = df_edges8_backward[(df_edges8_backward['molecular formula'].str.strip() == 'H3N')]
#products8_backward_H3N  = df_edges8_backward_H3N['mz_targets'] - df_edges8_backward_H3N['mz_sources']

#df_edges6_CHN = df_edges6[(df_edges6['molecular formula'].str.strip() == 'CHN')]
#products6_CHN  = df_edges6_CHN['mz_targets'] - df_edges6_CHN['mz_sources']
#df_edges6_backward_CHN = df_edges6_backward[(df_edges6_backward['molecular formula'].str.strip() == 'CHN')]
#products6_backward_CHN  = df_edges6_backward_CHN['mz_targets'] - df_edges6_backward_CHN['mz_sources']
#df_edges7_CHN = df_edges7[(df_edges7['molecular formula'].str.strip() == 'CHN')]
#products7_CHN  = df_edges7_CHN['mz_targets'] - df_edges7_CHN['mz_sources']
#df_edges7_backward_CHN = df_edges7_backward[(df_edges7_backward['molecular formula'].str.strip() == 'CHN')]
#products7_backward_CHN  = df_edges7_backward_CHN['mz_targets'] - df_edges7_backward_CHN['mz_sources']
#df_edges8_CHN = df_edges8[(df_edges8['molecular formula'].str.strip() == 'CHN')]
#products8_CHN  = df_edges8_CHN['mz_targets'] - df_edges8_CHN['mz_sources']
#df_edges8_backward_CHN = df_edges8_backward[(df_edges8_backward['molecular formula'].str.strip() == 'CHN')]
#products8_backward_CHN  = df_edges8_backward_CHN['mz_targets'] - df_edges8_backward_CHN['mz_sources']

#df_edges6_CH3N = df_edges6[(df_edges6['molecular formula'].str.strip() == 'CH3N')]
#products6_CH3N  = df_edges6_CH3N['mz_targets'] - df_edges6_CH3N['mz_sources']
#df_edges6_backward_CH3N = df_edges6_backward[(df_edges6_backward['molecular formula'].str.strip() == 'CH3N')]
#products6_backward_CH3N  = df_edges6_backward_CH3N['mz_targets'] - df_edges6_backward_CH3N['mz_sources']
#df_edges7_CH3N = df_edges7[(df_edges7['molecular formula'].str.strip() == 'CH3N')]
#products7_CH3N  = df_edges7_CH3N['mz_targets'] - df_edges7_CH3N['mz_sources']
#df_edges7_backward_CH3N = df_edges7_backward[(df_edges7_backward['molecular formula'].str.strip() == 'CH3N')]
#products7_backward_CH3N  = df_edges7_backward_CH3N['mz_targets'] - df_edges7_backward_CH3N['mz_sources']
#df_edges8_CH3N = df_edges8[(df_edges8['molecular formula'].str.strip() == 'CH3N')]
#products8_CH3N  = df_edges8_CH3N['mz_targets'] - df_edges8_CH3N['mz_sources']
#df_edges8_backward_CH3N = df_edges8_backward[(df_edges8_backward['molecular formula'].str.strip() == 'CH3N')]
#products8_backward_CH3N  = df_edges8_backward_CH3N['mz_targets'] - df_edges8_backward_CH3N['mz_sources']

#df_edges6_H2O = df_edges6[(df_edges6['molecular formula'].str.strip() == 'H2O')]
#products6_H2O  = df_edges6_H2O['mz_targets'] - df_edges6_H2O['mz_sources']
#df_edges6_backward_H2O = df_edges6_backward[(df_edges6_backward['molecular formula'].str.strip() == 'H2O')]
#products6_backward_H2O  = df_edges6_backward_H2O['mz_targets'] - df_edges6_backward_H2O['mz_sources']
#df_edges7_H2O = df_edges7[(df_edges7['molecular formula'].str.strip() == 'H2O')]
#products7_H2O  = df_edges7_H2O['mz_targets'] - df_edges7_H2O['mz_sources']
#df_edges7_backward_H2O = df_edges7_backward[(df_edges7_backward['molecular formula'].str.strip() == 'H2O')]
#products7_backward_H2O  = df_edges7_backward_H2O['mz_targets'] - df_edges7_backward_H2O['mz_sources']
#df_edges8_H2O = df_edges8[(df_edges8['molecular formula'].str.strip() == 'H2O')]
#products8_H2O  = df_edges8_H2O['mz_targets'] - df_edges8_H2O['mz_sources']
#df_edges8_backward_H2O = df_edges8_backward[(df_edges8_backward['molecular formula'].str.strip() == 'H2O')]
#products8_backward_H2O  = df_edges8_backward_H2O['mz_targets'] - df_edges8_backward_H2O['mz_sources']

#df_edges6_CO2 = df_edges6[(df_edges6['molecular formula'].str.strip() == 'CO2')]
#products6_CO2  = df_edges6_CO2['mz_targets'] - df_edges6_CO2['mz_sources']
#df_edges6_backward_CO2 = df_edges6_backward[(df_edges6_backward['molecular formula'].str.strip() == 'CO2')]
#products6_backward_CO2  = df_edges6_backward_CO2['mz_targets'] - df_edges6_backward_CO2['mz_sources']
#df_edges7_CO2 = df_edges7[(df_edges7['molecular formula'].str.strip() == 'CO2')]
#products7_CO2  = df_edges7_CO2['mz_targets'] - df_edges7_CO2['mz_sources']
#df_edges7_backward_CO2 = df_edges7_backward[(df_edges7_backward['molecular formula'].str.strip() == 'CO2')]
#products7_backward_CO2  = df_edges7_backward_CO2['mz_targets'] - df_edges7_backward_CO2['mz_sources']
#df_edges8_CO2 = df_edges8[(df_edges8['molecular formula'].str.strip() == 'CO2')]
#products8_CO2  = df_edges8_CO2['mz_targets'] - df_edges8_CO2['mz_sources']
#df_edges8_backward_CO2 = df_edges8_backward[(df_edges8_backward['molecular formula'].str.strip() == 'CO2')]
#products8_backward_CO2  = df_edges8_backward_CO2['mz_targets'] - df_edges8_backward_CO2['mz_sources']

#df_edges6_CH2O = df_edges6[(df_edges6['molecular formula'].str.strip() == 'CH2O')]
#products6_CH2O  = df_edges6_CH2O['mz_targets'] - df_edges6_CH2O['mz_sources']
#df_edges6_backward_CH2O = df_edges6_backward[(df_edges6_backward['molecular formula'].str.strip() == 'CH2O')]
#products6_backward_CH2O  = df_edges6_backward_CH2O['mz_targets'] - df_edges6_backward_CH2O['mz_sources']
#df_edges7_CH2O = df_edges7[(df_edges7['molecular formula'].str.strip() == 'CH2O')]
#products7_CH2O  = df_edges7_CH2O['mz_targets'] - df_edges7_CH2O['mz_sources']
#df_edges7_backward_CH2O = df_edges7_backward[(df_edges7_backward['molecular formula'].str.strip() == 'CH2O')]
#products7_backward_CH2O  = df_edges7_backward_CH2O['mz_targets'] - df_edges7_backward_CH2O['mz_sources']
#df_edges8_CH2O = df_edges8[(df_edges8['molecular formula'].str.strip() == 'CH2O')]
#products8_CH2O  = df_edges8_CH2O['mz_targets'] - df_edges8_CH2O['mz_sources']
#df_edges8_backward_CH2O = df_edges8_backward[(df_edges8_backward['molecular formula'].str.strip() == 'CH2O')]
#products8_backward_CH2O  = df_edges8_backward_CH2O['mz_targets'] - df_edges8_backward_CH2O['mz_sources']

#df_edges6_CHNO = df_edges6[(df_edges6['molecular formula'].str.strip() == 'CHNO')]
#products6_CHNO  = df_edges6_CHNO['mz_targets'] - df_edges6_CHNO['mz_sources']
#df_edges6_backward_CHNO = df_edges6_backward[(df_edges6_backward['molecular formula'].str.strip() == 'CHNO')]
#products6_backward_CHNO  = df_edges6_backward_CHNO['mz_targets'] - df_edges6_backward_CHNO['mz_sources']
#df_edges7_CHNO = df_edges7[(df_edges7['molecular formula'].str.strip() == 'CHNO')]
#products7_CHNO  = df_edges7_CHNO['mz_targets'] - df_edges7_CHNO['mz_sources']
#df_edges7_backward_CHNO = df_edges7_backward[(df_edges7_backward['molecular formula'].str.strip() == 'CHNO')]
#products7_backward_CHNO  = df_edges7_backward_CHNO['mz_targets'] - df_edges7_backward_CHNO['mz_sources']
#df_edges8_CHNO = df_edges8[(df_edges8['molecular formula'].str.strip() == 'CHNO')]
#products8_CHNO  = df_edges8_CHNO['mz_targets'] - df_edges8_CHNO['mz_sources']
#df_edges8_backward_CHNO = df_edges8_backward[(df_edges8_backward['molecular formula'].str.strip() == 'CHNO')]
#products8_backward_CHNO  = df_edges8_backward_CHNO['mz_targets'] - df_edges8_backward_CHNO['mz_sources']
#direction_H2_6    =    products6_H2.shape[0]    /    products6_backward_H2.shape[0]
#direction_O_6    =    products6_O.shape[0]    /    products6_backward_O.shape[0]
#direction_CO_6    =    products6_CO.shape[0]    /    products6_backward_CO.shape[0]
#direction_H3N_6    =    products6_H3N.shape[0]    /    products6_backward_H3N.shape[0]
#direction_CHN_6    =    products6_CHN.shape[0]    /    products6_backward_CHN.shape[0]
#direction_CH3N_6    =    products6_CH3N.shape[0]    /    products6_backward_CH3N.shape[0]
#direction_H2O_6    =    products6_H2O.shape[0]    /    products6_backward_H2O.shape[0]
#direction_CO2_6    =    products6_CO2.shape[0]    /    products6_backward_CO2.shape[0]
#direction_CH2O_6    =    products6_CH2O.shape[0]    /    products6_backward_CH2O.shape[0]
#direction_CHNO_6    =    products6_CHNO.shape[0]    /    products6_backward_CHNO.shape[0]

#direction_H2_7    =    products7_H2.shape[0]    /    products7_backward_H2.shape[0]
#direction_O_7    =    products7_O.shape[0]    /    products7_backward_O.shape[0]
#direction_CO_7    =    products7_CO.shape[0]    /    products7_backward_CO.shape[0]
#direction_H3N_7    =    products7_H3N.shape[0]    /    products7_backward_H3N.shape[0]
#direction_CHN_7    =    products7_CHN.shape[0]    /    products7_backward_CHN.shape[0]
#direction_CH3N_7    =    products7_CH3N.shape[0]    /    products7_backward_CH3N.shape[0]
#direction_H2O_7    =    products7_H2O.shape[0]    /    products7_backward_H2O.shape[0]
#direction_CO2_7    =    products7_CO2.shape[0]    /    products7_backward_CO2.shape[0]
#direction_CH2O_7    =    products7_CH2O.shape[0]    /    products7_backward_CH2O.shape[0]
#direction_CHNO_7    =    products7_CHNO.shape[0]    /    products7_backward_CHNO.shape[0]

#direction_H2_8    =    products8_H2.shape[0]    /    products8_backward_H2.shape[0]
#direction_O_8    =    products8_O.shape[0]    /    products8_backward_O.shape[0]
#direction_CO_8    =    products8_CO.shape[0]    /    products8_backward_CO.shape[0]
#direction_H3N_8    =    products8_H3N.shape[0]    /    products8_backward_H3N.shape[0]
#direction_CHN_8    =    products8_CHN.shape[0]    /    products8_backward_CHN.shape[0]
#direction_CH3N_8    =    products8_CH3N.shape[0]    /    products8_backward_CH3N.shape[0]
#direction_H2O_8    =    products8_H2O.shape[0]    /    products8_backward_H2O.shape[0]
#direction_CO2_8    =    products8_CO2.shape[0]    /    products8_backward_CO2.shape[0]
#direction_CH2O_8    =    products8_CH2O.shape[0]    /    products8_backward_CH2O.shape[0]
#direction_CHNO_8    =    products8_CHNO.shape[0]    /    products8_backward_CHNO.shape[0]

#print('Direction_H2_6:        '+str(direction_H2_6))
#print('Direction_O_6:        '+str(direction_O_6))
#print('Direction_CO_6:        '+str(direction_CO_6))
#print('Direction_H3N_6:        '+str(direction_H3N_6))
#print('Direction_CHN_6:        '+str(direction_CHN_6))
#print('Direction_CH3N_6:        '+str(direction_CH3N_6))
#print('Direction_H2O_6:        '+str(direction_H2O_6))
#print('Direction_CO2_6:        '+str(direction_CO2_6))
#print('Direction_CH2O_6:        '+str(direction_CH2O_6))
#print('Direction_CHNO_6:        '+str(direction_CHNO_6))

#print('Direction_H2_7:        '+str(direction_H2_7))
#print('Direction_O_7:        '+str(direction_O_7))
#print('Direction_CO_7:        '+str(direction_CO_7))
#print('Direction_H3N_7:        '+str(direction_H3N_7))
#print('Direction_CHN_7:        '+str(direction_CHN_7))
#print('Direction_CH3N_7:        '+str(direction_CH3N_7))
#print('Direction_H2O_7:        '+str(direction_H2O_7))
#print('Direction_CO2_7:        '+str(direction_CO2_7))
#print('Direction_CH2O_7:        '+str(direction_CH2O_7))
#print('Direction_CHNO_7:        '+str(direction_CHNO_7))

#print('Direction_H2_8:        '+str(direction_H2_8))
#print('Direction_O_8:        '+str(direction_O_8))
#print('Direction_CO_8:        '+str(direction_CO_8))
#print('Direction_H3N_8:        '+str(direction_H3N_8))
#print('Direction_CHN_8:        '+str(direction_CHN_8))
#print('Direction_CH3N_8:        '+str(direction_CH3N_8))
#print('Direction_H2O_8:        '+str(direction_H2O_8))
#print('Direction_CO2_8:        '+str(direction_CO2_8))
#print('Direction_CH2O_8:        '+str(direction_CH2O_8))
#print('Direction_CHNO_8:        '+str(direction_CHNO_8))


#index = ['H2', 'O', 'CO', 'H3N', 'CHN', 'CH3N', 'H2O', 'CO2', 'CH2O', 'CHNO']

#forward = [products6_H2.shape[0], products6_O.shape[0], products6_CO.shape[0], products6_H3N.shape[0], products6_CHN.shape[0], products6_CH3N.shape[0], products6_H2O.shape[0], products6_CO2.shape[0], products6_CH2O.shape[0], products6_CHNO.shape[0]]
#backward = [products6_backward_H2.shape[0], products6_backward_O.shape[0], products6_backward_CO.shape[0], products6_backward_H3N.shape[0], products6_backward_CHN.shape[0], products6_backward_CH3N.shape[0], products6_backward_H2O.shape[0], products6_backward_CO2.shape[0], products6_backward_CH2O.shape[0], products6_backward_CHNO.shape[0]]
##spacer = a = [0] * len(forward)
#df = pd.DataFrame({'Forward+Backward': np.add(forward, backward), 'Forward': forward,  'Backward': backward}, index=index)

#plt.figure()
#df.plot(kind='bar', rot=0,  width=.8, color={"Forward+Backward": "lightgray", "Forward": "black", 'Backward': 'w'},  edgecolor='black')
##ax = df.plot.bar(rot=0, stacked=True)
##~~~~~~~~~~~~~~~~~~~~~~~
#SAMPLE_START = '3_1_02'
#SAMPLE_END   = '3_1_1'
##~~~~~~~~~~~~~~~~~~~~~~~
#plt.title(SAMPLE_START+r'$\rightarrow$'+SAMPLE_END+'  (H$_2$O:CH$_3$OH:NH$_3$), Component '+str(FILTER_COMPONENTS))
#plt.xticks()
#plt.yticks()
#plt.xlabel('Transformations')
#plt.ylabel('Frequency')
#plt.savefig(name+'_Component'+str(FILTER_COMPONENTS)+'_'+SAMPLE_START+SAMPLE_END+'    trafo_bar_chart'+'.png')

#forward = [products7_H2.shape[0], products7_O.shape[0], products7_CO.shape[0], products7_H3N.shape[0], products7_CHN.shape[0], products7_CH3N.shape[0], products7_H2O.shape[0], products7_CO2.shape[0], products7_CH2O.shape[0], products7_CHNO.shape[0]]
#backward = [products7_backward_H2.shape[0], products7_backward_O.shape[0], products7_backward_CO.shape[0], products7_backward_H3N.shape[0], products7_backward_CHN.shape[0], products7_backward_CH3N.shape[0], products7_backward_H2O.shape[0], products7_backward_CO2.shape[0], products7_backward_CH2O.shape[0], products7_backward_CHNO.shape[0]]
#df = pd.DataFrame({'Forward+Backward': np.add(forward, backward), 'Forward': forward,  'Backward': backward}, index=index)

#plt.figure()
#df.plot(kind='bar', rot=0,  width=.8, color={"Forward+Backward": "lightgray", "Forward": "black", 'Backward': 'w'},  edgecolor='black')
##ax = df.plot.bar(rot=0, stacked=True)
##~~~~~~~~~~~~~~~~~~~~~~~
#SAMPLE_START = '3_1_1'
#SAMPLE_END   = '10_1_1'
##~~~~~~~~~~~~~~~~~~~~~~~
#plt.title(SAMPLE_START+r'$\rightarrow$'+SAMPLE_END+'  (H$_2$O:CH$_3$OH:NH$_3$), Component '+str(FILTER_COMPONENTS))
#plt.xticks()
#plt.yticks()
#plt.xlabel('Transformations')
#plt.ylabel('Frequency')
#plt.savefig(name+'_Component'+str(FILTER_COMPONENTS)+'_'+SAMPLE_START+SAMPLE_END+'    trafo_bar_chart'+'.png')

#forward = [products8_H2.shape[0], products8_O.shape[0], products8_CO.shape[0], products8_H3N.shape[0], products8_CHN.shape[0], products8_CH3N.shape[0], products8_H2O.shape[0], products8_CO2.shape[0], products8_CH2O.shape[0], products8_CHNO.shape[0]]
#backward = [products8_backward_H2.shape[0], products8_backward_O.shape[0], products8_backward_CO.shape[0], products8_backward_H3N.shape[0], products8_backward_CHN.shape[0], products8_backward_CH3N.shape[0], products8_backward_H2O.shape[0], products8_backward_CO2.shape[0], products8_backward_CH2O.shape[0], products8_backward_CHNO.shape[0]]
#df = pd.DataFrame({'Forward+Backward': np.add(forward, backward), 'Forward': forward,  'Backward': backward}, index=index)

#plt.figure()
#df.plot(kind='bar', rot=0,  width=.8, color={"Forward+Backward": "lightgray", "Forward": "black", 'Backward': 'w'},  edgecolor='black')
##ax = df.plot.bar(rot=0, stacked=True)
##~~~~~~~~~~~~~~~~~~~~~~~
#SAMPLE_START = '10_1_1'
#SAMPLE_END   = '3_1_5'
##~~~~~~~~~~~~~~~~~~~~~~~
#plt.title(SAMPLE_START+r'$\rightarrow$'+SAMPLE_END+'  (H$_2$O:CH$_3$OH:NH$_3$), Component '+str(FILTER_COMPONENTS))
#plt.xticks()
#plt.yticks()
#plt.xlabel('Transformations')
#plt.ylabel('Frequency')
#plt.savefig(name+'_Component'+str(FILTER_COMPONENTS)+'_'+SAMPLE_START+SAMPLE_END+'    trafo_bar_chart'+'.png')

##plt.show()





























################################
##          DIRECTIONS OF TRAFOS                           linkers
################################

#########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#name_comp0    = 'EDGES_[H2 ], [O ], [CO ], [H3N ]_min_set  Component [0]       nodes_export'
#name_comp1    = 'EDGES_[H2 ], [O ], [CO ], [H3N ]_min_set  Component [1]       nodes_export'
#########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#df_comp0 = pd.read_csv(name_comp0 + '.txt', sep='\t', header=0)
#df_comp1 = pd.read_csv(name_comp1 + '.txt', sep='\t', header=0)

#df_comp0 = df_comp0                # from comp0 to comp1
#df_comp1 = df_comp1

#df_comp0_id = df_comp0[['id','Label']]
#df_comp0_id = df_comp0_id.rename(columns={'Label':'X'})
#df_edges_cc_id_source = df_edges_cc.rename(columns={'Source':'id'})
#df_source_in_comp0 = df_edges_cc_id_source.merge(df_comp0_id, how='left', on='id')              
#df_source_in_comp0  = df_source_in_comp0[df_source_in_comp0['X'].str.strip() == 'x']
#df_source_in_comp0 = df_source_in_comp0.rename(columns={'id':'Source'})
#df_source_in_comp0 = df_source_in_comp0.drop(columns=['X'])

#df_comp1_id = df_comp1[['id','Label']]
#df_comp1_id = df_comp0_id.rename(columns={'Label':'X'})
#df_edges_cc_id_target = df_source_in_comp0.rename(columns={'Target':'id'})
#df_target_in_comp1 = df_edges_cc_id_target.merge(df_comp1_id, how='left', on='id')              
#df_target_in_comp1  = df_target_in_comp1[df_target_in_comp1['X'].str.strip() == 'x']
#df_target_in_comp1 = df_target_in_comp1.rename(columns={'id':'Target'})
#df_target_in_comp1 = df_target_in_comp1.drop(columns=['X'])

###############  sollte passen.

#df_edges_CH2O2  = df_target_in_comp1[df_target_in_comp1['molecular_formula'].str.strip() == 'CH2O2']                                 # CHN  H3N  CHNO  O  CH3N  CO2  CO  H2O  CH2O  H2
#df_edges_CH4O  = df_target_in_comp1[df_target_in_comp1['molecular_formula'].str.strip() == 'CH4O']
#df_edges_CH5N  = df_target_in_comp1[df_target_in_comp1['molecular_formula'].str.strip() == 'CH5N']

#products_CH2O2  = df_edges_CH2O2['mz_targets'] - df_edges_CH2O2['mz_sources']
#products_CH4O  = df_edges_CH4O['mz_targets'] - df_edges_CH4O['mz_sources']
#products_CH5N  = df_edges_CH5N['mz_targets'] - df_edges_CH5N['mz_sources']

#df_comp0 = df_comp1                # from comp0 to comp1
#df_comp1 = df_comp0

#df_comp0_id = df_comp0[['id','Label']]
#df_comp0_id = df_comp0_id.rename(columns={'Label':'X'})
#df_edges_cc_id_source = df_edges_cc.rename(columns={'Source':'id'})
#df_source_in_comp0 = df_edges_cc_id_source.merge(df_comp0_id, how='left', on='id')              
#df_source_in_comp0  = df_source_in_comp0[df_source_in_comp0['X'].str.strip() == 'x']
#df_source_in_comp0 = df_source_in_comp0.rename(columns={'id':'Source'})
#df_source_in_comp0 = df_source_in_comp0.drop(columns=['X'])

#df_comp1_id = df_comp1[['id','Label']]
#df_comp1_id = df_comp0_id.rename(columns={'Label':'X'})
#df_edges_cc_id_target = df_source_in_comp0.rename(columns={'Target':'id'})
#df_target_in_comp1 = df_edges_cc_id_target.merge(df_comp1_id, how='left', on='id')              
#df_target_in_comp1  = df_target_in_comp1[df_target_in_comp1['X'].str.strip() == 'x']
#df_target_in_comp1 = df_target_in_comp1.rename(columns={'id':'Target'})
#df_target_in_comp1 = df_target_in_comp1.drop(columns=['X'])

#df_edges_CH2O2_backward  = df_target_in_comp1[df_target_in_comp1['molecular_formula'].str.strip() == 'CH2O2']                                 # CHN  H3N  CHNO  O  CH3N  CO2  CO  H2O  CH2O  H2
#df_edges_CH4O_backward  = df_target_in_comp1[df_target_in_comp1['molecular_formula'].str.strip() == 'CH4O']
#df_edges_CH5N_backward  = df_target_in_comp1[df_target_in_comp1['molecular_formula'].str.strip() == 'CH5N']

##products_CH2O2_backward  = df_edges_CH2O2_backward['mz_targets'] - df_edges_CH2O2_backward['mz_sources']
##products_CH4O_backward  = df_edges_CH4O_backward['mz_targets'] - df_edges_CH4O_backward['mz_sources']
##products_CH5N_backward  = df_edges_CH5N_backward['mz_targets'] - df_edges_CH5N_backward['mz_sources']

#direction_CH2O2    =    df_edges_CH2O2.shape[0]    /    df_edges_CH2O2_backward.shape[0]
#direction_CH4O    =    df_edges_CH4O.shape[0]    /    df_edges_CH4O_backward.shape[0]
#direction_CH5N    =    df_edges_CH5N.shape[0]    /    df_edges_CH5N_backward.shape[0]

#print(df_edges_CH2O2.shape[0])
#print(df_edges_CH2O2_backward.shape[0])
#print(df_edges_CH4O.shape[0])
#print(df_edges_CH4O_backward.shape[0])
#print(df_edges_CH5N.shape[0])
#print(df_edges_CH5N_backward.shape[0])
#print(direction_CH2O2)
#print(direction_CH4O)
#print(direction_CH5N)

########index = ['CH2O2', 'CH4O', 'CH5N']
########forward = [df_edges_CH2O2.shape[0], df_edges_CH4O.shape[0], df_edges_CH5N.shape[0]]
########backward = [df_edges_CH2O2_backward.shape[0], df_edges_CH4O_backward.shape[0], df_edges_CH5N_backward.shape[0]]
########df = pd.DataFrame({'Forward+Backward': np.add(forward, backward), 'Forward': forward,  'Backward': backward}, index=index)

########plt.figure()
########df.plot(kind='bar', rot=0,  width=.8, color={"Forward+Backward": "lightgray", "Forward": "black", 'Backward': 'w'},  edgecolor='black')
#########ax = df.plot.bar(rot=0, stacked=True)
#########plt.title(SAMPLE_START+r'$\rightarrow$'+SAMPLE_END+'  (H$_2$O:CH$_3$OH:NH$_3$), Component '+str(FILTER_COMPONENTS))
########plt.xticks()
########plt.yticks()
########plt.xlabel('Transformations')
########plt.ylabel('Frequency')
#########plt.savefig(name+'_Component'+str(FILTER_COMPONENTS)+'_'+SAMPLE_START+SAMPLE_END+'    trafo_bar_chart'+'.png')
#########plt.savefig(name+'_linkers_trafo_bar_chart'+'.png')

#######~~~~~~~~~~~~~~~~~~~~~~~
#TRAFO = 'CH$_2$O$_2$'
#######~~~~~~~~~~~~~~~~~~~~~~~
#forward = [df_edges_CH2O2.shape[0]]
#backward = [df_edges_CH2O2_backward.shape[0]]
#df = pd.DataFrame({'Backward': backward, 'Forward': forward, 'Forward+Backward': np.add(forward, backward)})
#plt.figure()
#df.plot(kind='barh',  width=.8, color={"Forward+Backward": "lightgray", "Forward": "black", 'Backward': 'w'},  edgecolor='black',figsize=(12, 2))
#plt.title(r'Transformation '+TRAFO, fontsize=15)
#plt.xticks()
#plt.yticks()
#plt.xlabel('Frequency')
#plt.savefig(name+'_'+TRAFO+'_linkers_trafo_bar_chart'+'.svg')

#######~~~~~~~~~~~~~~~~~~~~~~~
#TRAFO = 'CH$_4$O'
#######~~~~~~~~~~~~~~~~~~~~~~~
#forward = [df_edges_CH4O.shape[0]]
#backward = [df_edges_CH4O_backward.shape[0]]
#df = pd.DataFrame({'Backward': backward, 'Forward': forward, 'Forward+Backward': np.add(forward, backward)})
#plt.figure()
#df.plot(kind='barh',  width=.8, color={"Forward+Backward": "lightgray", "Forward": "black", 'Backward': 'w'},  edgecolor='black',figsize=(12, 2))
#plt.title(r'Transformation '+TRAFO, fontsize=15)
#plt.xticks()
#plt.yticks()
#plt.xlabel('Frequency')
#plt.savefig(name+'_'+TRAFO+'_linkers_trafo_bar_chart'+'.svg')

#######~~~~~~~~~~~~~~~~~~~~~~~
#TRAFO = 'CH$_5$N'
#######~~~~~~~~~~~~~~~~~~~~~~~
#forward = [df_edges_CH5N.shape[0]]
#backward = [df_edges_CH5N_backward.shape[0]]
#df = pd.DataFrame({'Backward': backward, 'Forward': forward, 'Forward+Backward': np.add(forward, backward)})
#plt.figure()
#df.plot(kind='barh',  width=.8, color={"Forward+Backward": "lightgray", "Forward": "black", 'Backward': 'w'},  edgecolor='black',figsize=(12, 2))
#plt.title(r'Transformation '+TRAFO, fontsize=15)
#plt.xticks()
#plt.yticks()
#plt.xlabel('Frequency')
#plt.savefig(name+'_'+TRAFO+'_linkers_trafo_bar_chart'+'.svg')
































############################
##    BETWEENNESS CENTRALITY
############################
#bc = nx.betweenness_centrality(G_cc, normalized=True)
#df_bc = pd.DataFrame.from_dict(bc, orient='index', columns=['Betweeness centrality']).rename_axis('id').reset_index()
#df_nodes_bc = df_nodes.merge(df_bc, how='left', on='id')

#node_sizes = df_bc['Betweeness centrality'].tolist()
#node_sizes = [10**(x) for x in node_sizes]
##node_sizes = [10**( x / max(node_sizes) )*5000 for x in node_sizes]
#df_node_sizes = pd.DataFrame(node_sizes)

#plt.figure()
#df_node_sizes.plot()
#plt.savefig(name+'    '+color_type+'    betCen_freq(normalization)'+'.png')

#plt.figure()
#nx.draw(G_cc, pos = nx.nx_pydot.graphviz_layout(G_cc), with_labels=False, node_color=node_colors, width=.03, font_size=5, node_size=node_sizes) 
#plt.savefig(name+'    '+color_type+'    network_bc'+'.png')











#############################
##                EXPORT DATA
#############################

#df_export = df_nodes_cc
#df_export.to_csv(name+'  Component '+str(FILTER_COMPONENTS)+'       nodes_export.txt', header="'\t'.join(list(df_export.columns))", index=None, sep='\t')

#df_export = df_edges_cc
#df_export.to_csv(name+'  Component '+str(FILTER_COMPONENTS)+'       edges_export.txt', header="'\t'.join(list(df_export.columns))", index=None, sep='\t')
































  
  
  
#print('Time code = '+str(datetime.datetime.now() - begin_time)+' [h:min:sec]')
  
  





