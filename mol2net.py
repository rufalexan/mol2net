##﻿python3                                                   





'''
=========================================
          mol2net

Network analysis from a list of molecules
=========================================

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

import datetime
begin_time = datetime.datetime.now()
import os
import glob
import string
import chemparse
import re
import numpy as np
import pandas as pd
from collections import namedtuple
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import pydot
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.ticker as ticker
from mycolorpy import colorlist as mcp
import math


m_h = 1.007825
m_c = 12.
m_n = 14.003074
m_o = 15.994915
m_s = 31.972071
masses = (m_h, m_c, m_n, m_o, m_s)
valence_h = 1
valence_c = 4
valence_n = 3
valence_o = 2
valence_s = 2

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
trafos          =       ['H2','O','CO','NH3']
##trafoWeights    =       [ 1,   5,   1,   2]
trafoWeights    =       [1]*len(trafos)
columnNames     =       ['molecular_formula','abundance_int']             ##  for  molFormula  &  abundance/intensity
nComponents     =       1
sample          =       'Ice composition H2O_CH3OH_NH3'
only_uniques    =       'yes'
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

######################################
##                              TRAFOS
######################################
cf = re.compile('([A-Z][a-z]?)(\d?\d?\d?)')
all_elements = []
for xxxx in trafos:
    for el, count in re.findall(cf, xxxx):
        all_elements += [el]

atoms = set(all_elements)
counts = [[""] + trafos]
for e in atoms:
    count = [e]
    for xxxx in trafos:    
        d = dict(re.findall(cf, xxxx))
        n = d.get(e, 0)
        if n == '': n = 1
        count += [int(n)]
    counts += [count]

A = np.array(counts)
A_t = (A.transpose())

x2 = pd.DataFrame(A_t)
new_header = x2.iloc[0]
x2 = x2[1:]
x2.columns = new_header
x2 = x2.rename(columns={'': 'Molecular Formula'})
for i in ['H', 'C', 'N', 'O', 'S']:
    if i in x2.columns:
        x2[i] = x2[i].astype(int)
    else:
        x2[i] = 0    

h = x2['H'].to_numpy()
c = x2['C'].to_numpy()
n = x2['N'].to_numpy()    
o = x2['O'].to_numpy()
s = x2['S'].to_numpy()
elements = np.transpose(np.array([h, c, n, o, s]))
massExact = np.dot(elements, masses)
x2['Mass (exact)'] = massExact.round(6)
weights = pd.DataFrame(trafoWeights,columns = ['trafoWeights'])
weights['trafos'] = trafos
x3 = x2.join(weights.set_index(['trafos'], verify_integrity=True ),on=['Molecular Formula'], how='left')

xTRAFOS = x3
mdTRAFOS = xTRAFOS['Mass (exact)'].round(6).to_numpy()

######################################
##                          DETECTIONS
######################################
path = os.getcwd()
files = glob.glob(os.path.join("*.tsv"))
d = {}
for f in files:
    key = f
    name = f.replace(".tsv", "")
    x = globals()[f"orig_{name}"] = pd.read_csv(f, sep='\t', skiprows=0)
    if 'id' in x.columns:
        pass
    else:
        x.insert(0, 'id', range(1, 1 + len(x)))
#    
    try:
        x = x[['id',sample,columnNames[0],columnNames[1]]].rename(columns={columnNames[0]: 'Molecular Formula',columnNames[1]:'Intensity'})
    except:
        x = x[['id',columnNames[0],columnNames[1]]].rename(columns={columnNames[0]: 'Molecular Formula',columnNames[1]:'Intensity'})
#
    molFormulas = x['Molecular Formula']
    a = []
    for i in molFormulas:
        tmp_a = chemparse.parse_formula(i)
        a.append(tmp_a)
#
    x_a = pd.DataFrame(a).fillna(0)
    x2 = pd.concat([x,x_a], axis=1)
#
    if 'H' in x2.columns:
        x2['H'] = x2['H'].astype(int)
    else:
        x2['H'] = 0    
#
    if 'C' in x2.columns:
        x2['C'] = x2['C'].astype(int)
    else:
        x2['C'] = 0    
#
    if 'N' in x2.columns:
        x2['N'] = x2['N'].astype(int)
    else:
        x2['N'] = 0    
#
    if 'O' in x2.columns:
        x2['O'] = x2['O'].astype(int)
    else:
        x2['O'] = 0    
#
    if 'S' in x2.columns:
        x2['S'] = x2['S'].astype(int)
    else:
        x2['S'] = 0    
#
    h = x2['H'].to_numpy()
    c = x2['C'].to_numpy()
    o = x2['O'].to_numpy()
    n = x2['N'].to_numpy()    
    s = x2['S'].to_numpy()
    elements = np.transpose(np.array([h, c, n, o, s]))
    massExact = np.dot(elements, masses)
    x2['Mass (exact)'] = massExact.round(6)
    x2['mu'] = mu = ( (h*valence_h + c*valence_c + n*valence_n + o*valence_o + s*valence_s) / 2 ) - (h + c + o + n + s) + 1
    x2['DBE'] = dbe = c - (h/2) + (n/2) + 1
    x2['HC'] = h/c
    x2['NC'] = n/c
    x2['OC'] = o/c
    x2['SC'] = s/c
    x2['FILTER'] =  1 * ( (mu > 0) )
    atom_list = ['H','C','N','O','S']
    for i in range(len(atom_list)-1):
        x2['Molecular Formula'] = x2['Molecular Formula'].replace({atom_list[i]+'0'}, {atom_list[i]}, regex=True)
        x2['Molecular Formula'] = x2['Molecular Formula'].replace({atom_list[i]+'1'+atom_list[i+1]}, {atom_list[i]}, regex=True)
        x2['Molecular Formula'] = x2['Molecular Formula'].replace({'S'+'0'}, {''}, regex=True)
        x2['Molecular Formula'] = x2['Molecular Formula'].replace({'S'+'1'}, {'S'}, regex=True)
#
    x2['Mass (frequency)'] = x2['Mass (exact)'].map(dict(x2['Mass (exact)'].value_counts()))
    try:
        if only_uniques == 'yes':
            x3 = x2.drop_duplicates(subset=['Mass (exact)'], keep='first')
        else:
            pass
    except:
        x3 = x2
#
    xNODES0 = xDET = x3[(x3['FILTER']==1)]
    globals()[f"{name}"] = xDET
    mass2 = mass = xDET['Mass (exact)'].round(6).to_numpy()
    data = xDET[['id','Mass (exact)']]
#
    ##########################################     x x    
    ##                              MD MATCHES
    ##########################################
    md_matches = namedtuple('md_matches', 'md_matches hits')
    new_list = []
    for md in mdTRAFOS:
	    for element in mass:
		    if element+md in mass2:
			    tmp_new_list = md_matches(element, element+md)
			    new_list.append(tmp_new_list)
    #
    matches = np.array(new_list).reshape(len(new_list),2)
    matches = matches[np.argsort(matches[:,0])]
    sources = data.rename(columns={'id': 'Source', 'Mass (exact)': 'Mass (source)'})
    targets = data.rename(columns={'id': 'Target', 'Mass (exact)': 'Mass (target)'})
    matches2 = pd.DataFrame({'Mass (source)': matches[:, 0], 'Mass (target)': matches[:, 1]})
    source_match = matches2.merge(sources, how='left', on=['Mass (source)'])
    target_match = source_match.merge(targets, how='left', on=['Mass (target)'])
    target_match['Mass difference'] = mass_diff = np.around(target_match['Mass (target)'] - target_match['Mass (source)'],6)
    xxxx = xTRAFOS[['Mass (exact)','Molecular Formula']].rename(columns={'Mass (exact)': 'Mass difference'})
    x = target_match.merge(xxxx, how='left', on=['Mass difference'])
    x['type'] = type = np.full((len(mass_diff),1), 'Undirected')
    x['Label'] = Label = np.full((len(mass_diff),1), 'x')
    xxxx = xTRAFOS[['Mass (exact)','trafoWeights']]
    x2 = x.join(xxxx.set_index(['Mass (exact)'], verify_integrity=True ),on=['Mass difference'], how='left').rename(columns={'trafoWeights': 'Weight'}).sort_values(by=['Source'])
    x2['Mass difference (frequency)'] = x2['Mass difference'].map(dict(x2['Mass difference'].value_counts()))
    #
    xEDGES0 = x2
    #
    md_freq = x2[['Molecular Formula','Mass difference (frequency)']]
    md_freq = md_freq.drop_duplicates(subset=['Molecular Formula'], keep='first')
    md_freq = md_freq.sort_values(by=['Mass difference (frequency)'], ascending=False)
    print('\n\n\n\n')
    print(name+', '+str(len(xNODES0))+' nodes,', str(len(xEDGES0))+' edges')
    print(md_freq)
    print('\n\n\n\n')
#
    ##########################################     x x    
    ##                             NETWORK ANA
    ##########################################
    G0 = nx.from_pandas_edgelist(xEDGES0, 'Source', 'Target', create_using=nx.Graph())
    G0cc = list(nx.connected_components(G0))
    d = {name:k for k,comp in enumerate(G0cc) for name in comp}               #dict(enumerate(G0cc))
    df_comp = pd.DataFrame.from_dict(d, orient='index', columns=['component']).rename_axis('id').reset_index()
    nodes2 = xNODES0.merge(df_comp, how='left', on='id')
    nans = nodes2[nodes2.isna().any(axis=1)]
    xNODES0 = nodes2[~nodes2['component'].isnull()]
    #
    df_source = xEDGES0.rename(columns={'Source':'id'})
    edges2 = df_source.merge(df_comp, how='left', on='id')
    edges2 = edges2.rename(columns={'component':'component source','id':'Source'})
    df_target = edges2.rename(columns={'Target':'id'})
    edges3 = df_target.merge(df_comp, how='left', on='id')
    edges3 = edges3.rename(columns={'component':'component target','id':'Target','mass_difference':'mass difference','molecular_formula':'molecular formula'})
    #
    df_comp_freq = pd.DataFrame(xNODES0['component'].value_counts()).rename(columns={'component': 'freq_component'})
    xNODES0['df_comp_freq'] = xNODES0['component'].map(dict(xNODES0['component'].value_counts()))
    df_comp_freq_norm = pd.DataFrame(xNODES0['component'].value_counts(normalize=True)).rename(columns={'component': 'freq_component_norm [0-1]'})
    df_comp_freq_comb = pd.concat([df_comp_freq,df_comp_freq_norm], axis=1)
    print('\n\n\n\n')
    print(name+' (all components)')
    print(df_comp_freq_comb)
    print('\n\n\n\n')
    #
    FILTER_COMPONENTS  =        list(df_comp_freq_comb.index)[:nComponents]
    xEDGES = edges3[(edges3['component source'].isin(FILTER_COMPONENTS))]
    xNODES = xNODES0[(xNODES0['component'].isin(FILTER_COMPONENTS))]
    G = nx.from_pandas_edgelist(xEDGES, 'Source', 'Target', create_using=nx.Graph())
    a = pd.DataFrame(G.nodes(data = True)).rename(columns={0:'id'}).drop([1], axis=1)
    a2 = a.join(xNODES.set_index('id'), on='id')
    a2['id_G0_nodes'] = range(1,1+len(a2))
    xNODES = a2
#
    degrees = pd.DataFrame(G.degree(), columns=["id", "Degree"])
    degrees['id'] = degrees['id'].astype(int)
    xNODES = xNODES.join(degrees.set_index('id'), on='id') 
    xNODES['Degree counts'] = xNODES['Degree'].map(dict(xNODES['Degree'].value_counts()))
#
##    bc = nx.betweenness_centrality(G, normalized=True)                                                                                                                                                                                             ## takes a while....
##    bc2 = pd.DataFrame.from_dict(bc, orient='index', columns=['betCen']).rename_axis('id').reset_index()
##    xNODES = xNODES.merge(bc2, how='left', on='id')
#
    try:
        color_type = sample                                                                                                                                                                                                                                   ##   DAS NOCH VERALLG0EMEINERN ??
    except:
        color_type = 'Degree'
#
    colors = mcp.gen_color(cmap='gist_rainbow',n=len(xNODES[color_type].unique()))                    ##  cmap  AUTUMN tab20b coolwarm gist_rainbow tab20c brg      https://matplotlib.org/stable/tutorials/colors/colormaps.html
    ##colors = ['#15B01A',    'blue', 'orange',     'lightgray',   'red', 'green']
    samples0 = xNODES[color_type].unique().tolist()
    try:
        samples = sorted(samples0, key=lambda x: int("".join([i for i in x if i.isdigit()])))
    except:
        samples = sorted(samples0)
    #
    sample_colors = pd.DataFrame({color_type: samples, 'color': colors})
    xNODES = xNODES.merge(sample_colors, how='left', on=color_type)
    xNODES = xNODES.sort_values(by=['id_G0_nodes'], ascending=True)                              ## re-ordering/re-indexing according to graph G0 (id_G0_nodes)
    node_colors = xNODES['color'].tolist()
#
    node_size = 'Degree'
    xNODES = xNODES.sort_values(by=['id_G0_nodes'], ascending=True)
    xNODES['node_size'] = (xNODES[node_size] - np.min(xNODES[node_size])) / (np.max(xNODES[node_size]) - np.min(xNODES[node_size]))
    node_sizes = xNODES['node_size'].tolist()
    node_sizes = [i * 1e2 for i in node_sizes]
#
###########    l = ['Degree','Time point']
###########    for i0 in l:
###########        x = xNODES[['id',i0]]
###########        for i in ['Source','Target']:
###########            x2 = xEDGES.join(x.set_index(['id'], verify_integrity=True ),on=[i], how='left').rename(columns={i0:i0+' ('+str(i)+')'})
###########            x3 = x2[i0+' ('+str(i)+')']
###########            xEDGES[i0+' ('+str(i)+')'] = x3
#
#
#
##    #########################################
##    ##                    DEGREE DISTRIBUTION
##    #########################################
##    x_name = 'degree'
##    y_name = 'degree counts'
##    text   = 'id'
##    x = xNODES[x_name].tolist()
##    y = xNODES[y_name].tolist()
##    text = xNODES[text].tolist()
###
##    plt.figure()
##    plt.title(name+'\n, Comp'+str(FILTER_COMPONENTS)+', Trafos'+str(trafos), wrap=True, fontsize=12)
##    plt.scatter(x,y,c='blue')
####    for i in range(len(x)):
####        plt.annotate(text[i], (x[i], y[i] + 0.2), fontsize=10)
###
####    plt.xscale('log',base=10) 
####    plt.yscale('log',base=10) 
##    plt.xlabel(x_name)
##    plt.ylabel(y_name)
####    plt.savefig(name+'_'+str(FILTER_COMPONENTS)+'_'+str(trafos)+'    degreeDistri'+'.png')
##    plt.show()
#
    #########################################
    ##                                 LAYOUT
    #########################################
    plt.figure()
    plt.title(name, wrap=True, fontsize=12)
##    nx.draw(G, pos = nx.nx_pydot.graphviz_layout(G), with_labels=False, node_color=node_colors, node_size=20, width=.03, font_size=5)
##    nx.draw(G, pos = nx.nx_pydot.graphviz_layout(G), with_labels=False, node_color=node_colors, node_size=node_sizes, width=.05, alpha=1)
##    nx.draw(G, pos = nx.spring_layout(G, scale=1), with_labels=False, node_color=node_colors, node_size=node_sizes, width=.05, alpha=1)
##    nx.draw(G, pos = nx.spring_layout(G, k=.5/math.sqrt(G.order()), scale=10), with_labels=False, node_color=node_colors, node_size=node_sizes, width=.05, alpha=1)
    nx.draw(G, pos = nx.nx_pydot.graphviz_layout(G), with_labels=False, node_color=node_colors, node_size=2, width=.01, alpha=1)
    for i in range(len(sample_colors)):
        plt.plot([], [], sample_colors.values[i][1], marker='o', markersize=10, label=sample_colors.values[i][0])      
#
    plt.legend()                      ## loc='upper left'
##    plt.text(.6, -1.1,str(G.number_of_nodes())+' nodes, '+str(G.number_of_edges())+' edges, node size ~ '+str(node_size), fontsize=8, wrap=True)
    plt.text(.6, -1.1,str(G.number_of_nodes())+' nodes, '+str(G.number_of_edges())+' edges', fontsize=8, wrap=True)
    plt.savefig(name+'_Comp'+str(FILTER_COMPONENTS)+'    '+str(color_type)+'    network'+'.png')
##    plt.show()
#    
##    ##################################################
##    ##    CONTINOUS COLOR GRADIENT (e.g. element maps)
##    ##################################################
##    fontsize=20
##    l = ['H', 'C', 'N', 'O', 'S', 'HC', 'NC', 'OC', 'SC', 'Mass (exact)', 'Intensity', 'DBE']
###
##    for i in l:
##        color_type = i
##        carac = xNODES[['id', color_type]]
##        carac = carac.set_index('id')
##        carac = carac.reindex(G.nodes())
##        carac[color_type]=pd.Categorical(carac[color_type])
##        carac[color_type].cat.codes
##        nodes = G.nodes()
##        plt.figure()
##        plt.title(color_type, wrap=True, fontsize=25)
####        nx.draw(G, pos = nx.nx_pydot.graphviz_layout(G), with_labels=False, node_color=carac[color_type].cat.codes, cmap=plt.cm.jet, node_size=3, width=.1, font_size=10)
####        pos = nx.spring_layout(G)
##        pos = nx.nx_pydot.graphviz_layout(G)
####        pos = nx.fruchterman_reingold_layout(G)
####        pos = nx.nx_pydot.graphviz_layout(G)
##        nx.draw(G, pos, with_labels=False, node_color=node_colors, node_size=2, width=.01, alpha=1)
##        ec = nx.draw_networkx_edges(G, pos, width=.01, alpha=0.2)
##    #    nc = nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=carac[color_type], with_labels=False, node_size=5, cmap=plt.cm.jet) 
##        nc = nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=carac[color_type], node_size=3, cmap=plt.cm.jet) 
##    ##    plt.colorbar(nc)
##        cb = plt.colorbar(nc, orientation='vertical')
##        nc.figure.axes[0].tick_params(axis="both", labelsize=21)           ## change the label size
##        nc.figure.axes[1].tick_params(axis="y", labelsize=21)              ## change the tick label size of colorbar
##        plt.axis('off')
##        plt.savefig(name+', Comp'+str(FILTER_COMPONENTS)+'    '+color_type+'    network'+'.png')
##        #plt.show()





















##    #########################################   
##    ##                   edges vs edge counts
##    #########################################
##    edge_edgeCounts = xEDGES.drop_duplicates(subset=['mass difference'], keep='first')
##    edge_edgeCounts = edge_edgeCounts.sort_values(by=['mass difference'], ascending = True)
##    x_name = 'mass difference'
##    y_name = 'mass difference counts'
##    text   = 'trafo'
##    x = edge_edgeCounts[x_name].tolist()
##    y = edge_edgeCounts[y_name].tolist()
##    text = edge_edgeCounts[text].tolist()
###
##    plt.figure()
##    plt.title(name, wrap=True, fontsize=12)
##    plt.scatter(x,y,c='black')
##    for i in range(len(x)):
##        plt.annotate(text[i], (x[i], y[i] + 0.2), fontsize=10)
###
##    plt.xlabel('mass (transformation)')
##    plt.ylabel('transformation counts')
##    plt.savefig(name+'_Comp'+str(FILTER_COMPONENTS)+'    '+str(color_type)+'    edge_edgeCounts'+'.png')
##    ##plt.show()
##    
##    #########################################   
##    ##                    degree vs abundance
##    #########################################
##    xNODES2 = xNODES[ (xNODES['sample'] == 'detections') ]
##    x_name = 'degree'
##    y_name = 'abundance'
##    text   = 'molecular_formula'
##    x = xNODES2[x_name].tolist()
##    y = xNODES2[y_name].tolist()
##    ann = xNODES2[text].tolist()

##    fig, ax = plt.subplots()
##    plt.title(name, wrap=True, fontsize=12)
##    ax.scatter(x,y,c='black')
##    for i, txt in enumerate(ann):
##        ax.annotate(txt, (x[i], y[i]))
###
##    plt.yscale('log',base=2) 
##    plt.xlabel(x_name)
##    ##plt.xlabel(y_name)
##    plt.ylabel('log('+str(abundance)+')')
##    plt.savefig(name+'_Comp'+str(FILTER_COMPONENTS)+'    '+str(color_type)+'    degree_abundance'+'.png')
##    ##plt.show()


























#################################
####          COMPONENT HISTOGRAM
#################################
##plt.figure()
##pd.DataFrame(df_nodes_cc['component'].value_counts(normalize=True)).sort_index().plot(kind='bar', rot=0, ylabel='', legend=False, color='black', width=1)
###pd.DataFrame(df_nodes_cc['component'].value_counts(normalize=False)).sort_index().plot(kind='bar', rot=0, ylabel='', legend=False, color='black', width=1)
##plt.title(name+', Component '+str(FILTER_COMPONENTS), wrap=True)          
##plt.xticks(np.arange(0, len(df_comp_freq), 20), fontsize=15)                      ## ice
####plt.xticks(np.arange(0, len(df_comp_freq), 100), fontsize=15)                      ## paris
##plt.yticks(fontsize=15)
##plt.yscale('log')
##plt.xlabel('Component number', fontsize=15)
###plt.ylabel('Frequency (log scale)')
##plt.ylabel('Frequency (log scale, normalized)', fontsize=15)
###plt.ylabel('Frequency (normalized)')
##plt.savefig(name+', Comp'+str(FILTER_COMPONENTS)+ '    hist_components'+'.svg')







#######################################################
#######                                   component map
#######################################################
##fontsize=20
##color_type = 'component'
##carac = df_nodes_cc[['id', color_type]]
##carac = carac.set_index('id')
##carac = carac.reindex(G.nodes())
##carac[color_type]=pd.Categorical(carac[color_type])
##carac[color_type].cat.codes
##nodes = G.nodes()
##plt.figure()
##plt.title(color_type, wrap=True, fontsize=20)
##nx.draw(G, pos = nx.nx_pydot.graphviz_layout(G), with_labels=False, node_color=carac[color_type].cat.codes, cmap=plt.cm.jet, node_size=3, width=.1, font_size=10)
###pos = nx.spring_layout(G)
###pos = nx.fruchterman_reingold_layout(G)
##pos = nx.nx_pydot.graphviz_layout(G)
##ec = nx.draw_networkx_edges(G, pos, alpha=0.2)
###    nc = nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=carac[color_type], with_labels=False, node_size=5, cmap=plt.cm.jet) 
##nc = nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=carac[color_type], node_size=3, cmap=plt.cm.jet) 
####plt.colorbar(nc)

##cb = plt.colorbar(nc, orientation='vertical').set_label(label=color_type, size=fontsize+2)
##nc.figure.axes[0].tick_params(axis="both", labelsize=21)           ## change the label size
##nc.figure.axes[1].tick_params(axis="y", labelsize=21)              ## change the tick label size of colorbar

##plt.axis('off')
##plt.savefig(name+', Comp'+str(FILTER_COMPONENTS)+'    '+color_type+'    network'+'.png')
####plt.show()





###########################
###     CLUSTER TRANSITIONS                              fig 2: min_set     directions between nh3_poor to medium, etc.. 
###         WITH DIRECTIONS
###########################
###~~~~~~~~~~~~~~~~~~~~~~~
##FROM_TO = 'Ice composition H2O_CH3OH_NH3'
###~~~~~~~~~~~~~~~~~~~~~~~
##df_samples = xNODES[['id', FROM_TO]]
##df_edges4 = xEDGES.rename(columns={'Source':'id'})
##df_edges4 = df_edges4.merge(df_samples, how='left', on='id')
##df_edges4 = df_edges4.rename(columns={FROM_TO: FROM_TO+str('_Source'), 'id':'Source'})
##df_edges5 = df_edges4.rename(columns={'Target':'id'})
##df_edges5 = df_edges5.merge(df_samples, how='left', on='id')
##df_edges5 = df_edges5.rename(columns={FROM_TO: FROM_TO+str('_Target'), 'id':'Target'})

###~~~~~~~~~~~~~~~~~~~~~~~
##SAMPLE_START = '3_1_0.2'
##SAMPLE_END   = '3_1_1'

###SAMPLE_START = '3_1_1'
###SAMPLE_END   = '10_1_1'

###SAMPLE_START = '10_1_1'
###SAMPLE_END   = '3_1_5'
###~~~~~~~~~~~~~~~~~~~~~~~
##df_edges6 = df_edges5[(df_edges5[FROM_TO+str('_Source')] == SAMPLE_START) & (df_edges5[FROM_TO+str('_Target')] == SAMPLE_END)]
##df_edges6_backward= df_edges5[(df_edges5[FROM_TO+str('_Source')] == SAMPLE_END) & (df_edges5[FROM_TO+str('_Target')] == SAMPLE_START)]
##products6  = df_edges6['Mass (target)'] - df_edges6['Mass (source)']
##products6_backward  = df_edges6_backward['Mass (target)'] - df_edges6_backward['Mass (source)']

###trafo_names = ['H2', 'CO', 'O', 'H3N']           ### automatize later
###trafo_names = df_edges['molecular_formula'].unique().tolist()       ## doesn't work.. maybe because of order below
###trafo_names = ['H2', 'O', 'CO', 'H3N', 'CHN', 'CH3N', 'H2O', 'CO2', 'CH2O', 'CHNO']
###for x in trafo_names:
##for x in trafos:
##    TRAFO = x
##    globals()['df_edges6_%s' % x] = df_edges6[(df_edges6['Molecular Formula'].str.strip() == TRAFO)]
##    globals()['products6_%s' % x]  = globals()['df_edges6_%s' % x]['Mass (target)'] - globals()['df_edges6_%s' % x]['Mass (source)']   
##    globals()['df_edges6_backward_%s' % x] = df_edges6_backward[(df_edges6_backward['Molecular Formula'].str.strip() == TRAFO)]
##    globals()['products6_backward_%s' % x] = globals()['df_edges6_backward_%s' % x]['Mass (target)'] - globals()['df_edges6_backward_%s' % x]['Mass (source)']

###forward = [products6_H2.shape[0], products6_O.shape[0], products6_CO.shape[0], products6_H3N.shape[0]]
###backward = [products6_backward_H2.shape[0], products6_backward_O.shape[0], products6_backward_CO.shape[0], products6_backward_H3N.shape[0]]
###df = pd.DataFrame({'Forward+Backward': np.add(forward, backward), 'Forward': forward,  'Backward': backward}, index=l)
##forward = [products6_H2.shape[0], products6_O.shape[0], products6_CO.shape[0], products6_H3N.shape[0], products6_CHN.shape[0], products6_CH3N.shape[0], products6_H2O.shape[0], products6_CO2.shape[0], products6_CH2O.shape[0], products6_CHNO.shape[0]]               ##  to be generalized... xx     
##backward = [products6_backward_H2.shape[0], products6_backward_O.shape[0], products6_backward_CO.shape[0], products6_backward_H3N.shape[0], products6_backward_CHN.shape[0], products6_backward_CH3N.shape[0], products6_backward_H2O.shape[0], products6_backward_CO2.shape[0], products6_backward_CH2O.shape[0], products6_backward_CHNO.shape[0]]
##df = pd.DataFrame({'Forward': forward,  'Backward': backward}, index=trafo_names)

##fontsize = 12
###plt.figure()
###df.plot.bar(rot=0, stacked=True)
###df.plot(kind='bar', rot=0,  width=.8, color={"Forward+Backward": "lightgray", "Forward": "black", 'Backward': 'w'},  edgecolor='black')
##df.plot(kind='bar', rot=0,  width=.8, color={"Forward": "black", 'Backward': 'w'},  edgecolor='black', stacked=True, fontsize=fontsize)
##plt.title(SAMPLE_START+r'$\rightarrow$'+SAMPLE_END+'    ('+str(name2)+', Component ' +str(FILTER_COMPONENTS)+')', fontsize=fontsize+2)
###plt.xticks()
###plt.yticks()
##plt.xlabel('Transformations', fontsize=fontsize+2)
##plt.ylabel('Frequency', fontsize=fontsize+2)
##plt.legend(fontsize=fontsize)
##plt.savefig(name+'_Component'+str(FILTER_COMPONENTS)+'_'+SAMPLE_START+SAMPLE_END+'    trafo_bar_chart'+'.png')
####plt.show()













##xEDGES['xxxx'] = xEDGES['Source']*1.123456 + xEDGES['Target']*1.123456








#xNODES.plot.scatter('id','node_size')
#plt.show()

##x = x
##x.to_csv('xxxxxxxxxxxxxxx.txt', header="'\t'.join(list(x.columns))", index=None, sep='\t')




##xEXPORT.to_csv(name+'_Comp'+str(FILTER_COMPONENTS)+'    C'+str(C_min)+'-C'+str(C_max)+'    '+str(trafos)+'    '+str(no_generations)+' generations    max_path_length '+str(max_path_length)+'.txt', sep='\t', index=False)

print('Time code = '+str(datetime.datetime.now() - begin_time)+' [h:min:sec]')





