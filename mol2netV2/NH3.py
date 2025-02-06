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
#####################PREAMBLE
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
trafos          =       ['H2', 'O','CO'] #modified e-mail/why? #['H2','O','CO','NH3','CHN','H2O','CO2','CH2O','CHNO']
##trafoWeights    =       [ 1,   5,   1,   2]
trafoWeights    =       [1]*len(trafos)
columnNames     =       ['molecular_formula','abundance_int']             ##  for  molFormula  &  abundance/intensity
sample          =       'sample'                   #previously 'Ice composition H2O_CH3OH_NH3'
nComponents     =       100
only_uniques    =       'yes'
color_type      =       sample
node_size       =       'Degree'
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
files = glob.glob(os.path.join("*.csv")) #.csv
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
    x2['FILTER'] =  1 * ( (mu >= 0) )
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
    FILTER_COMPONENTS  =        list(df_comp_freq_comb.index)[:nComponents] #Z: you can select manually, we will automate this later
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


###########################################
####                    DEGREE DISTRIBUTION
###########################################
x_name = 'Degree'
y_name = 'Degree counts'
text   = 'id'
x = xNODES[x_name].tolist()
y = xNODES[y_name].tolist()
text = xNODES[text].tolist()
###
plt.figure()
#plt.title(name+'\n, Comp'+str(FILTER_COMPONENTS)+', Trafos'+str(trafos), wrap=True, fontsize=12)
plt.scatter(x,y,c='blue')
for i in range(len(x)):
    #plt.annotate(text[i], (x[i], y[i] + 0.2), fontsize=10)

    plt.xscale('log',base=10) 
    plt.yscale('log',base=10) 
plt.xlabel(x_name)
plt.ylabel(y_name)
#plt.savefig(name+'_'+str(FILTER_COMPONENTS)+'_'+str(trafos)+'    degreeDistri'+'.png')
plt.show()



#################################
####          COMPONENT HISTOGRAM
#################################
plt.figure()
pd.DataFrame(xNODES0['component'].value_counts(normalize=True)).plot(kind='bar', rot=0, ylabel='', legend=False, color='black', width=1)
#pd.DataFrame(xNODES0['component'].value_counts(normalize=False)).sort_index().plot(kind='bar', rot=0, ylabel='', legend=False, color='black', width=1)
#plt.title(name+', Component '+str(FILTER_COMPONENTS), wrap=True)     
#plt.title('Minimal set [H2, O, CO, NH3]')
plt.xticks(np.arange(0, len(df_comp_freq), 20), fontsize=15)                      ## ice
####plt.xticks(np.arange(0, len(df_comp_freq), 100), fontsize=15)                      ## paris
plt.yticks(fontsize=15)
plt.yscale('log')
plt.xlabel('Component number', fontsize=15)
###plt.ylabel('Frequency (log scale)')
plt.ylabel('Frequency (log scale, normalized)', fontsize=15)
###plt.ylabel('Frequency (normalized)')
#plt.savefig(name+', Comp'+str(FILTER_COMPONENTS)+ '    hist_components'+'.svg')
#plt.savefig(name+', Comp'+str(FILTER_COMPONENTS)+ '    hist_components'+'.png')

plt.show()
#Normalized color gradient 

# Define the color map
cmap = cm.get_cmap('viridis')  # You can choose other colormaps from matplotlib

# Normalize component values to the range [0, 1]
norm = plt.Normalize(xNODES['component'].min(), xNODES['component'].max())

# Create a scalar mappable for color mapping
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Create a list of colors for each node based on the 'component' column
node_colors_xNODES = [cmap(norm(component)) for component in xNODES['component']]

# Plot the graph with node colors
plt.figure()
plt.title("Network representation", wrap=True, fontsize=12)
nx.draw(G, pos=nx.nx_pydot.graphviz_layout(G), with_labels=False, node_color=node_colors_xNODES, node_size=2, width=.01, alpha=1)

# Plot the colorbar
plt.colorbar(sm, label='Component')

plt.text(.6, -1.1, f"{G.number_of_nodes()} nodes, {G.number_of_edges()} edges", fontsize=8, wrap=True)

plt.savefig('NetworkRepresentation_clean_NH3'+'.png')

plt.show()

# List of elements for visualization
l = ['H', 'C', 'N', 'O', 'HC', 'NC', 'OC', 'Mass (exact)', 'DBE']

# Create an empty DataFrame to store the sum of standard deviations for each subgraph
sum_std_df = pd.DataFrame(index=xNODES['df_comp_freq'].unique(), columns=l)

# Iterate over elements for visualization
for color_type in l:
    carac = xNODES[['id', color_type]]
    carac = carac.set_index('id')
    carac = carac.reindex(G.nodes())
    carac[color_type] = pd.Categorical(carac[color_type])
    carac[color_type].cat.codes
    nodes = G.nodes() 
    #in the code above, we just color nodes based on the value of the element in question (in l)

    plt.figure()
    plt.title(color_type, wrap=True, fontsize=25)

    pos = nx.nx_pydot.graphviz_layout(G)
    nx.draw(G, pos, with_labels=False, node_color=node_colors, node_size=2, width=.01, alpha=1)
    ec = nx.draw_networkx_edges(G, pos, width=.01, alpha=0.2)
    nc = nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=carac[color_type], node_size=3, cmap=plt.cm.jet)

    cb = plt.colorbar(nc, orientation='vertical')
    nc.figure.axes[0].tick_params(axis="both", labelsize=21)
    nc.figure.axes[1].tick_params(axis="y", labelsize=21)
    plt.axis('off')
    plt.savefig('fCSoff_Supporting Figure 7 - ' + color_type + '_minimal_set_NH3_all_components' + '.pdf')
    plt.show()

    # Calculate the standard deviation of H for each subgraph
    subgraph_std = xNODES.groupby('df_comp_freq')[color_type].std()

    # Normalize the standard deviations to use for coloring
    normalized_std = (subgraph_std - subgraph_std.min()) / (subgraph_std.max() - subgraph_std.min())

    # Create a dictionary to map subgraphs to their normalized standard deviations
    subgraph_color_dict = dict(zip(normalized_std.index, normalized_std.values))

    # Map the normalized standard deviations to each subgraph in the graph
    node_colors = [subgraph_color_dict.get(df_comp_freq, 0) for df_comp_freq in xNODES['df_comp_freq']]

    # Store the normalized standard deviations in the sum_std_df DataFrame
    sum_std_df[color_type] = normalized_std

# Sum the standard deviations for each subgraph across all elements
sum_std_df['sum_std'] = sum_std_df.sum(axis=1)

# Plot the graph with the new colormap based on the sum of standard deviations
plt.figure()
plt.title('Sum of Standard Deviations in Each Subgraph', wrap=True, fontsize=25)

pos = nx.nx_pydot.graphviz_layout(G)
nx.draw(G, pos, with_labels=False, node_color=node_colors, node_size=3, width=.01, alpha=1)
ec = nx.draw_networkx_edges(G, pos, width=.01, alpha=0.2)
nc = nx.draw_networkx_nodes(G, pos, nodelist=G.nodes(), node_color=node_colors, node_size=3, cmap=plt.cm.jet)

# Add colorbar
cb = plt.colorbar(nc, orientation='vertical')
nc.figure.axes[0].tick_params(axis="both", labelsize=21)
nc.figure.axes[1].tick_params(axis="y", labelsize=21)

plt.axis('off')
plt.savefig('fCSoff_Sum_of_Standard_Deviations_NH3.pdf')
plt.show()
#top 25-percentile code
n=len(subgraph_std)
n2=int(np.around(n/4))
size_subgraph_std=subgraph_std.tail(n2)

xNODES['colored_subgraph'] = -1  # Default value, you can change this to any other default value
xNODES['colored_subgraph_component'] = -1  # Default value, you can change this to any other default value

# Find the 5 subgraphs with the smallest standard deviation
top5_subgraphs = size_subgraph_std.nsmallest(5).index

# Create a dictionary to map subgraphs to their index for the top 5
top5_subgraph_dict = dict(zip(top5_subgraphs, range(5)))

# Map the top 5 subgraphs to orange color and update 'colored_subgraph' and 'colored_subgraph_component'
node_colors_top5 = ['orange' if df_comp_freq in top5_subgraphs else 'gray' for df_comp_freq in xNODES['df_comp_freq']]
xNODES.loc[xNODES['df_comp_freq'].isin(top5_subgraphs), 'colored_subgraph'] = 'orange'
xNODES.loc[xNODES['df_comp_freq'].isin(top5_subgraphs), 'colored_subgraph_component'] = xNODES['component']

# Plot the graph with the top 5 subgraphs in orange
plt.figure()
#plt.title('Top 5 Subgraphs with Smallest Standard Deviation', wrap=True, fontsize=25)

pos = nx.nx_pydot.graphviz_layout(G)
nx.draw(G, pos, with_labels=False, node_color=node_colors_top5, node_size=3, width=.01, alpha=1)
ec = nx.draw_networkx_edges(G, pos, width=.01, alpha=0.2)
nc = nx.draw_networkx_nodes(G, pos, nodelist=G.nodes(), node_color=node_colors_top5, node_size=3, cmap=plt.cm.jet)

# Add colorbar
cb = plt.colorbar(nc, orientation='vertical')
nc.figure.axes[0].tick_params(axis="both", labelsize=21)  # Change label size
nc.figure.axes[1].tick_params(axis="y", labelsize=21)  # Change tick label size of colorbar

plt.axis('off')
plt.savefig('fCSoff_Top_5_Subgraphs_with_Smallest_Std_NH3.png')
plt.show()
# Initialize 'colored_subgraph' column with 'gray'
xNODES['colored_subgraph'] = 'gray'

# Find the 5 subgraphs with the smallest standard deviation
top5_subgraphs = size_subgraph_std.nsmallest(5).index

# Define shades of orange for the top 5 subgraphs
#shades_of_orange = ['#FFA500', '#FF8C00', '#FFA07A', '#FF7F50', '#FFD700']
shades_of_orange = ['#DFAA4A', '#FF8C00', '#FFA07A', '#FF7F50', '#FFD700']

#

# Assign each component a different shade of orange
for i, comp in enumerate(top5_subgraphs):
    # Calculate the color for this component
    color = shades_of_orange[i]
    
    # Update 'colored_subgraph' column for nodes in this component to the calculated color
    xNODES.loc[xNODES['df_comp_freq'] == comp, 'colored_subgraph'] = color

# Plot the graph with the top 5 subgraphs in different shades of orange
plt.figure()
pos = nx.nx_pydot.graphviz_layout(G)
nx.draw(G, pos, with_labels=False, node_color=xNODES['colored_subgraph'], node_size=3, width=.01, alpha=1)
plt.axis('off')
plt.savefig('fCSoff_Top_5_Subgraphs_with_Smallest_Std_NH3.png')
plt.show()
fCxNODES_annotated = xNODES.copy()
fCxNODES_annotated.to_csv("fCxNODES_annotated.tsv")
# Get unique components with color 'gray'
gray_components = set(xNODES[xNODES['colored_subgraph'] == 'gray']['component'].unique())

# Filter out gray components from FILTER_COMPONENTS
FILTER_COMPONENTS_filtered = np.array([comp for comp in FILTER_COMPONENTS if comp not in gray_components])

# Now FILTER_COMPONENTS_filtered contains only the IDs of components that do not have the color 'gray' under 'colored_subgraph'
plt.figure()

# Filter nodes based on FILTER_COMPONENTS_filtered
filtered_nodes = xNODES[xNODES['component'].isin(FILTER_COMPONENTS_filtered)]

# Filter colors for the filtered nodes
filtered_colors = []
for _, node_data in filtered_nodes.iterrows():
    node_id = node_data['id']
    color = node_data['colored_subgraph']
    if color != 'gray':
        filtered_colors.append(color)
    else:
        filtered_colors.append('#808080')  # Use #808080 for gray nodes

# Draw the graph with filtered nodes and colors
nx.draw(G, pos=nx.nx_pydot.graphviz_layout(G), nodelist=filtered_nodes['id'], with_labels=False, node_color=filtered_colors, node_size=1, width=.05, alpha=1)

# Display statistics next to each component
for comp in FILTER_COMPONENTS_filtered:
    # Filter data for the current component
    component_data = filtered_nodes[filtered_nodes['component'] == comp]
    
    # Find a node in the component and get its position
    node_id = component_data.iloc[0]['id']
    x, y = pos[node_id]

    # Calculate statistics for the component
    stats = {
        'C_min': component_data['C'].min(),
        'C_max': component_data['C'].max(),
        'H_min': component_data['H'].min(),
        'H_max': component_data['H'].max(),
        'N_min': component_data['N'].min(),
        'N_max': component_data['N'].max(),
        'O_min': component_data['O'].min(),
        'O_max': component_data['O'].max(),
        'S_min': component_data['S'].min(),
        'S_max': component_data['S'].max(),
        'mz_range': f"{component_data['Mass (exact)'].min():.2f}-{component_data['Mass (exact)'].max():.2f}",
        'DBE_range': f"{component_data['DBE'].min():.2f}-{component_data['DBE'].max():.2f}",
        'num_nodes': len(component_data),
    }

    text = f"Component: {comp}\n\
            C({stats['C_min']} - {stats['C_max']})H({stats['H_min']} - {stats['H_max']})N({stats['N_min']} - {stats['N_max']})O({stats['O_min']} - {stats['O_max']})\n\
            S: {stats['S_min']} - {stats['S_max']}\n\
            mz range: {stats['mz_range']}\n\
            DBE range: {stats['DBE_range']}\n\
            Mol. For.: {stats['num_nodes']}"

    plt.text(x-1000, y+100, text, fontsize=8, ha='left', va='center') #play with these to change the position of the text

# Other plot settings
for i in range(len(sample_colors)):
    plt.plot([], [], sample_colors.values[i][1], marker='o', markersize=10, label=sample_colors.values[i][0])

plt.legend()
plt.savefig('fc_S8_CO_Annotated_Soff_Range_Filtered_Filtered_v2' +'.pdf')
plt.show()
