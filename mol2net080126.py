'''
###################################################
#               mol2net                           #
#                                                 #
#    Network analysis from a list of molecules    #
###################################################   

AUTHOR: Alexander Ruf (rufalexan@gmail.com)

If you use this code for your work, please cite the 
corresponding paper:

-- 10.1021/acs.analchem.2c01271

-- the github repository github.com/rufalexan/mol2net 

-- corresponding Zenodo DOI:
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
import time
begin_time = datetime.datetime.now()
import os
import glob
import string
import chemparse
import re
import numpy as np
import pandas as pd
from collections import namedtuple, defaultdict
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
__all__ = ['nx', 'drawNet']   ## allow using nx in draw function
import pydot
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.ticker as ticker
from mycolorpy import colorlist as mcp
import math
from datetime import datetime
now = datetime.now()
date = now.strftime("%Y%m%d %H%M%S") + f".{int(now.microsecond / 1000):03d}"
from tqdm import tqdm
m_h, m_c, m_n, m_o, m_s = 1.007825, 12., 14.003074, 15.994915, 31.972071
masses_elements = (m_h, m_c, m_n, m_o, m_s)
valence_h, valence_c, valence_n, valence_o, valence_s = 1, 4, 3, 2, 2
mpl.rcParams['font.size'] = 12
import seaborn as sns


  

def molForm2graph(  fileName,molFormula,abundance , trafos   ,  comps=None,sample=None,trafoWeight=None  ):     # _opt  trafoWeight = [1,5,1,2]
    fnName = '_molForm2graph'
    import time
    start_time = time.time()
    print(f"\n\n{'#'*60}")
    print(f"==== {fnName} ====  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    #
    global xTRAFOS, mdTRAFOS, name, orig , xNODES0, xDET, xEDGES0, G0 , df_comp_freq, xEDGES, xNODES, FILTER_COMPONENTS, md_freq, G, sample2, trafos2, comp_diversity_pivot, molFormula2, comps2, m_xc
    #
    if comps is None:
        comps = list(range(10))
    else:
        comps = comps
    #
    comps2 = comps
    #
    trafos2 = trafos
    molFormula2 = molFormula
    if trafoWeight is None:
            trafoWeight = [1] * len(trafos)
    #
    def get_node_sizes(xNODES, scale_by="degree", nodeScale=1.0, base_size=300):
        if scale_by not in xNODES.columns:
            raise ValueError(f"Column {scale_by} not in xNODES")
            #
        values = xNODES[scale_by].fillna(0).values.astype(float)
        if values.max() != values.min():
            norm = (values - values.min()) / (values.max() - values.min())
        else:
            norm = np.ones_like(values)
            #
        return (norm * base_size * nodeScale).tolist()
        #
    globals()['get_node_sizes'] = get_node_sizes
    #####################################################################################  _trafos
    def parse_formula(formula):
        def parse_group(subformula, multiplier=1):
            counts = defaultdict(int)
            matches = re.findall(r'([A-Z][a-z]?)(-?\d*)', subformula)
            for element, count in matches:
                count = int(count) if count else 1
                counts[element] += count * multiplier
            return counts
        #
        stack = []
        element_counts = defaultdict(int)
        i = 0
        while i < len(formula):
            if formula[i] == '(':
                stack.append(i)  # Store index of '('
            elif formula[i] == ')':
                start = stack.pop()
                i += 1
                multiplier = 1
                num_match = re.match(r'-?\d+', formula[i:])
                if num_match:
                    multiplier = int(num_match.group())
                    i += len(num_match.group())
                group_counts = parse_group(formula[start + 1:i - 1], multiplier)
                for element, count in group_counts.items():
                    element_counts[element] += count
            else:
                num_match = re.match(r'([A-Z][a-z]?)(-?\d*)', formula[i:])
                if num_match:
                    element, count = num_match.groups()
                    count = int(count) if count else 1
                    element_counts[element] += count
                    i += len(num_match.group()) - 1
            i += 1
        #
        return dict(element_counts)
    #
    def parse_multiple_trafos(trafos):
        data = []
        for formula in trafos:
            element_counts = parse_formula(formula)
            row = {"Formula": formula, **element_counts}
            data.append(row)
        return pd.DataFrame(data).fillna(0).astype(dtype=int, copy=False, errors='ignore')
    #
    x2 = parse_multiple_trafos(trafos)
    x2 = x2.rename(columns={'Formula': 'Molecular Formula'})
    if 'S' not in x2.columns:
        x2['S'] = 0
    #
    x3 = x2.drop(['Molecular Formula'], axis=1)
    for i in x3.columns.to_list():
        if i in x2.columns:
            x2[i] = x2[i].astype(int)
        else:
            x2[i] = 0    
    #
        globals()[f"{i.lower()}"] = x2[i].to_numpy()
    #
    for elem in ['H','C','N','O','S']:
        if elem not in x2.columns:
            x2[elem] = 0  # fill missing elements with 0
    #
    elements = x2[['H', 'C', 'N', 'O', 'S']].to_numpy(dtype=float)   # shape (N, 5)
    massExact = elements @ np.array(masses_elements)                 # shape (N,)
    x2['Mass (exact)'] = massExact.round(6)
    #
    weights = pd.DataFrame(trafoWeight,columns = ['trafoWeight'])
    weights['trafos'] = trafos
    x3 = x2.join(weights.set_index(['trafos'], verify_integrity=True ),on=['Molecular Formula'], how='left')
    xTRAFOS = x3
    mdTRAFOS = xTRAFOS['Mass (exact)'].round(6).to_numpy()
    #####################################################################################  _nodes
    path = os.getcwd()
    files = glob.glob(os.path.join(fileName))
    d = {}
    for f in files:
        key = f
        name = f.replace(".dat", "")
        x = orig = globals()[f"orig_{name}"] = pd.read_csv(f, sep='\t', skiprows=0, encoding = 'unicode_escape')
        #
        if abundance not in x.columns:
            x[abundance] = 1
        #
        if 'id' in x.columns:
            pass
        else:
            x.insert(0, 'id', range(1, 1 + len(x)))
            #    
        try:
            x = x[['id',sample,molFormula,abundance]].rename(columns={molFormula: 'Molecular Formula',abundance:'Abundance'})
        except:
            x = x[['id',molFormula,abundance]].rename(columns={molFormula: 'Molecular Formula',abundance:'Abundance'})
            #
        molFormulas = x['Molecular Formula']
        a = []
        for i in molFormulas:
            tmp_a = chemparse.parse_formula(i)
            a.append(tmp_a)
        #
        x_a = pd.DataFrame(a).fillna(0)
        if 'S' not in x_a.columns:
            x_a['S'] = 0
        #
        x2 = pd.concat([x,x_a], axis=1)
        #
        for i in x_a.columns.to_list():
            if i in x2.columns:
                x2[i] = x2[i].astype(int)
            else:
                x2[i] = 0    
        #
            globals()[f"{i.lower()}"] = x2[i].to_numpy()
        #
        elements = np.transpose(np.array([h, c, n, o, s]))
        massExact = np.dot(elements, masses_elements)
        #
        def assign_chemical_family(row):
            family = ''
            for elem in ['C', 'H', 'N', 'O', 'S']:
                if row.get(elem, 0) > 0:
                    family += elem
            return family
        x2['chemicalFamily'] = x2.apply(assign_chemical_family, axis=1)
        #
        x2['Mass (exact)'] = massExact.round(6)
        x2['mu'] = mu = ( (h*valence_h + c*valence_c + n*valence_n + o*valence_o + s*valence_s) / 2 ) - (h + c + o + n + s) + 1
        x2['DBE'] = dbe = c - (h/2) + (n/2) + 1
        x2['H/C'] = h/c
        x2['N/C'] = n/c
        x2['O/C'] = o/c
        x2['S/C'] = s/c
        x2['S/O'] = np.nan_to_num(s / np.where(o == 0, np.nan, o), nan=0)
        m_xc = .5         # m is the fraction of oxygen in an astrophysical ice residue and is estimated to be 0.5 for similar experiments (Danger et al. 2016)  __x
        x2[f'X_C(m={m_xc})'] = pd.Series((2*c + n - h - 2*m_xc*o) / np.where(np.abs(dbe - m_xc*o) < 1e-10, np.nan, dbe - m_xc*o) + 1).fillna(0)
        x2['AMD'] = abs(massExact - massExact.round(0)) 
        x2['KMD'] = (massExact * 14/14.01565).round(0) - (massExact * 14/14.01565)
        x2['FILTER'] =  1 * ( (mu >= 0) )
        x2['log(Abundance)'] = np.log10(x2['Abundance'])
        atom_list = ['H','C','N','O','S']
        for i in range(len(atom_list)-1):
            x2['Molecular Formula'] = x2['Molecular Formula'].replace({atom_list[i]+'0'}, {atom_list[i]}, regex=True)
            x2['Molecular Formula'] = x2['Molecular Formula'].replace({atom_list[i]+'1'+atom_list[i+1]}, {atom_list[i]}, regex=True)
            x2['Molecular Formula'] = x2['Molecular Formula'].replace({'N'+'0'}, {''}, regex=True)
            x2['Molecular Formula'] = x2['Molecular Formula'].replace({'N'+'1'}, {'N'}, regex=True)
            x2['Molecular Formula'] = x2['Molecular Formula'].replace({'S'+'0'}, {''}, regex=True)
            x2['Molecular Formula'] = x2['Molecular Formula'].replace({'S'+'1'}, {'S'}, regex=True)
        #
        x2['Mass counts'] = x2['Mass (exact)'].map(dict(x2['Mass (exact)'].value_counts()))
        try:
            x3 = x2.drop_duplicates(subset=['Mass (exact)'], keep='first')
        except:
            x3 = x2
        #
        if sample is not None and sample in x2.columns:     # Assign sample (as color label) depending on whether 'sample' is provided
            sample = sample
        else:
            sample = 'chemicalFamily'
        #
        sample2 = sample
        #
        x3 = x3[(x3['FILTER']==1)]
        xNODES0 = xDET = x3
        globals()[f"{name}"] = xDET
        mass = xDET['Mass (exact)'].round(6).to_numpy()
        data = xDET[['id','Mass (exact)']]
        #####################################################################################  _edges
        md_matches = namedtuple('md_matches', 'md_matches hits')                
        new_list = []
        for md in mdTRAFOS:
            for element in mass:
                if element+md in mass:
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
        xxxx = xTRAFOS[['Mass (exact)','trafoWeight']]
        x2 = x.join(xxxx.set_index(['Mass (exact)'], verify_integrity=True ),on=['Mass difference'], how='left').rename(columns={'trafoWeight': 'Weight'}).sort_values(by=['Source'])
        #
        xEDGES0 = x2
        ###################################################################################  _networkAna 
        G0 = nx.from_pandas_edgelist(xEDGES0, 'Source', 'Target', create_using=nx.Graph())
        nodes_in_graph = set(G0.nodes)
        xNODES0 = xNODES0[xNODES0["id"].isin(nodes_in_graph)].copy()
        #
        components = sorted(nx.connected_components(G0), key=len, reverse=True)         # FILTER G0 --> G
        node_to_comp = {}
        for comp_id, comp_nodes in enumerate(components):
            for node in comp_nodes:
                node_to_comp[node] = comp_id
        #
        xNODES0["Component"] = xNODES0["id"].map(node_to_comp)   # Attach that to the node dataframe (keeps everything in one place)
        FILTER_COMPONENTS = comps             # keeps these comps, ascending numbers sorted by size
        #
        nodes_mask = xNODES0["Component"].isin(FILTER_COMPONENTS)
        xNODES   = xNODES0.loc[nodes_mask]
        node_id_to_comp = xNODES.set_index("id")["Component"]    # Build a Series mapping node id --> Component
        xNODES = xNODES.copy()
        xNODES['Component counts'] = xNODES["Component"].map(dict(xNODES["Component"].value_counts()))
        #
        keeper_nodes = set(xNODES["id"])
        G = G0.subgraph(keeper_nodes).copy()
        degree_dict = dict(G.degree())
        #
        xNODES['degree'] = xNODES['id'].map(degree_dict)
        xNODES['degree counts'] = xNODES['degree'].map(xNODES['degree'].value_counts())
        #
        edges_mask = ( xEDGES0["Source"].isin(keeper_nodes) & xEDGES0["Target"].isin(keeper_nodes) )
        xEDGES   = xEDGES0.loc[edges_mask]
        xEDGES = xEDGES.copy()
        xEDGES["Component_source"] = xEDGES["Source"].map(node_id_to_comp)   # Map to edges
        xEDGES["Component_target"] = xEDGES["Target"].map(node_id_to_comp)
        xEDGES = xEDGES.rename(columns={'Molecular Formula': 'Transformation', 'Mass difference': 'Transformation mass'})
        xEDGES['Transformation counts'] = xEDGES["Transformation"].map(dict(xEDGES["Transformation"].value_counts()))
        xEDGES['degree_source'] = xEDGES['Source'].map(degree_dict)
        xEDGES['degree_target'] = xEDGES['Target'].map(degree_dict)
        xEDGES['degree_source counts'] = xEDGES['degree_source'].map(xEDGES['degree_source'].value_counts())
        xEDGES['degree_target counts'] = xEDGES['degree_target'].map(xEDGES['degree_target'].value_counts())
        #
        comp_counts = xNODES0['Component'].value_counts()
        comp_rank = {comp: rank+0 for rank, comp in enumerate(comp_counts.index)}
        df_comp_freq = comp_counts.rename_axis("Component").reset_index(name="node counts")
        df_comp_freq["node counts normalized"] = df_comp_freq["node counts"] / df_comp_freq["node counts"].sum() 
        df_comp_freq["Component"] = df_comp_freq["Component"].map(comp_rank)
        #
        a = xEDGES[['Transformation','Transformation counts']]
        md_freq = a.drop_duplicates(keep='first').sort_values(by=['Transformation counts'], ascending=False)
        #
        top_nodes = xNODES[xNODES["Component"].isin(list(FILTER_COMPONENTS))]
        comp_diversity = top_nodes.groupby(["Component", sample2]).size().reset_index(name="count")
        comp_diversity_pivot = comp_diversity.pivot(index="Component", columns=sample2, values="count").fillna(0).astype(int)
        components_nonzero = comp_diversity_pivot[(comp_diversity_pivot != 0).all(axis=1)].index.to_list()        
        #
        print(f"{name}\nSamples: {xNODES[sample].unique().tolist()}\ntrafos: {trafos2}")
        print(f"G0: {len(xNODES0)} nodes, {len(xEDGES0)} edges ({len(df_comp_freq)} (all) components)")
        print(f"G: {len(xNODES)} nodes, {len(xEDGES)} edges (component(s) {FILTER_COMPONENTS})")
        print('-'*60)
        print(df_comp_freq)
        print('-'*60)
        print(md_freq)
        print('                                       (For component(s)' +str(FILTER_COMPONENTS) +')')
        print('-'*60)
        print(comp_diversity_pivot)
        print("Components with maximum sample diversity:", components_nonzero)
        print('#'*60)
        print('\n'*4)






def drawNet(nodeColors=None, layout=nx.nx_pydot.graphviz_layout, nodeSize=1.0,
            label=False, labelSize=2, edgeColors=None, edgeWidth=1.0,
            nodeAlpha=0.5, edgeAlpha=0.8):
    """
    Draw the molecular network graph using global mol2net variables.
    Automatically colors nodes based on `sample2` (e.g., chemicalFamily)
    using a predefined palette for known families, and a contrasting colormap
    for unknown ones.
    """
    import matplotlib.pyplot as plt
    import networkx as nx
    import time
    from datetime import datetime
    import pandas as pd
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors

    fnName = '_drawNet'
    start_time = time.time()
    print('#' * 60)
    print(f"==== {fnName} ====  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    global G, xNODES, sample2, FILTER_COMPONENTS, trafos2, xEDGES, name, date
    default_color = '#e0e0e0'

    # ============================
    # NODE COLORS
    # ============================
    if nodeColors is None:
        unique_samples = sorted(xNODES[sample2].unique())

        # 1️⃣ Predefined color mapping for known families
        predefined_colors = {
            'CHO': '#2c81d1',       # blue
            'CHNO': '#ff7f0e',      # orange
            'CHNOS': 'red',     # red
            'CHOS': '#2ca02c',      # green
            'CHNS': '#bcbd22',      # yellow
            'CH': '#7f7f7f',        # gray
            'C': '#000000',         # black
            'CHS': '#98df8a',       # light green
            'CHN': '#9467bd',       # light purple
        }

        # 2️⃣ Assign contrasting colors to unknowns
        missing_samples = [s for s in unique_samples if s not in predefined_colors]
        if missing_samples:
            n_missing = len(missing_samples)
            if n_missing <= 10:
                cmap = cm.get_cmap('gist_rainbow', n_missing)     #tab10
            else:
                cmap = cm.get_cmap('gist_rainbow', n_missing)
            new_colors = [mcolors.to_hex(cmap(i)) for i in range(n_missing)]
            for s, c in zip(missing_samples, new_colors):
                predefined_colors[s] = c

        # 3️⃣ Build nodeColors DataFrame
        colors = [predefined_colors.get(s, default_color) for s in unique_samples]
        nodeColors = pd.DataFrame({sample2: unique_samples, 'color': colors})

    elif isinstance(nodeColors, list) and len(nodeColors) == 2:
        nodeColors = pd.DataFrame(list(zip(*nodeColors)), columns=[sample2, 'color'])
    else:
        unique_samples = xNODES[sample2].unique()
        colors = [default_color] * len(unique_samples)
        nodeColors = pd.DataFrame({sample2: unique_samples, 'color': colors})

    # ============================
    # NODE MAPPING
    # ============================
    node_colors, node_sizes, labels_dict = [], [], {}
    degree_max = xNODES['degree'].max() if 'degree' in xNODES.columns else 1

    for node in G.nodes():
        row = xNODES[xNODES['id'] == node]
        if not row.empty:
            sample_value = row[sample2].values[0]
            color_row = nodeColors[nodeColors[sample2] == sample_value]
            color = color_row['color'].values[0] if not color_row.empty else default_color
            node_colors.append(color)
            deg = row['degree'].values[0] if 'degree' in row.columns else 1
            node_sizes.append(max(2, deg / degree_max * 40 * nodeSize))
            if label:
                labels_dict[node] = str(row['Molecular Formula'].values[0])
        else:
            node_colors.append(default_color)
            node_sizes.append(2)
            if label:
                labels_dict[node] = str(node)

    # ============================
    # EDGE COLORS
    # ============================
    edge_color_list, edge_legend_entries = [], {}
    if edgeColors is None:
        edge_color_list = [default_color] * len(xEDGES)
    else:
        for _, row in xEDGES.iterrows():
            trafo = row['Transformation']
            color = edgeColors.get(trafo, default_color)
            edge_color_list.append(color)
            if trafo in edgeColors:
                edge_legend_entries[trafo] = color

    # ============================
    # LAYOUT
    # ============================
    pos = layout(G)
    layout_name = layout.__name__

    # ============================
    # DRAW NETWORK
    # ============================
    plt.figure(figsize=(10, 8))
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=nodeAlpha)
    nx.draw_networkx_edges(G, pos, edge_color=edge_color_list, width=edgeWidth, alpha=edgeAlpha)
    if label:
        nx.draw_networkx_labels(G, pos, labels=labels_dict, font_size=labelSize)

    # ============================
    # LEGENDS
    # ============================
    for i in range(len(nodeColors)):
        c = nodeColors.iloc[i]['color']
        if c != default_color:
            plt.plot([], [], marker='o', color=c, markersize=10,
                     label=nodeColors.iloc[i][sample2])

    if edgeColors is not None:
        for trafo, color in edge_legend_entries.items():
            if color != default_color:
                plt.plot([], [], color=color, linewidth=2, label=trafo)

    plt.legend(loc='upper left', bbox_to_anchor=(-0.02, 1.02), fontsize=10)
    plt.axis('off')
    plt.subplots_adjust(top=0.92, bottom=0.05)

    # ============================
    # METADATA TEXT
    # ============================
    plt.text(
        0.98, 0.02,
        f"{name}\nG: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges\n"
        f"component(s): {FILTER_COMPONENTS}\ntrafos: {trafos2},\nlayout: {layout_name}",
        fontsize=10, wrap=True, ha='right', va='bottom', transform=plt.gcf().transFigure
    )

    # ============================
    # SAVE & REPORT
    # ============================
    outfile = f"{date}  {name}  Comp{FILTER_COMPONENTS}  {trafos2}  drawNet.png"
    plt.savefig(outfile, bbox_inches='tight', dpi=200)

    elapsed_total = time.time() - start_time
    print(f"{fnName} finished in {elapsed_total:.2f}s, saved as {outfile}")
    print(f"Time per node: {elapsed_total / G.number_of_nodes():.6f} s/node")
    print(f"Time per edge: {elapsed_total / G.number_of_edges():.6f} s/edge")
    print('#' * 60)
    print('\n' * 4)











def compDistri(G, layout=nx.nx_pydot.graphviz_layout, nodeScale=1.0):
    """
    Plot component distribution and network layout with node sizes scaled consistently.
    """
    fnName = '_compDistri'
    import time
    start_time = time.time()
    print(f"\n==== {fnName} ====  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    plt.figure(figsize=(11,5))

    # --- Plot 1: Component distribution ---
    plt.subplot(1,2,1)                                                  
    plt.bar(df_comp_freq['Component'].to_numpy(), 
            df_comp_freq['node counts normalized'].to_numpy(), 
            color="black")
    plt.yscale('log', base=10) 
    plt.xlabel('Component', fontsize=20)
    plt.ylabel('# Nodes (normalized)', fontsize=20)

    # --- Plot 2: Network layout ---
    plt.subplot(1,2,2)
    fontsize=20
    color_label = 'Component'
    carac = xNODES0[['id', color_label]].set_index('id')
    carac = carac.reindex(G.nodes())
    carac[color_label] = pd.Categorical(carac[color_label])

    nodes = G.nodes()
    pos = layout(G)
    layout_name = layout.__name__

    # --- Unified node scaling using get_node_sizes ---
    xNODES_plot = xNODES0.set_index('id').loc[list(nodes)].copy()
    xNODES_plot['degree'] = [G.degree(n) for n in nodes]  # optional: include degree if needed
    node_sizes = get_node_sizes(
        xNODES_plot, 
        scale_by='degree',       # you can also use color_label if you want color categories to scale
        nodeScale=nodeScale, 
        base_size=300
    )

    nx.draw_networkx_edges(G, pos=pos, alpha=0.2, width=0.1)
    nx.draw_networkx_nodes(
        G, pos=pos, nodelist=nodes,
        node_color=carac[color_label].cat.codes,
        cmap=plt.cm.jet,
        node_size=node_sizes
    )

    cb = plt.colorbar(
        plt.cm.ScalarMappable(
            cmap=plt.cm.jet,
            norm=plt.Normalize(vmin=carac[color_label].cat.codes.min(), vmax=carac[color_label].cat.codes.max())
        ),
        orientation='vertical'
    )
    cb.set_label(label=color_label, size=fontsize+2)
    plt.axis('off')

    plt.tight_layout()
    plt.text(
        0.75, 0.1, 
        f"{name}\n{G.number_of_nodes()} nodes, {G.number_of_edges()} edges,\nlayout: {layout_name},\ntrafos: {trafos2}",
        fontsize=8, wrap=True, transform=plt.gcf().transFigure
    )
    plt.savefig(f"{date}  {name}  Comp{FILTER_COMPONENTS}  {fnName}   {trafos2} {sample2}.png")
#    plt.show()







def compDistri2(nodeColors=None):
    """
    Plot component distributions and sample diversity using the same color scheme as drawNet.
    Ensures correct color mapping for chemical families.
    """
    import time
    import matplotlib.pyplot as plt
    import seaborn as sns
    from datetime import datetime
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    import pandas as pd

    fnName = '_compDistri2'
    start_time = time.time()
    print(f"\n\n{'#'*60}")
    print(f"==== {fnName} ====  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    print('#'*60)

    default_color = '#e0e0e0'

    # ============================
    # DEFINE COLORS
    # ============================
    unique_samples = sorted(xNODES[sample2].unique())

    # 1️⃣ Predefined color mapping for known families
    predefined_colors = {
        'CHO': '#2c81d1',       # blue
        'CHNO': '#ff7f0e',      # orange
        'CHNOS': 'red',     # red
        'CHOS': '#2ca02c',      # green
        'CHNS': '#bcbd22',      # yellow
        'CH': '#7f7f7f',        # gray
        'C': '#000000',         # black
        'CHS': '#98df8a',       # light green
        'CHN': '#9467bd',       # light purple
    }

    # Assign contrasting colors to unknown families
    missing_samples = [s for s in unique_samples if s not in predefined_colors]
    if missing_samples:
        n_missing = len(missing_samples)
        if n_missing <= 10:
            cmap = cm.get_cmap('tab10', n_missing)
        else:
            cmap = cm.get_cmap('rainbow', n_missing)
        new_colors = [mcolors.to_hex(cmap(i)) for i in range(n_missing)]
        for s, c in zip(missing_samples, new_colors):
            predefined_colors[s] = c

    # ============================
    # BUILD nodeColors DF
    # ============================
    colors = [predefined_colors.get(s, default_color) for s in unique_samples]
    nodeColors = pd.DataFrame({sample2: unique_samples, 'color': colors})

    # ============================
    # PREPARE PALETTE FOR SEABORN
    # ============================
    # Get unique samples in the actual melted data
    df_melted = comp_diversity_pivot.reset_index().melt(
        id_vars='Component', var_name=sample2, value_name='Count'
    )
    unique_samples_in_data = sorted(df_melted[sample2].unique())
    palette = {s: predefined_colors.get(s, default_color) for s in unique_samples_in_data}

    # ============================
    # PLOTTING
    # ============================
    plt.figure(figsize=(13,6))

    # Plot 1: Component distribution (log scale)
    plt.subplot(1,2,1)
    plt.bar(df_comp_freq['Component'].to_numpy(), 
            df_comp_freq['node counts normalized'].to_numpy(), 
            color="black")
    plt.yscale('log', base=10) 
    plt.xlabel('Component', fontsize=14)
    plt.ylabel('# Nodes (normalized)', fontsize=14)

    # Plot 2: Ion mode diversity per component
    plt.subplot(1,2,2)
    sns.barplot(
        data=df_melted, x='Component', y='Count',
        hue=sample2, palette=palette,
        hue_order=unique_samples_in_data  # <-- ensures correct mapping
    )
    plt.xlabel('Component', fontsize=14)
    plt.ylabel('# Nodes', fontsize=14)
    plt.legend(title=sample2, fontsize=10, title_fontsize=11)

    # Adjust layout and add metadata text
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2)
    plt.text(
        0.65, 0.02, 
        f"{name}\nG0: {len(xNODES0)} nodes, {len(xEDGES0)} edges ({len(df_comp_freq)} (all) components)\n"
        f"G: {len(xNODES)} nodes, {len(xEDGES)} edges (component(s) {FILTER_COMPONENTS})\n"
        f"trafos: {', '.join(trafos2)}", 
        fontsize=8, wrap=True, transform=plt.gcf().transFigure
    )

    # Save figure
    fname = f"{date}  {name}  Comp{FILTER_COMPONENTS}  {fnName}   {trafos2} {sample2}.png"
    plt.savefig(fname, dpi=300)
    elapsed_total = time.time() - start_time
    print(f"{fnName} finished in {elapsed_total:.2f}s, saved as {fname}")














def elementMaps(elements, nodeScale=1.0, nCols=4, layout=nx.nx_pydot.graphviz_layout):
    import time, math, matplotlib.pyplot as plt, networkx as nx

    fnName = '_elementMaps'
    start_time = time.time()
    print(f"\n{'#'*60}")
    print(f"==== {fnName} ====  {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}")
    daten2 = now.strftime("%Y%m%d %H%M%S") + f".{int(now.microsecond / 1000):03d}"

    global G, xNODES, FILTER_COMPONENTS, trafos2, name, date

    # Check data
    if 'G' not in globals() or G is None:
        print("Graph G not defined. Run molForm2graph first.")
        return
    if xNODES is None or xNODES.empty:
        print("xNODES missing or empty.")
        return

    # Layout and figure setup
    nRows = math.ceil(len(elements) / nCols)
    fig = plt.figure(figsize=(5 * nCols, 5 * nRows))
    pos = layout(G)
    layout_name = layout.__name__

    # Align node table with graph
    xNODES_plot = xNODES.set_index('id').reindex(G.nodes()).copy()
    if 'degree' not in xNODES_plot.columns:
        xNODES_plot['degree'] = 1

    node_sizes = get_node_sizes(
        xNODES_plot,
        scale_by='degree',
        nodeScale=nodeScale,
        base_size=10
    )

    # Draw each element panel
    for n, elem in enumerate(elements):
        ax = plt.subplot(nRows, nCols, n + 1)

        nx.draw_networkx_edges(G, pos, width=len(xNODES_plot) / 10000, alpha=0.2)

        if elem in xNODES_plot.columns:
            vals = xNODES_plot[elem].to_numpy()
        else:
            vals = [0] * len(G)

        nx.draw_networkx_nodes(
            G, pos,
            node_color=vals,
            node_size=node_sizes,
            cmap=plt.cm.jet,
            vmin=min(vals),
            vmax=max(vals)
        )

        sm = plt.cm.ScalarMappable(
            cmap=plt.cm.jet,
            norm=plt.Normalize(vmin=min(vals), vmax=max(vals))
        )
        sm.set_array([])
        plt.colorbar(sm, orientation='horizontal', ax=ax)

        plt.title(elem, fontsize=18)
        ax.set_frame_on(False)
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

    # Add metadata text
    plt.text(
        0.98, 0.02,
        f"{name}\nG: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges\n"
        f"component(s): {FILTER_COMPONENTS}\nlayout: {layout_name}\n"
        f"node size ~ degree\ntrafos: {trafos2}",
        fontsize=10, ha='right', va='bottom', transform=plt.gcf().transFigure
    )

    # Save and report
    plt.tight_layout()
    outfile = f"{date}  {name}  Comp{FILTER_COMPONENTS}   {trafos2} {fnName}.png"
    plt.savefig(outfile, bbox_inches='tight', dpi=200)
    plt.close()

    elapsed_total = time.time() - start_time
    print(f"{fnName} finished in {elapsed_total:.1f}s, saved as {outfile}")
    print(f"Time per node: {elapsed_total / G.number_of_nodes():.6f} s/node")
    print(f"Time per edge: {elapsed_total / max(1, G.number_of_edges()):.6f} s/edge")
    print('#' * 60)








def TrafoCounts(FROM_TO, SAMPLE_START, SAMPLE_END):
    import pandas as pd
    import matplotlib.pyplot as plt
    import time

    fnName = '_TrafoCounts'
    start_time = time.time()

    print(f"\n\n{'#'*60}")
    print(f"==== {fnName} ====  {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    global xNODES, xEDGES, date, name, FILTER_COMPONENTS

    # Attach source and target sample info
    df_samples = xNODES[['id', FROM_TO]]
    df_edges = xEDGES.rename(columns={'Source': 'id'})
    df_edges = df_edges.merge(df_samples, how='left', on='id')
    df_edges = df_edges.rename(columns={FROM_TO: FROM_TO + '_Source', 'id': 'Source'})

    df_edges = df_edges.rename(columns={'Target': 'id'})
    df_edges = df_edges.merge(df_samples, how='left', on='id')
    df_edges = df_edges.rename(columns={FROM_TO: FROM_TO + '_Target', 'id': 'Target'})

    # Forward and backward subsets
    df_forward = df_edges[
        (df_edges[FROM_TO + '_Source'] == SAMPLE_START) &
        (df_edges[FROM_TO + '_Target'] == SAMPLE_END)
    ]
    df_backward = df_edges[
        (df_edges[FROM_TO + '_Source'] == SAMPLE_END) &
        (df_edges[FROM_TO + '_Target'] == SAMPLE_START)
    ]

    # Extract unique trafos
    trafos = df_edges['Transformation'].unique()

    # Sort trafos by average mass difference
    avg_mass_diff = {}
    for t in trafos:
        mask = df_edges['Transformation'] == t
        diffs = df_edges.loc[mask, 'Mass (target)'] - df_edges.loc[mask, 'Mass (source)']
        avg_mass_diff[t] = diffs.mean() if len(diffs) > 0 else float('inf')

    trafos_sorted = sorted(trafos, key=lambda x: avg_mass_diff[x])

    # Count forward/backward edges per transformation
    forward_counts = [df_forward[df_forward['Transformation'] == t].shape[0] for t in trafos_sorted]
    backward_counts = [df_backward[df_backward['Transformation'] == t].shape[0] for t in trafos_sorted]

    # Data for plotting
    df_plot = pd.DataFrame({'Forward': forward_counts, 'Backward': backward_counts}, index=trafos_sorted)

    # Plot
    fontsize = 12
    df_plot.plot(
        kind='bar',
        rot=0,
        width=0.8,
        color={'Forward': 'black', 'Backward': 'white'},
        edgecolor='black',
        stacked=True,
        fontsize=fontsize
    )

    plt.title(f"{SAMPLE_START} → {SAMPLE_END}", fontsize=fontsize + 2)
    plt.xlabel('Transformation', fontsize=fontsize + 2)
    plt.ylabel('Frequency', fontsize=fontsize + 2)
    plt.legend(fontsize=fontsize)
    plt.tight_layout()

    # Save figure in mol2net style
    outfile = f"{date}  {name}  Comp{FILTER_COMPONENTS}   {trafos2} {fnName}.png"
    plt.savefig(outfile, bbox_inches='tight', dpi=300)
    plt.close()

    print(f"{fnName} finished in {time.time() - start_time:.1f}s, saved as {outfile}")
    print('#' * 60)

def trafoFreq(md_freq):             
    fnName = '_trafoFreq'
    plt.figure(figsize=(8,6))
    plt.bar(md_freq['Transformation'].to_numpy(), md_freq['Transformation (frequency)'].to_numpy(), color="black" )
    plt.xlabel('Transformation', fontsize=20)
    plt.ylabel('# transformations', fontsize=20)
    plt.legend('',frameon=False)
    plt.title(name+'\n, Comp'+str(FILTER_COMPONENTS), wrap=True, fontsize=12)
    plt.savefig(date +'  ' +name+'  Comp'+str(FILTER_COMPONENTS)+'  '+str(fnName)+'  '+str(trafos2)+'  '+str(sample2)+'.png')
#    plt.show()

def clusterTransitions(mn, sample2, sample_from, sample_to, component_source=None):
    now=datetime.now() ; date = now.strftime("%Y%m%d %H%M%S") + f".{int(now.microsecond / 1000):03d}"
    #global clusterTransitions
    node_info = mn.xNODES.set_index("id")[[mn.sample2, "Component"]]
    edges = mn.xEDGES.copy()
    edges = edges.merge(node_info, how="left", left_on="Source", right_index=True)
    edges = edges.rename(columns={sample2: "source_sample", "Component": "source_component"})
    edges = edges.merge(node_info, how="left", left_on="Target", right_index=True)
    edges = edges.rename(columns={sample2: "target_sample", "Component": "target_component"})
    mask_forward = (edges["source_sample"] == sample_from) & (edges["target_sample"] == sample_to)
    if component_source is not None:
        mask_forward &= (edges["source_component"] == component_source)
    forward_edges = edges[mask_forward]
    forward_counts = forward_edges['Transformation'].value_counts()
    mask_backward = (edges["source_sample"] == sample_to) & (edges["target_sample"] == sample_from)
    if component_source is not None:
        mask_backward &= (edges["source_component"] == component_source)
    backward_edges = edges[mask_backward]
    backward_counts = backward_edges['Transformation'].value_counts()
    all_transformations = set(forward_counts.index).union(set(backward_counts.index))
    transformation_counts = pd.DataFrame(index=all_transformations)
    transformation_counts['Forward'] = forward_counts
    transformation_counts['Backward'] = backward_counts
    clusterTransitions = transformation_counts = transformation_counts.fillna(0).astype(int)
    #
    if transformation_counts.empty:
        print("##########################\nNO TRANSITION TRAFOS FOUND\n##########################")
        return
        #
    fig, ax = plt.subplots(figsize=(6,5))
    transformation_counts.sort_index().plot(kind='bar',     stacked=True,     color=["black", "white"],     edgecolor='black',     ax=ax ) 
    ax.set_ylabel('Counts')
    ax.set_xlabel('Transformation')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, ha='center')
    ax.legend(title='Direction')
    fig.tight_layout()
    fig.subplots_adjust(top=.95, bottom=.15)  # increase top and bottom padding
    plt.title(f"{sample_from} $\\rightarrow$ {sample_to}", fontsize=12, wrap=True)
    plt.text(.01,.01, f"{name}, component {component_source}", fontsize=10, wrap=True, transform=plt.gcf().transFigure)
    plt.savefig(f"{date}  {name}  Comp{FILTER_COMPONENTS}  {trafos2}  count_sample_transformations  {sample2}.png")
    #plt.show()
    return clusterTransitions































































































'''  __x

def degDistri(xNODES):      
    fnName = '_degDistr'
    x_name = 'Degree'
    y_name = 'Degree counts'
    text   = 'id'
    x = xNODES[x_name].tolist()
    y = xNODES[y_name].tolist()
    text = xNODES[text].tolist()
    #
    plt.figure(figsize=(6,5))
    plt.title(name+'\n, Comp'+str(FILTER_COMPONENTS)+', Trafos'+str(trafos), wrap=True, fontsize=12)
    plt.scatter(x,y,c='black', s=100)
##    for i in range(len(x)):
##        plt.annotate(text[i], (x[i], y[i] + 0.2), fontsize=10)
    #
      
##    plt.xscale('log',base=10) 
##    plt.yscale('log',base=10) 
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    plt.savefig(date +'  ' +name+'  Comp'+str(FILTER_COMPONENTS)+'  '+str(fnName)+'  '+str(sample2)+'.png')
##    plt.show()

def degDistriLog(xNODES):      
    fnName = '_degDistrLog'
    x_name = 'Degree'
    y_name = 'Degree counts'
    text   = 'id'
    x = xNODES[x_name].tolist()
    y = xNODES[y_name].tolist()
    text = xNODES[text].tolist()
    #
    plt.figure(figsize=(6,5))
    plt.title(name+'\n, Comp'+str(FILTER_COMPONENTS)+', Trafos'+str(trafos), wrap=True, fontsize=12)
    plt.scatter(x,y,c='black', s=100)
##    for i in range(len(x)):
##        plt.annotate(text[i], (x[i], y[i] + 0.2), fontsize=10)
    #
    plt.xscale('log',base=10) 
    plt.yscale('log',base=10) 
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    plt.savefig(date +'  ' +name+'  Comp'+str(FILTER_COMPONENTS)+'  '+str(fnName)+'  '+str(sample2)+'.png')
##    plt.show()

def edgeCounts(xEDGES):       
    fnName = '_edgeCounts'
    a  = xEDGES.drop_duplicates(subset=['Mass difference'], keep='first')
    a2 = a[['Molecular Formula','Transformation (frequency)']]
    a2.index = a2['Molecular Formula']
    a2 = a2.drop(['Molecular Formula'],axis=1)
    #
    a2.plot.bar(rot=0, color = 'black', legend=False)
    plt.ylabel('Transformation counts')
    plt.savefig(date +'  ' +name+'  Comp'+str(FILTER_COMPONENTS)+'  '+str(fnName)+'  '+str(sample2)+'.png')
    ##    plt.show()











def  _x  ###################################################   


def annotation(G,nodeColor):                   ######################   _annotation                                                 




colorFamilies = ['Reds','Oranges','Greens','Blues','Purples','bone','pink','spring','copper','Wistia']
colorFamilies = colorFamilies[0:len(trafos)]
####################################   _splittedTrafos
z = []
for i in trafos:
    splittedTrafos = [x for x in trafos if x != i]             
    trafoWeight = [1]*len(splittedTrafos)
    trafos(splittedTrafos,trafoWeight)
    #
    fileName        =       "bennu.dat"       
    molFormula      =       'Composition'
    abundance       =       'Intensity'
    sample          =       'Sample'                       
    nComps          =       1
    node_size       =       'Degree'
    molForm2graph(fileName,molFormula,abundance,sample,nComps,node_size)
    ####################################  _std
    l = ['H', 'C', 'N', 'O', 'S',         'HC', 'NC', 'OC', 'SC', 'Mass (exact)',         'Abundance', 'DBE', 'X_C', 'AMD', 'KMD']           
    a = xNODES[l]
    b=[]
    for i in xNODES['Component'].unique().tolist():
        a = xNODES[xNODES['Component'] == i]
        a2 = pd.concat([a, a.describe().loc[["std"]]])
        a3 = a2.iloc[[-1]]
        a2['sum_std'] = a3.sum(axis=1)
        b.append({'Component':i,'sum_std':a2.iloc[-1]['sum_std']})
    #    
    a = pd.DataFrame.from_dict(b)
    b = xNODES[['df_comp_freq','Component']]
    c = a.merge(b, how='left', on='Component')
    d = c.drop_duplicates(subset=['Component'], keep='first')
    d = d.sort_values(by=['sum_std'], ascending=True)
    n = 95
    e = d[ ( d['df_comp_freq'] >= np.percentile(d['df_comp_freq'], n) ) & ( d['sum_std'] <= np.percentile(d['sum_std'], n) ) ]
    ##d['chooceComp'] = d['df_comp_freq'] - d['sum_std']
    f = e[0:5]
    #
    z.append(f['Component'].tolist())

z = np.concatenate(z)
#
colors=[]
for j in colorFamilies:
    c = mcp.gen_color(cmap=j,n=len(f))    ## https://matplotlib.org/stable/tutorials/colors/colormaps.html
    colors.append(c)

colors = np.concatenate(colors)
g3 = pd.DataFrame({'Component': z, 'colorSplitted': colors})
h = xNODES.merge(g3, how='left', on='Component')
h['colorSplitted'] = h['colorSplitted'].fillna('lightgray')
xNODES = h
xNODES = xNODES.drop_duplicates(subset=['id'], keep='first')
colorSplitted = xNODES['colorSplitted'].to_list()

sampleColorSplitted = xNODES[['Sample','colorSplitted']]











###################################   _annotation


plt.figure(figsize=(11,5))
pos = nx.nx_pydot.graphviz_layout(G)
nx.draw(G, pos = pos, with_labels=False, node_color=colorSplitted, node_size=5, width=.01, alpha=1)
##for i in range(len(sampleColorSplitted)):
##    plt.plot([], [], sampleColorSplitted.values[i][1], marker='o', markersize=10, label=sampleColorSplitted.values[i][0])      


plt.show()



ValueError: s must be a scalar, or float array-like with the same size as x and y







for comp in FILTER_COMPONENTS_filtered:
    component_data = filtered_nodes[filtered_nodes['Component'] == comp]
    node_id = component_data.iloc[0]['id']
    x, y = pos[node_id]
    stats = {'C_min': component_data['C'].min(), 'C_max': component_data['C'].max(), 'H_min': component_data['H'].min(), 'H_max': component_data['H'].max(), 'N_min': component_data['N'].min(), 'N_max': component_data['N'].max(), 'O_min': component_data['O'].min(), 'O_max': component_data['O'].max(), 'S_min': component_data['S'].min(), 'S_max': component_data['S'].max(), 'mz_range': f"{component_data['Mass (exact)'].min():.2f}-{component_data['Mass (exact)'].max():.2f}", 'DBE_range': f"{component_data['DBE'].min():.2f}-{component_data['DBE'].max():.2f}", 'num_nodes': len(component_data)}
    text = f"Component: {comp}\n\          C({stats['C_min']} - {stats['C_max']})H({stats['H_min']} - {stats['H_max']})N({stats['N_min']} - {stats['N_max']})O({stats['O_min']} - {stats['O_max']})\n\             S: {stats['S_min']} - {stats['S_max']}\n\             mz range: {stats['mz_range']}\n\             DBE range: {stats['DBE_range']}\n\             Mol. For.: {stats['num_nodes']}" 
    plt.text(x-550, y, text, fontsize=8, ha='left', va='center')


plt.show()


#############################################################COMBINATION######################################################################################

##combine those dataframes 
##load xnodes
#xNODES_annotated_NH3 = pd.read_csv('NH3.csv')
#xNODES_annotated_CO = pd.read_csv('CO.csv')
#xNODES_annotated_O = pd.read_csv('O.csv')
#xNODES_annotated_H2 = pd.read_csv('H2.csv')
#plt.figure()  #adjust the figure size as needed figsize=(20, 20)
#plt.title("Selected components with chemical properties", wrap=True, fontsize=12)

##create a list of the specified dataframes
#annotated_dataframes = [xNODES_annotated_NH3, xNODES_annotated_CO, xNODES_annotated_O, xNODES_annotated_H2]

##create a dictionary to store node colors from the specified dataframes
#node_colors_dict = {}

##iterate over each specified dataframe
#for df in annotated_dataframes:
#    # Filter nodes with non-gray colors
#    non_gray_nodes = df[df['colored_subgraph'] != 'gray']
#    # Update node colors in the dictionary
#    node_colors_dict.update(dict(zip(non_gray_nodes['id'], non_gray_nodes['colored_subgraph'])))

##define function to get node color
#def get_node_color(node_id):
#    return node_colors_dict.get(node_id, '#808080')  # Default to gray if color not found

### Apply function to each node in xNODES to get its color
#xNODES['node_color'] = xNODES['id'].apply(get_node_color)

#### Plot the graph with the modified colors
#nx.draw(G, pos=nx.nx_pydot.graphviz_layout(G), with_labels=False, node_color=xNODES['node_color'], node_size=0.5, width=.01, alpha=1)

#### Plot the sample colors
##legend_entries = [('orange', 'NH3'), ('black', 'CO'), ('green', 'O'), ('blue', 'H2')]
##for color, label in legend_entries:
##    plt.scatter([], [], color=color, marker='o', s=100, label=label)
##    
##plt.legend()
#plt.text(.6, -1.1, str(G.number_of_nodes())+' nodes, '+str(G.number_of_edges())+' edges', fontsize=8, wrap=True)
#plt.savefig('full_layout_colored'+'.png')
######plt.show()
















##                        
##                        
##                        
##                        
##                        
##                        splittedTrafos = ['H2','O']
##                        trafoWeight = [1]*len(splittedTrafos)
##                        trafos(splittedTrafos,trafoWeight)
##                        #
##                        fileName        =       "bennu.dat"       
##                        molFormula      =       'Composition'
##                        abundance       =       'Intensity'
##                        sample          =       'Sample'                       
##                        nComps          =       1
##                        node_size       =       'Degree'
##                        molForm2graph(fileName,molFormula,abundance,sample,nComps,node_size)
##                        ####################################  _std
##                        l = ['H', 'C', 'N', 'O', 'S',         'HC', 'NC', 'OC', 'SC', 'Mass (exact)',         'Abundance', 'DBE', 'X_C', 'AMD', 'KMD']           
##                        a = xNODES0[l]
##                        b=[]
##                        for i in xNODES0['Component'].unique().tolist():
##                            a = xNODES[xNODES0['Component'] == i]
##                            a2 = pd.concat([a, a.describe().loc[["std"]]])
##                            a3 = a2.iloc[[-1]]
##                            a2['sum_std'] = a3.sum(axis=1)
##                            b.append({'Component':i,'sum_std':a2.iloc[-1]['sum_std']})
##                        #    
##                        a = pd.DataFrame.from_dict(b)
##                        b = xNODES0[['df_comp_freq','Component']]
##                        c = a.merge(b, how='left', on='Component')
##                        d = c.drop_duplicates(subset=['Component'], keep='first')
##                        d = d.sort_values(by=['sum_std'], ascending=True)
##                        n = 95
##                        e = d[ ( d['df_comp_freq'] >= np.percentile(d['df_comp_freq'], n) ) & ( d['sum_std'] <= np.percentile(d['sum_std'], n) ) ]
##                        ##d['chooceComp'] = d['df_comp_freq'] - d['sum_std']
##                        f = e[0:5]
##                        #
##                        colors = mcp.gen_color(cmap='Reds',n=len(f))    ## https://matplotlib.org/stable/tutorials/colors/colormaps.html
##                        g = pd.DataFrame({'Component': f['Component'].tolist(), 'color': colors})
##                        h = xNODES0.merge(g, how='left', on='Component')
##                        
##                        
##                        
##                        h['color'] = h['color'].fillna('lightgray')
##                        xNODES = h
##                        








^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             ####################SPLITTING##############################################
^##                                                             
^##                                                             
^##                                                             
^##                                                             ## x x      
^##                                                             
^##                                                             ## This code prints the element maps for a given trafo set
^##                                                             ## the code then uses a definition that we discussed with alex to define 5 components of interest
^##                                                             ## you can then use the 5 components of interest to study certain properties of your choice
^##                                                             
^##                                                             ## Make sure to generate the network first! With your trafoset of choice
^##                                                             ## this example uses trafos          =       ['H2','O','CO'] #so minus NH3
^##                                                             ## Increase the number of components to have a greater accuracy
^##                                                             ## (it makes sense to have more than 5... since you want to select 5 in the end!! :) )
^##                                                             
^##                                                             ## PART 1 : Element maps and Standard Deviation
^##                                                             
^##                                                             
^##                                                             ## List of elements for visualization
^##                                                             #l = ["H", "C", "N", "O", "HC", "NC", "OC", "Mass (exact)", "DBE"]
^##                                                             ##l = ["H", "C"]
^##                                                             ## Create an empty DataFrame to store the sum of standard deviations for each subgraph
^##                                                             #sum_std_df = pd.DataFrame(index=xNODES["df_comp_freq"].unique(), columns=l)
^##                                                             
^##                                                             ## Iterate over elements for visualization
^##                                                             #for color_label in l:
^##                                                             #    carac = xNODES[["id", color_label]]
^##                                                             #    carac = carac.set_index("id")
^##                                                             #    carac = carac.reindex(G.nodes())
^##                                                             #    carac[color_label] = pd.Categorical(carac[color_label])
^##                                                             #    carac[color_label].cat.codes
^##                                                             #    nodes = G.nodes()
^##                                                             #    # in the code above, we just color nodes based on the value of the element in question (in l)
^##                                                             ##
^##                                                             ##    plt.figure()
^##                                                             ##    plt.title(color_label, wrap=True, fontsize=25)
^##                                                             ###
^##                                                             ##    pos = nx.nx_pydot.graphviz_layout(G)
^##                                                             ##    nx.draw(G,pos, with_labels=False,node_color=node_colors,node_size=2,width=0.01,alpha=1)
^##                                                             ##    ec = nx.draw_networkx_edges(G, pos, width=0.01, alpha=0.2)
^##                                                             ##    nc = nx.draw_networkx_nodes(G, pos,nodelist=nodes,node_color=carac[color_label],node_size=3,cmap=plt.cm.jet)
^##                                                             ###
^##                                                             ##    cb = plt.colorbar(nc, orientation="vertical")
^##                                                             ##    nc.figure.axes[0].tick_params(axis="both", labelsize=21)
^##                                                             ##    nc.figure.axes[1].tick_params(axis="y", labelsize=21)
^##                                                             ##    plt.axis("off")
^##                                                             ##    plt.savefig("Elements"+str(trafos)+color_label+".png" )
^##                                                             ###    plt.show()
^##                                                             ##
^##                                                             #    # Calculate the standard deviation of H for each subgraph
^##                                                             #    subgraph_std = xNODES.groupby("df_comp_freq")[color_label].std()
^##                                                             ##
^##                                                             #    # Normalize the standard deviations to use for coloring
^##                                                             #    normalized_std = (subgraph_std - subgraph_std.min()) / (subgraph_std.max() - subgraph_std.min())
^##                                                             ##
^##                                                             #    # Create a dictionary to map subgraphs to their normalized standard deviations
^##                                                             #    subgraph_color_dict = dict(zip(normalized_std.index, normalized_std.values))
^##                                                             ##
^##                                                             #    # Map the normalized standard deviations to each subgraph in the graph
^##                                                             #    node_colors = [subgraph_color_dict.get(df_comp_freq, 0) for df_comp_freq in xNODES["df_comp_freq"]]
^##                                                             ##
^##                                                             #    # Store the normalized standard deviations in the sum_std_df DataFrame
^##                                                             #    sum_std_df[color_label] = normalized_std
^##                                                             ##
^##                                                             ## Sum the standard deviations for each subgraph across all elements
^##                                                             #sum_std_df["sum_std"] = sum_std_df.sum(axis=1)
^##                                                             
^##                                                             ### Plot the graph with the new colormap based on the sum of standard deviations
^##                                                             ##plt.figure()
^##                                                             ##plt.title("Minimal Set - "+sub, wrap=True, fontsize=20)
^##                                                             
^##                                                             ##pos = nx.nx_pydot.graphviz_layout(G)
^##                                                             ##nx.draw(G, pos, with_labels=False, node_color=node_colors, node_size=3, width=0.01, alpha=1)
^##                                                             ##ec = nx.draw_networkx_edges(G, pos, width=0.01, alpha=0.2)
^##                                                             ##nc = nx.draw_networkx_nodes(G, pos, nodelist=G.nodes(), node_color=
^##                                                             ##node_colors, node_size=3, cmap=plt.cm.jet)
^##                                                             
^##                                                             ### Add colorbar
^##                                                             ##cb = plt.colorbar(nc, orientation="vertical")
^##                                                             ##nc.figure.axes[0].tick_params(axis="both", labelsize=21)
^##                                                             ##nc.figure.axes[1].tick_params(axis="y", labelsize=21)
^##                                                             
^##                                                             ##plt.axis("off")
^##                                                             ###plt.savefig(str(trafos)+str(FILTER_COMPONENTS)+'Sum_of_Standard_Deviations.png')
^##                                                             ##plt.show()
^##                                                             
^##                                                             #############plot component map for subset ######
^##                                                             #fontsize=17
^##                                                             #color_label = 'component'
^##                                                             #carac = xNODES[['id', color_label]]
^##                                                             #carac = carac.set_index('id')
^##                                                             #carac = carac.reindex(G.nodes())
^##                                                             #carac[color_label]=pd.Categorical(carac[color_label])
^##                                                             #carac[color_label].cat.codes
^##                                                             #nodes = G.nodes()
^##                                                             #plt.figure()
^##                                                             #plt.title("Minimal Set - "+sub, wrap=True, fontsize=17)
^##                                                             #pos = nx.nx_pydot.graphviz_layout(G)
^##                                                             #nx.draw(G, pos, with_labels=False, node_color=carac[color_label].cat.codes, cmap=plt.cm.jet, node_size=3, width=.1, font_size=10)
^##                                                             #############pos = nx.spring_layout()
^##                                                             ###########pos = nx.fruchterman_reingold_layout(G0)
^##                                                             #ec = nx.draw_networkx_edges(G, pos, alpha=0.2)
^##                                                             ######################    nc = nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=carac[color_label], with_labels=False, node_size=5, cmap=plt.cm.jet) 
^##                                                             #nc = nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=carac[color_label], node_size=3, cmap=plt.cm.jet) 
^##                                                             #####################plt.colorbar(nc)
^##                                                             
^##                                                             #cb = plt.colorbar(nc, orientation='vertical').set_label(label=color_label, size=fontsize+2)
^##                                                             #nc.figure.axes[0].tick_params(axis="both", labelsize=21)           ## change the label size
^##                                                             #nc.figure.axes[1].tick_params(axis="y", labelsize=21)              ## change the tick label size of colorbar
^##                                                             
^##                                                             #plt.axis('off')
^##                                                             #plt.savefig(str(trafos)+str(FILTER_COMPONENTS)+'dividednetwork.png')
^##                                                             ##plt.show()
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             ###############COMPONENTS WITH SMALLEST STD #######################################
^##                                                             ## Now you have a plot showing the standard deviations of each component in your splitted set
^##                                                             ## we now want to pick the 5 most interesting components
^##                                                             
^##                                                             
^##                                                             ## Size criteria
^##                                                             ## top 25-percentile code
^##                                                             #n = len(subgraph_std)
^##                                                             #n2 = int(np.around(n / 4))
^##                                                             #size_subgraph_std = subgraph_std.tail(n2)
^##                                                             
^##                                                             #xNODES["colored_subgraph"] = -1  # Default value, you can change this to any other default value
^##                                                             #xNODES["colored_subgraph_component"] = -1  # Default value, you can change this to any other default value
^##                                                             
^##                                                             ## Find the 5 subgraphs with the smallest standard deviation
^##                                                             #top5_subgraphs = size_subgraph_std.nsmallest(5).index
^##                                                             
^##                                                             ## Create a dictionary to map subgraphs to their index for the top 5
^##                                                             #top5_subgraph_dict = dict(zip(top5_subgraphs, range(5)))
^##                                                             
^##                                                             ## map the top 5 subgraphs to orange color and update 'colored_subgraph' and 'colored_subgraph_component'
^##                                                             #node_colors_top5 = ["orange" if df_comp_freq in top5_subgraphs else "gray" for df_comp_freq in xNODES["df_comp_freq"]]
^##                                                             #xNODES.loc[xNODES["df_comp_freq"].isin(top5_subgraphs), "colored_subgraph"] = "orange"
^##                                                             #xNODES.loc[xNODES["df_comp_freq"].isin(top5_subgraphs), "colored_subgraph_component"] = xNODES["component"]
^##                                                             
^##                                                             ##################### Plot the graph with the top 5 subgraphs in orange
^##                                                             #plt.figure()
^##                                                             #plt.title("Minimal Set -"+sub+"\nBest components", wrap=True, fontsize=17)
^##                                                             ### plt.title('Top 5 Subgraphs with Smallest Standard Deviation', wrap=True, fontsize=25)
^##                                                             
^##                                                             #pos = nx.nx_pydot.graphviz_layout(G)
^##                                                             #nx.draw(G,pos,with_labels=False,node_color=node_colors_top5,node_size=3,width=0.01,alpha=1,)
^##                                                             #ec = nx.draw_networkx_edges(G, pos, width=0.01, alpha=0.2)
^##                                                             #nc = nx.draw_networkx_nodes(G,pos,nodelist=G.nodes(),node_color=node_colors_top5,node_size=3,cmap=plt.cm.jet,)
^##                                                             
^##                                                             ### colorbar
^##                                                             ##cb = plt.colorbar(nc, orientation="vertical")
^##                                                             ##nc.figure.axes[0].tick_params(axis="both", labelsize=21)  # Change label size
^##                                                             ##nc.figure.axes[1].tick_params(axis="y", labelsize=21)  # Change tick label size of colorbar
^##                                                             
^##                                                             #plt.axis("off")
^##                                                             #plt.savefig(str(trafos)+str(FILTER_COMPONENTS)+'Top_5_Subgraphs_with_Smallest_Std.png')
^##                                                             ##plt.show()
^##                                                             
^##                                                             ## PART 3: Saving and where the components are Top_5_Subgraphs_with_Smallest_Std_NH3_colored
^##                                                             
^##                                                             ## To store the information for each component, you can copy the xNODES DataFrame
^##                                                             #xNODES_annotated_NH3 = xNODES.copy()  # here we store the entire DataFrame
^##                                                             
^##                                                             ## in this dataframe, we also have the information of the 5 selected _minimal_set_NH3_all_components
^##                                                             ## you can retrieve them by color xNODES.loc[xNODES['colored_subgraph'] == '#orange' ]
^##                                                             ## and yes...  I use the color here to find the components again! We know that the colored ones
^##                                                             ## are the ones we selected
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             ################## SHADE OF COLORS ##########################################
^##                                                             
^##                                                             ## What if you wanted to differentiate between the 5 components?
^##                                                             ## they have the same color until now
^##                                                             ## We can change that here
^##                                                             ## If you want to color each one of the 5 components of choice with a different colorn
^##                                                             ## you can use this code
^##                                                             
^##                                                             ## initialize 'colored_subgraph' column with 'gray'
^##                                                             #xNODES["colored_subgraph"] = "gray"
^##                                                             
^##                                                             ## find the 5 subgraphs with the smallest standard deviation
^##                                                             #top5_subgraphs = size_subgraph_std.nsmallest(5).index
^##                                                             
^##                                                             ### define shades of color for the top 5 subgraphs
^##                                                             ### Pick 5 (clearly distinct) shades to color your components
^##                                                             #shades_of_orange = ["#DFAA4A", "#FF8C00", "#FFA07A", "#FF7F50", "#FFD700"] #minimal subset - NH3
^##                                                             #shades_of_blue = ['#0000FF', '#4169E1', '#6495ED', '#4682B4', '#87CEEB'] #minimal subset - CO
^##                                                             #shades_of_green = ['#008000', '#00B300', '#00CC66', '#93C572', '#00FFA6'] #minimal subset - O
^##                                                             #shades_of_purple = ['#B400B4', '#A27BA2', '#FF99FF', '#FF33FF', '#FF99CC'] # minimal subset - H2
^##                                                             
^##                                                             ##
^##                                                             
^##                                                             ## Assign each component a different shade of orange
^##                                                             #for i, comp in enumerate(top5_subgraphs):
^##                                                             #    # Calculate the color for this component
^##                                                             ##    color = shades_of_orange[i]
^##                                                             ##    color = shades_of_blue[i]
^##                                                             ##    color = shades_of_green[i]
^##                                                             #    color = shades_of_purple[i]
^##                                                             ##
^##                                                             #    # Update 'colored_subgraph' column for nodes in this component to the calculated color
^##                                                             #    xNODES.loc[xNODES["df_comp_freq"] == comp, "colored_subgraph"] = color
^##                                                             
^##                                                             ## Plot the graph with the top 5 subgraphs in different shades of color
^##                                                             #plt.figure()
^##                                                             #plt.title("Minimal Set -"+sub+"\nBest Components", wrap=True, fontsize=17)
^##                                                             #pos = nx.nx_pydot.graphviz_layout(G)
^##                                                             #nx.draw(G,pos,with_labels=False,node_color=xNODES["colored_subgraph"],node_size=3,width=0.01,alpha=1,)
^##                                                             #plt.axis("off")
^##                                                             #plt.savefig('Top_5_Subgraphs_with_Smallest_Std_colored'+str(FILTER_COMPONENTS)+str(trafos)+'.png')
^##                                                             ##plt.show()
^##                                                             
^##                                                             ## Example of finding the nodes that we selected
^##                                                             ## the nodes that form the same component have the same color
^##                                                             
^##                                                             #xNODES.loc[xNODES["colored_subgraph"] == "#FFA500"]
^##                                                             
^##                                                             ## make sure to save the dataframe!!! if you generate another network, it'll be stored in xNODES
^##                                                             ## and that will overwrite our work!!!!!!!
^##                                                             #xNODES_annotated_NH3 = xNODES.copy()  # here we store the entire DataFrame
^##                                                             #xNODES_annotated_NH3.to_csv(sub+'.csv', index=False)
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             ####### ANNOTATION ##
^##                                                             ##get unique components with color 'gray'
^##                                                             #gray_components = set(xNODES[xNODES['colored_subgraph'] == 'gray']['component'].unique())
^##                                                             ##filter out gray components from FILTER_COMPONENTS
^##                                                             #FILTER_COMPONENTS_filtered = np.array([comp for comp in FILTER_COMPONENTS if comp not in gray_components])
^##                                                             ##now FILTER_COMPONENTS_filtered contains only the IDs of components that do not have the color 'gray' under 'colored_subgraph'
^##                                                             #plt.figure()
^##                                                             #plt.title("Minimal Set - "+sub, wrap=True, fontsize=17)
^##                                                             ##filter nodes based on FILTER_COMPONENTS_filtered
^##                                                             #filtered_nodes = xNODES[xNODES['component'].isin(FILTER_COMPONENTS_filtered)]
^##                                                             ##filter colors for the filtered nodes
^##                                                             #filtered_colors = []
^##                                                             #for _, node_data in filtered_nodes.iterrows():
^##                                                             #    node_id = node_data['id']
^##                                                             #    color = node_data['colored_subgraph']
^##                                                             #    if color != 'gray':
^##                                                             #        filtered_colors.append(color)
^##                                                             #    else:
^##                                                             #        filtered_colors.append('#808080')  # Use #808080 for gray nodes
^##                                                             ##draw the graph with filtered nodes and colors
^##                                                             #pos = nx.nx_pydot.graphviz_layout(G)
^##                                                             #nx.draw(G, pos, nodelist=filtered_nodes['id'], with_labels=False, node_color=filtered_colors, node_size=1, width=.05, alpha=1)
^##                                                             ##display statistics next to each component
^##                                                             #for comp in FILTER_COMPONENTS_filtered:
^##                                                             #    #filter data for the current component
^##                                                             #    component_data = filtered_nodes[filtered_nodes['component'] == comp]
^##                                                             #    #find a node in the component and get its position
^##                                                             #    node_id = component_data.iloc[0]['id']
^##                                                             #    x, y = pos[node_id]
^##                                                             #    #stats for the component
^##                                                             #    stats = {
^##                                                             #        'C_min': component_data['C'].min(),
^##                                                             #        'C_max': component_data['C'].max(),
^##                                                             #        'H_min': component_data['H'].min(),
^##                                                             #        'H_max': component_data['H'].max(),
^##                                                             #        'N_min': component_data['N'].min(),
^##                                                             #        'N_max': component_data['N'].max(),
^##                                                             #        'O_min': component_data['O'].min(),
^##                                                             #        'O_max': component_data['O'].max(),
^##                                                             #        'S_min': component_data['S'].min(),
^##                                                             #        'S_max': component_data['S'].max(),
^##                                                             #        'mz_range': f"{component_data['Mass (exact)'].min():.2f}-{component_data['Mass (exact)'].max():.2f}",
^##                                                             #        'DBE_range': f"{component_data['DBE'].min():.2f}-{component_data['DBE'].max():.2f}",
^##                                                             #        'num_nodes': len(component_data),
^##                                                             #    }
^##                                                             #    text = f"Component: {comp}\n\
^##                                                             #            C({stats['C_min']} - {stats['C_max']})H({stats['H_min']} - {stats['H_max']})N({stats['N_min']} - {stats['N_max']})O({stats['O_min']} - {stats['O_max']})\n\
^##                                                             #            S: {stats['S_min']} - {stats['S_max']}\n\
^##                                                             #            mz range: {stats['mz_range']}\n\
^##                                                             #            DBE range: {stats['DBE_range']}\n\
^##                                                             #            Mol. For.: {stats['num_nodes']}"
^##                                                             ##
^##                                                             ##    plt.text(x-50, y, text, fontsize=8, ha='left', va='center')
^##                                                             #    plt.text(x-550, y, text, fontsize=8, ha='left', va='center')
^##                                                             ##plot setting
^##                                                             #for i in range(len(sample_colors)):
^##                                                             #    plt.plot([], [], sample_colors.values[i][1], marker='o', markersize=10, label=sample_colors.values[i][0])
^##                                                             ##
^##                                                             ##plt.legend()
^##                                                             #plt.savefig('Smallest_Std_annotated'+str(FILTER_COMPONENTS)+str(trafos)+'.png') 
^##                                                             ####plt.show()
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             #############################################################COMBINATION######################################################################################
^##                                                             
^##                                                             ##combine those dataframes 
^##                                                             ##load xnodes
^##                                                             #xNODES_annotated_NH3 = pd.read_csv('NH3.csv')
^##                                                             #xNODES_annotated_CO = pd.read_csv('CO.csv')
^##                                                             #xNODES_annotated_O = pd.read_csv('O.csv')
^##                                                             #xNODES_annotated_H2 = pd.read_csv('H2.csv')
^##                                                             #plt.figure()  #adjust the figure size as needed figsize=(20, 20)
^##                                                             #plt.title("Selected components with chemical properties", wrap=True, fontsize=12)
^##                                                             
^##                                                             ##create a list of the specified dataframes
^##                                                             #annotated_dataframes = [xNODES_annotated_NH3, xNODES_annotated_CO, xNODES_annotated_O, xNODES_annotated_H2]
^##                                                             
^##                                                             ##create a dictionary to store node colors from the specified dataframes
^##                                                             #node_colors_dict = {}
^##                                                             
^##                                                             ##iterate over each specified dataframe
^##                                                             #for df in annotated_dataframes:
^##                                                             #    # Filter nodes with non-gray colors
^##                                                             #    non_gray_nodes = df[df['colored_subgraph'] != 'gray']
^##                                                             #    # Update node colors in the dictionary
^##                                                             #    node_colors_dict.update(dict(zip(non_gray_nodes['id'], non_gray_nodes['colored_subgraph'])))
^##                                                             
^##                                                             ##define function to get node color
^##                                                             #def get_node_color(node_id):
^##                                                             #    return node_colors_dict.get(node_id, '#808080')  # Default to gray if color not found
^##                                                             
^##                                                             ### Apply function to each node in xNODES to get its color
^##                                                             #xNODES['node_color'] = xNODES['id'].apply(get_node_color)
^##                                                             
^##                                                             #### Plot the graph with the modified colors
^##                                                             #nx.draw(G, pos=nx.nx_pydot.graphviz_layout(G), with_labels=False, node_color=xNODES['node_color'], node_size=0.5, width=.01, alpha=1)
^##                                                             
^##                                                             #### Plot the sample colors
^##                                                             ##legend_entries = [('orange', 'NH3'), ('black', 'CO'), ('green', 'O'), ('blue', 'H2')]
^##                                                             ##for color, label in legend_entries:
^##                                                             ##    plt.scatter([], [], color=color, marker='o', s=100, label=label)
^##                                                             ##    
^##                                                             ##plt.legend()
^##                                                             #plt.text(.6, -1.1, str(G.number_of_nodes())+' nodes, '+str(G.number_of_edges())+' edges', fontsize=8, wrap=True)
^##                                                             #plt.savefig('full_layout_colored'+'.png')
^##                                                             ######plt.show()
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             
^##                                                             














def splittingLines(G,node_colors):                                                                    



fnName = '_splittingLines'
plt.figure(figsize=(11,5))
pos = nx.nx_pydot.graphviz_layout(G)
nx.draw(G, pos = pos, with_labels=False, node_color=node_colors, node_size=2, width=.01, alpha=1)
for i in range(len(sample_colors)):
    plt.plot([], [], sample_colors.values[i][1], marker='o', markersize=10, label=sample_colors.values[i][0])      

x_coords, y_coords = zip(*pos.values())
max_x, max_y = max(x_coords), max(y_coords)

plt.title(name, wrap=True, fontsize=12)
plt.legend()                     
plt.text(.6, -1.1,str(G.number_of_nodes())+' nodes, '+str(G.number_of_edges())+' edges', fontsize=8, wrap=True)


#add horizontal X axis with linearly spaced ticks
num_ticks_x = 5
x_ticks = np.linspace(0, max_x, num_ticks_x)
plt.xticks(x_ticks, [f'{tick:.2f}' for tick in x_ticks])
for x, label in zip(x_ticks, [f'{tick:.2f}' for tick in x_ticks]):
    plt.text(x, -0.02 * max_y, label, ha='center', va='bottom', color='black', fontsize=8)

plt.axhline(0, color='black', linestyle='--', linewidth=1, xmax=max_x)

#add vertical Y axis with linearly spaced ticks
num_ticks_y = 5
y_ticks = np.linspace(0, max_y, num_ticks_y)
plt.yticks(y_ticks, [f'{tick:.2f}' for tick in y_ticks])
for y, label in zip(y_ticks, [f'{tick:.2f}' for tick in y_ticks]):
    plt.text(-0.02 * max_x, y, label, ha='right', va='center', color='black', fontsize=8)

plt.axvline(0, color='black', linestyle='--', linewidth=1, ymax=max_y)



plt.show()








#Draw a line through points M1 and M2

component [1.0]
#M1_coords = (3371.0, 2919.0)
#M2_coords = (4046.0, 5109.0)

#component [0.0]
M1_coords = (801.0, 2580.0)
M2_coords = (1000.0, 2640.0)

#component [26.0]
#first lower line
M1_coords = (2606.0, 2340.0)
M2_coords = (3909.0, 2850.0)

#second upper line
M3_coords = (1166.0, 1681.0)
M4_coords = (3498.0, 3700.0)
#M3_coords = (2700, 3025.24)


###Find the closest nodes in pos to M1, M2, and M3
###once we draw a line, we want to find the nodes that the line go through, this is done by
###carefully picking points that are on top of nodes, and then letting the code find those nodes again using pos

#M1_node = min(pos, key=lambda node: ((pos[node][0] - M1_coords[0])**2 + (pos[node][1] - M1_coords[1])**2)**0.5)
#M2_node = min(pos, key=lambda node: ((pos[node][0] - M2_coords[0])**2 + (pos[node][1] - M2_coords[1])**2)**0.5)
M3_node = min(pos, key=lambda node: ((pos[node][0] - M3_coords[0])**2 + (pos[node][1] - M3_coords[1])**2)**0.5)
M4_node = min(pos, key=lambda node: ((pos[node][0] - M4_coords[0])**2 + (pos[node][1] - M4_coords[1])**2)**0.5)
















plt.show()











plt.savefig(date +'  ' +name+'  Comp'+str(FILTER_COMPONENTS)+'  '+str(fnName)+'.png')



####################Extend lines to the end of the plot
###################now our lines have to be extended so that it covers all the plot, then we can define our regions
##slope1 = (pos[M1_node][1] - M2_coords[1]) / (pos[M1_node][0] - M2_coords[0])
##intercept1 = pos[M1_node][1] - slope1 * pos[M1_node][0]
##x_line1 = np.linspace(min(x_coords), max(x_coords), 100)
##y_line1 = slope1 * x_line1 + intercept1

###################Extend lines to the end of the plot 
##################now our lines have to be extended so that it covers all the plot, then we can define our regions
#slope1 = (pos[M1_node][1] - M2_coords[1]) / (pos[M1_node][0] - M2_coords[0])
#intercept1 = pos[M1_node][1] - slope1 * pos[M1_node][0]
##x_line1 = np.linspace(min(x_coords), max(x_coords), 100)
#x_line1 = np.arange(3000.0, 4300.0, 1)
#y_line1 = slope1 * x_line1 + intercept1

##slope2 = (pos[M3_node][1] - M4_coords[1]) / (pos[M3_node][0] - M4_coords[0])
##intercept2 = pos[M3_node][1] - slope2 * pos[M3_node][0]
##x_line2 = np.linspace(min(x_coords), max(x_coords), 100)
##y_line2 = slope2 * x_line2 + intercept2


########Assign 'Sample' values based on regions defined by the lines
########the regions are defined as follows:

##xNODES['sam']=xNODES['sample']
##for node, (x, y) in pos.items():
##	if y > slope2 * x + intercept2:
##		xNODES.loc[xNODES['id'] == node, 'sam'] = 'green'
###
##		plt.plot(pos[node][0], pos[node][1], marker='o', markersize=1, color='green')
###	elif y <= slope2 * x + intercept2 and y > slope1 * x + intercept1:
###		xNODES.loc[xNODES['id'] == node, 'sam'] = 'red'
###		plt.plot(pos[node][0], pos[node][1], marker='o', markersize=1, color='red')
##	elif y <= slope1 * x + intercept1:
##		xNODES.loc[xNODES['id'] == node, 'sam'] = 'blue'
##		plt.plot(pos[node][0], pos[node][1], marker='o', markersize=1, color='blue')

########Assign 'Sample' values based on regions defined by the lines
########the regions are defined as follows:

#xNODES['sam']=xNODES['sample']
#for node, (x, y) in pos.items():
#	if y > slope1 * x + intercept1:
#		xNODES.loc[xNODES['id'] == node, 'sam'] = 'red'
#		plt.plot(pos[node][0], pos[node][1], marker='o', markersize=1, color='red')
#	elif y <= slope1 * x + intercept1:
#		xNODES.loc[xNODES['id'] == node, 'sam'] = 'nonred'
#		plt.plot(pos[node][0], pos[node][1], marker='o', markersize=1, color='gray')

###xNODES['sam']=xNODES['sam'].astype(int)




##########plot your legend and additional text as before
###for i in range(len(sample_colors)):	
###    plt.plot([], [], sample_colors.values[i][1], marker='o', markersize=10, label=sample_colors.values[i][0])


###plt.show()
###plt.savefig(name+'LAYOUT'+str(FILTER_COMPONENTS)+str(trafos)+'.png')


###nx.draw(G, pos, with_labels=False, node_color=node_colors, node_size=2, width=.01, alpha=1)
##for i in range(len(sample_colors)):
##    plt.plot([], [], sample_colors.values[i][1], marker='o', markersize=10, label=sample_colors.values[i][0])      


#########################draw extended lines
#plt.plot(x_line1, y_line1, color='black', linestyle='--', linewidth=2)
##plt.plot(x_line2, y_line2, color='black', linestyle='--', linewidth=2)
#############################plt.plot(x_line2, y_line2, color='green', linestyle='--', linewidth=2)
###############plt.show()

####plt.legend()
####plt.text(.6, -1.1, str(G.number_of_nodes())+' nodes, '+str(G.number_of_edges())+' edges', fontsize=8, wrap=True)
#######plt.savefig(name+'LAYOUT'+str(FILTER_COMPONENTS)+str(trafos)+'.png')

####plt.savefig(name+'LAYOUT_probe'+str(FILTER_COMPONENTS)+str(trafos)+'.png')

#plt.savefig(name+'LAYOUT_colors'+str(FILTER_COMPONENTS)+str(trafos)+'.png')





































import mol2net as mn, time, networkx as nx;import time;start=time.time(); print(f"\n\n\n\n{'#'*60}\n{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))}\n{'#'*60}\n\n\n\n");

mn.molForm2graph(fileName='iceAurelien_mol2netAnalChem2022.dat' , molFormula='molecular_formula',abundance='abundance_int',sample="Ice composition H2O_CH3OH_NH3", comps=list(range(0,10)), trafos=['H2','O','CO','NH3']);
mn.molForm2graph(fileName='iceAurelien_mol2netAnalChem2022.dat' , molFormula='molecular_formula',abundance='abundance_int', comps=list(range(0,10)), trafos=['H2','O','CO','NH3']);
mn.molForm2graph(fileName='iceAurelien_mol2netAnalChem2022.dat' , molFormula='molecular_formula',abundance='abundance_int',sample="Ice composition H2O_CH3OH_NH3", comps=[0,1], trafos=['H2','O','CO','NH3']);
mn.molForm2graph(fileName='iceAurelien_mol2netAnalChem2022.dat' , molFormula='molecular_formula',abundance='abundance_int',comps=[0,1], trafos=['H2','O','CO','NH3']);

mn.compDistri2();
mn.drawNet(nodeColors=[['3_1_1','3_1_0.2','3_1_5','10_1_1','3_1_1 overirradiation','3_1_1 16h'],['green','blue','red','gray','orange','darkgreen']]);
mn.drawNet();
mn.elementMaps(elements= ['H','C','N','O' , 'H/C','Mass (exact)','N/C','O/C' ,'DBE','X_C','AMD','KMD' ],nCols=5); 
mn.elementMaps(elements= ['H','C','N','O' , 'H/C','Mass (exact)','N/C','O/C' ,'DBE','X_C','AMD','KMD' ]); 
mn.TrafoCounts(FROM_TO="Ice composition H2O_CH3OH_NH3",SAMPLE_START="3_1_1",SAMPLE_END="3_1_5");













make package  git
time line of code    progressbar  tqdm
add annotation
add splitting lines    trafoHistos
add distanceMatrix_formAnnot
dbe m fraction as a input variable
do graph metric map like element maps
compDistri2   adjust color input
{start} instead of {date} in file name
make sample and comps opt input para
X_C (m =)
correlate element mpas
trafoCounts  add comp number


chemicalFamily  CHN  CHNO  CHO
Component                     
0                54  2209   94
1               169   583    0
2                 0     4    0
3                 2     1    0
4                 0     3    0
5                 0     3    0
6                 3     0    0
7                 0     3    0
8                 0     0    2
9                 0     2    0
Components with maximum sample diversity: [0]       refine this..


 
















_xx  _opt (optional)  _nodes _edges __x (notes)                                                                      '''
