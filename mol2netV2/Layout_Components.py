# Layout with component representation
#this plots your network and colors the components with different colors

#Normalized color gradient 

#define the color map
cmap = cm.get_cmap('viridis')  #feel free to choose other colormaps from matplotlib! :D 

#normalize component values 
norm = plt.Normalize(xNODES['component'].min(), xNODES['component'].max())
#create a scalar mappable for color mapping
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

#list of colors for each node based on the 'component' column
node_colors_xNODES = [cmap(norm(component)) for component in xNODES['component']]

#plot the graph with node colors
plt.figure()
plt.title("Network representation", wrap=True, fontsize=12)
nx.draw(G, pos=nx.nx_pydot.graphviz_layout(G), with_labels=False, node_color=node_colors_xNODES, node_size=2, width=.01, alpha=1)

#plot the colorbar
plt.colorbar(sm, label='Component')

plt.text(.6, -1.1, f"{G.number_of_nodes()} nodes, {G.number_of_edges()} edges", fontsize=8, wrap=True)

plt.savefig('NetworkRepresentation'+'.png')

plt.show()

