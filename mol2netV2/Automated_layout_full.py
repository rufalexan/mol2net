# This code prints the element maps for a given trafo set
# the code then uses a definition that we discussed with alex to define 5 components of interest
# you can then use the 5 components of interest to study certain properties of your choice

# Make sure to generate the network first! With your trafoset of choice
# this example uses trafos          =       ['H2','O','CO'] #so minus NH3
# Increase the number of components to have a greater accuracy
# (it makes sense to have more than 5... since you want to select 5 in the end!! :) )

# PART 1 : Element maps and Standard Deviation


# List of elements for visualization
l = ["H", "C", "N", "O", "HC", "NC", "OC", "Mass (exact)", "DBE"]

# Create an empty DataFrame to store the sum of standard deviations for each subgraph
sum_std_df = pd.DataFrame(index=xNODES["df_comp_freq"].unique(), columns=l)

# Iterate over elements for visualization
for color_type in l:
    carac = xNODES[["id", color_type]]
    carac = carac.set_index("id")
    carac = carac.reindex(G.nodes())
    carac[color_type] = pd.Categorical(carac[color_type])
    carac[color_type].cat.codes
    nodes = G.nodes()
    # in the code above, we just color nodes based on the value of the element in question (in l)

    plt.figure()
    plt.title(color_type, wrap=True, fontsize=25)

    pos = nx.nx_pydot.graphviz_layout(G)
    nx.draw(
        G,
        pos,
        with_labels=False,
        node_color=node_colors,
        node_size=2,
        width=0.01,
        alpha=1,
    )
    ec = nx.draw_networkx_edges(G, pos, width=0.01, alpha=0.2)
    nc = nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=nodes,
        node_color=carac[color_type],
        node_size=3,
        cmap=plt.cm.jet,
    )

    cb = plt.colorbar(nc, orientation="vertical")
    nc.figure.axes[0].tick_params(axis="both", labelsize=21)
    nc.figure.axes[1].tick_params(axis="y", labelsize=21)
    plt.axis("off")
    plt.savefig(
        "Supporting Figure 7 - "
        + color_type
        + "_minimal_set_NH3_all_components"
        + ".png"
    )
    plt.show()

    # Calculate the standard deviation of H for each subgraph
    subgraph_std = xNODES.groupby("df_comp_freq")[color_type].std()

    # Normalize the standard deviations to use for coloring
    normalized_std = (subgraph_std - subgraph_std.min()) / (
        subgraph_std.max() - subgraph_std.min()
    )

    # Create a dictionary to map subgraphs to their normalized standard deviations
    subgraph_color_dict = dict(zip(normalized_std.index, normalized_std.values))

    # Map the normalized standard deviations to each subgraph in the graph
    node_colors = [
        subgraph_color_dict.get(df_comp_freq, 0)
        for df_comp_freq in xNODES["df_comp_freq"]
    ]

    # Store the normalized standard deviations in the sum_std_df DataFrame
    sum_std_df[color_type] = normalized_std

# Sum the standard deviations for each subgraph across all elements
sum_std_df["sum_std"] = sum_std_df.sum(axis=1)

# Plot the graph with the new colormap based on the sum of standard deviations
plt.figure()
plt.title("Sum of Standard Deviations in Each Subgraph", wrap=True, fontsize=25)

pos = nx.nx_pydot.graphviz_layout(G)
nx.draw(
    G, pos, with_labels=False, node_color=node_colors, node_size=3, width=0.01, alpha=1
)
ec = nx.draw_networkx_edges(G, pos, width=0.01, alpha=0.2)
nc = nx.draw_networkx_nodes(
    G, pos, nodelist=G.nodes(), node_color=node_colors, node_size=3, cmap=plt.cm.jet
)

# Add colorbar
cb = plt.colorbar(nc, orientation="vertical")
nc.figure.axes[0].tick_params(axis="both", labelsize=21)
nc.figure.axes[1].tick_params(axis="y", labelsize=21)

plt.axis("off")
plt.savefig("Sum_of_Standard_Deviations.png")
plt.show()

# PART 2: Selecting components of interest
# Now you have a plot showing the standard deviations of each component in your splitted set
# we now want to pick the 5 most interesting components


# Size criteria
# top 25-percentile code
n = len(subgraph_std)
n2 = int(np.around(n / 4))
size_subgraph_std = subgraph_std.tail(n2)

xNODES[
    "colored_subgraph"
] = -1  # Default value, you can change this to any other default value
xNODES[
    "colored_subgraph_component"
] = -1  # Default value, you can change this to any other default value

# Find the 5 subgraphs with the smallest standard deviation
top5_subgraphs = size_subgraph_std.nsmallest(5).index

# Create a dictionary to map subgraphs to their index for the top 5
top5_subgraph_dict = dict(zip(top5_subgraphs, range(5)))

# map the top 5 subgraphs to orange color and update 'colored_subgraph' and 'colored_subgraph_component'
node_colors_top5 = [
    "orange" if df_comp_freq in top5_subgraphs else "gray"
    for df_comp_freq in xNODES["df_comp_freq"]
]
xNODES.loc[xNODES["df_comp_freq"].isin(top5_subgraphs), "colored_subgraph"] = "orange"
xNODES.loc[
    xNODES["df_comp_freq"].isin(top5_subgraphs), "colored_subgraph_component"
] = xNODES["component"]

# Plot the graph with the top 5 subgraphs in orange
plt.figure()
# plt.title('Top 5 Subgraphs with Smallest Standard Deviation', wrap=True, fontsize=25)

pos = nx.nx_pydot.graphviz_layout(G)
nx.draw(
    G,
    pos,
    with_labels=False,
    node_color=node_colors_top5,
    node_size=3,
    width=0.01,
    alpha=1,
)
ec = nx.draw_networkx_edges(G, pos, width=0.01, alpha=0.2)
nc = nx.draw_networkx_nodes(
    G,
    pos,
    nodelist=G.nodes(),
    node_color=node_colors_top5,
    node_size=3,
    cmap=plt.cm.jet,
)

# colorbar
cb = plt.colorbar(nc, orientation="vertical")
nc.figure.axes[0].tick_params(axis="both", labelsize=21)  # Change label size
nc.figure.axes[1].tick_params(
    axis="y", labelsize=21
)  # Change tick label size of colorbar

plt.axis("off")
plt.savefig("Top_5_Subgraphs_with_Smallest_Std_NH3.png")
plt.show()

# PART 3: Saving and where the components are Top_5_Subgraphs_with_Smallest_Std_NH3_colored

# To store the information for each component, you can copy the xNODES DataFrame
xNODES_annotated_NH3 = xNODES.copy()  # here we store the entire DataFrame
# in this dataframe, we also have the information of the 5 selected _minimal_set_NH3_all_components
# you can retrieve them by color xNODES.loc[xNODES['colored_subgraph'] == '#orange' ]
# and yes...  I use the color here to find the components again! We know that the colored ones
# are the ones we selected


# PART 4: Different colors for each of the 5 components

# What if you wanted to differentiate between the 5 components?
# they have the same color until now
# We can change that here
# If you want to color each one of the 5 components of choice with a different colorn
# you can use this code

# initialize 'colored_subgraph' column with 'gray'
xNODES["colored_subgraph"] = "gray"

# find the 5 subgraphs with the smallest standard deviation
top5_subgraphs = size_subgraph_std.nsmallest(5).index

# define shades of orange for the top 5 subgraphs
# Pick 5 (clearly distinct) shades to color your components
shades_of_orange = ["#DFAA4A", "#FF8C00", "#FFA07A", "#FF7F50", "#FFD700"]

#

# Assign each component a different shade of orange
for i, comp in enumerate(top5_subgraphs):
    # Calculate the color for this component
    color = shades_of_orange[i]

    # Update 'colored_subgraph' column for nodes in this component to the calculated color
    xNODES.loc[xNODES["df_comp_freq"] == comp, "colored_subgraph"] = color

# Plot the graph with the top 5 subgraphs in different shades of orange
plt.figure()
pos = nx.nx_pydot.graphviz_layout(G)
nx.draw(
    G,
    pos,
    with_labels=False,
    node_color=xNODES["colored_subgraph"],
    node_size=3,
    width=0.01,
    alpha=1,
)
plt.axis("off")
plt.savefig("Top_5_Subgraphs_with_Smallest_Std_NH3_colored.png")
plt.show()

# Example of finding the nodes that we selected
# the nodes that form the same component have the same color

xNODES.loc[xNODES["colored_subgraph"] == "#FFA500"]

# make sure to save the dataframe!!! if you generate another network, it'll be stored in xNODES
# and that will overwrite our work!!!!!!!
xNODES_annotated_NH3 = xNODES.copy()  # here we store the entire DataFrame