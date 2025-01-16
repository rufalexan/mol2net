#in this code, we extract the coordinate positions of the nodes on the screen
#then we construct a way of selecting various clusters based on their posiiton in this coordinate system
#this is the current code used for the Paris dataset, it serves as an example and has yet to be generalized for other cases
#ask me (Zied) for further details. Until now, it is only valid for that dataset


plt.figure(figsize=(100, 100))
plt.title(name, wrap=True, fontsize=12)

nx.draw(G, pos, with_labels=False, node_color=node_colors, node_size=0.5, width=.01, alpha=1)

#extract x and y coordinates from pos
x_coords, y_coords = zip(*pos.values())

#find the maximum x and y coordinates
max_x, max_y = max(x_coords), max(y_coords)

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

#now we get to the fun part! we start by defining the lines we want to "cut" our graph with

#define points for lines
M1_coords = (5800, 6050.4)
M2_coords = (3900.60, 0)
M3_coords = (2700, 3025.24)
M4_x = 3400  #x-coordinate for the third vertical line
M5_x = 4750  #x-coordinate for the fourth vertical line

#these values were obtained by trial and error. You will have to play around with your network to find out which lines are the best to separate your graph into various clusters of interest

#find the closest nodes in pos to M1, M2, and M3
M1_node = min(pos, key=lambda node: ((pos[node][0] - M1_coords[0])**2 + (pos[node][1] - M1_coords[1])**2)**0.5)
M2_node = min(pos, key=lambda node: ((pos[node][0] - M2_coords[0])**2 + (pos[node][1] - M2_coords[1])**2)**0.5)
M3_node = min(pos, key=lambda node: ((pos[node][0] - M3_coords[0])**2 + (pos[node][1] - M3_coords[1])**2)**0.5)
#the idea here is to locate the nodes that are closest to those lines


#extend lines to the end of the plot
#slope1 = (pos[M1_node][1] - M2_coords[1]) / (pos[M1_node][0] - M2_coords[0]) 
slope1 = ((pos[M1_node][1] - M2_coords[1]) / (pos[M1_node][0] - M2_coords[0]))*1.1
 
intercept1 = pos[M1_node][1] - slope1 * pos[M1_node][0]
x_line1 = np.linspace(min(x_coords), max(x_coords), 100)
y_line1 = slope1 * x_line1 + intercept1

#slope2 = (pos[M3_node][1] - M2_coords[1]) / (pos[M3_node][0] - M2_coords[0])
slope2 = ((pos[M3_node][1] - M2_coords[1]) / (pos[M3_node][0] - M2_coords[0]))*0.7
intercept2 = pos[M3_node][1] - slope2 * pos[M3_node][0]
x_line2 = np.linspace(min(x_coords), max(x_coords), 100)
y_line2 = slope2 * x_line2 + intercept2


#now draw the extended lines! 
#draw extended lines
plt.plot(x_line1, y_line1, color='blue', linestyle='--', linewidth=2)
plt.plot(x_line2, y_line2, color='green', linestyle='--', linewidth=2)
plt.axvline(M4_x, color='purple', linestyle='--', linewidth=2)
plt.axvline(M5_x, color='red', linestyle='--', linewidth=2)


#Assign 'Sample' values based on regions defined by the lines
#for node, (x, y) in pos.items():
#    if y > slope2 * x + intercept2 and y > slope1 * x + intercept1:
#        # Over line 2 and over line 1 -> Sample 2
#        xNODES.loc[xNODES['id'] == node, 'Sample'] = '2'
#        plt.plot(pos[node][0], pos[node][1], marker='o', markersize=5, color='yellow')
#    elif y <= slope2 * x + intercept2 and x <= M4_x:
#        # Below line 2 and left of line 3 -> Sample 1
#        xNODES.loc[xNODES['id'] == node, 'Sample'] = '1'
#        plt.plot(pos[node][0], pos[node][1], marker='o', markersize=5, color='blue')
#    elif y <= slope2 * x + intercept2 and x > M4_x:
#        # Below line 2 and right of line 3 -> Sample 2
#        xNODES.loc[xNODES['id'] == node, 'Sample'] = '2'
#        plt.plot(pos[node][0], pos[node][1], marker='o', markersize=5, color='yellow')
#    elif y <= slope1 * x + intercept1 and x < M5_x:
#        # Below line 1 and left of 4 -> Sample 2
#        xNODES.loc[xNODES['id'] == node, 'Sample'] = '2'
#        plt.plot(pos[node][0], pos[node][1], marker='o', markersize=5, color='yellow')
#    elif y <= slope1 * x + intercept1 and x < M5_x:
#        # Below line 1 and right of 4 -> Sample 3
#        xNODES.loc[xNODES['id'] == node, 'Sample'] = '3'
#        plt.plot(pos[node][0], pos[node][1], marker='o', markersize=5, color='red')

#list of specific node IDs to assign to Sample 2
specific_sample_2_ids = [18118,18368,17663,18203,17254,17204,17640,17580,17894,18127,17999,17875,18013,18140,17916]  #replace with actual IDs

for node, (x, y) in pos.items():
    if node in specific_sample_2_ids:
        #specific nodes assigned to Sample 2
        xNODES.loc[xNODES['id'] == node, 'Sample'] = '2'
        plt.plot(pos[node][0], pos[node][1], marker='o', markersize=5, color='yellow')
    elif y > slope2 * x + intercept2 and y > slope1 * x + intercept1:
        #over line 2 and over line 1 -> Sample 2
        xNODES.loc[xNODES['id'] == node, 'Sample'] = '2'
        plt.plot(pos[node][0], pos[node][1], marker='o', markersize=5, color='yellow')
    elif y <= slope2 * x + intercept2 and x <= M4_x:
        #below line 2 and left of line 3 -> Sample 1
        xNODES.loc[xNODES['id'] == node, 'Sample'] = '1'
        plt.plot(pos[node][0], pos[node][1], marker='o', markersize=5, color='blue')
    elif y <= slope2 * x + intercept2 and x > M4_x:
        #below line 2 and right of line 3 -> Sample 2
        xNODES.loc[xNODES['id'] == node, 'Sample'] = '2'
        plt.plot(pos[node][0], pos[node][1], marker='o', markersize=5, color='yellow')
    elif y <= slope1 * x + intercept1 and x < M5_x:
        #below line 1 and left of 4 -> Sample 2
        xNODES.loc[xNODES['id'] == node, 'Sample'] = '2'
        plt.plot(pos[node][0], pos[node][1], marker='o', markersize=5, color='yellow')
    elif y <= slope1 * x + intercept1 and x >= M5_x:
        #below line 1 and right of 4 -> Sample 3
        xNODES.loc[xNODES['id'] == node, 'Sample'] = '3'
        plt.plot(pos[node][0], pos[node][1], marker='o', markersize=5, color='red')


        

xNODES['Sample'] = xNODES['Sample'].astype(int)

#plot legend and additional text
for i in range(len(sample_colors)):
    plt.plot([], [], sample_colors.values[i][1], marker='o', markersize=10, label=sample_colors.values[i][0])

plt.legend()
plt.text(.6, -1.1, str(G.number_of_nodes()) + ' nodes, ' + str(G.number_of_edges()) + ' edges', fontsize=8, wrap=True)
plt.show()
