#########################
#     CLUSTER TRANSITIONS                              fig 2: min_set     directions between nh3_poor to medium, etc..
#         WITH DIRECTIONS
#########################
#~~~~~~~~~~~~~~~~~~~~~~~
FROM_TO = 'Sample'
#~~~~~~~~~~~~~~~~~~~~~~~
df_samples = xNODES[['id', FROM_TO]]
df_edges4 = xEDGES.rename(columns={'Source':'id'})
df_edges4 = df_edges4.merge(df_samples, how='left', on='id')
df_edges4 = df_edges4.rename(columns={FROM_TO: FROM_TO+str('_Source'), 'id':'Source'})
df_edges5 = df_edges4.rename(columns={'Target':'id'})
df_edges5 = df_edges5.merge(df_samples, how='left', on='id')
df_edges5 = df_edges5.rename(columns={FROM_TO: FROM_TO+str('_Target'), 'id':'Target'})

#~~~~~~~~~~~~~~~~~~~~~~~
SAMPLE_START = '2'
SAMPLE_END   = '3'

#SAMPLE_START = '3_1_1'
#SAMPLE_END   = '10_1_1'

#SAMPLE_START = '10_1_1'
#SAMPLE_END   = '3_1_5'
#~~~~~~~~~~~~~~~~~~~~~~~
df_edges6 = df_edges5[(df_edges5[FROM_TO+str('_Source')] == SAMPLE_START) & (df_edges5[FROM_TO+str('_Target')] == SAMPLE_END)]
df_edges6_backward= df_edges5[(df_edges5[FROM_TO+str('_Source')] == SAMPLE_END) & (df_edges5[FROM_TO+str('_Target')] == SAMPLE_START)]
products6  = df_edges6['Mass (target)'] - df_edges6['Mass (source)']
products6_backward  = df_edges6_backward['Mass (target)'] - df_edges6_backward['Mass (source)']

for x in trafos:
    TRAFO = x
    globals()['df_edges6_%s' % x] = df_edges6[(df_edges6['Molecular Formula'].str.strip() == TRAFO)]
    globals()['products6_%s' % x]  = globals()['df_edges6_%s' % x]['Mass (target)'] - globals()['df_edges6_%s' % x]['Mass (source)']  
    globals()['df_edges6_backward_%s' % x] = df_edges6_backward[(df_edges6_backward['Molecular Formula'].str.strip() == TRAFO)]
    globals()['products6_backward_%s' % x] = globals()['df_edges6_backward_%s' % x]['Mass (target)'] - globals()['df_edges6_backward_%s' % x]['Mass (source)']

forward = []
for x in trafos:
    a = globals()['products6_%s' % x].shape[0]
    b = int(a)
    forward.append(b)

backward = []
for x in trafos:
    a = globals()['products6_backward_%s' % x].shape[0]
    b = int(a)
    backward.append(b)

df = pd.DataFrame({'Forward': forward,  'Backward': backward}, index=trafos)

fontsize = 12
df.plot(kind='bar', rot=0,  width=.8, color={"Forward": "black", 'Backward': 'w'},  edgecolor='black', stacked=True, fontsize=fontsize)
plt.title(SAMPLE_START+r'$\rightarrow$'+SAMPLE_END+'    ('+str(sample)+', Component ' +str(FILTER_COMPONENTS)+')', fontsize=fontsize+2)
plt.xlabel('Transformations', fontsize=fontsize+2)
plt.ylabel('Frequency', fontsize=fontsize+2)
plt.legend(fontsize=fontsize)
plt.savefig(name+'_Component'+str(FILTER_COMPONENTS)+'_'+SAMPLE_START+SAMPLE_END+'    trafo_bar_chart'+'.png')
plt.show()