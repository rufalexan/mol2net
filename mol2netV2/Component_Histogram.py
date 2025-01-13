#################################
####          COMPONENT HISTOGRAM
#################################
plt.figure()
pd.DataFrame(xNODES0['component'].value_counts(normalize=True)).plot(kind='bar', rot=0, ylabel='', legend=False, color='black', width=1)
#pd.DataFrame(xNODES0['component'].value_counts(normalize=False)).sort_index().plot(kind='bar', rot=0, ylabel='', legend=False, color='black', width=1)
#plt.title(name+', Component '+str(FILTER_COMPONENTS), wrap=True)     
plt.title('Minimal set [H2, O, CO, NH3]')
plt.xticks(np.arange(0, len(df_comp_freq), 20), fontsize=15)                      ## ice
####plt.xticks(np.arange(0, len(df_comp_freq), 100), fontsize=15)                      ## paris
plt.yticks(fontsize=15)
plt.yscale('log')
plt.xlabel('Component number', fontsize=15)
###plt.ylabel('Frequency (log scale)')
plt.ylabel('Frequency (log scale, normalized)', fontsize=15)
###plt.ylabel('Frequency (normalized)')
#plt.savefig(name+', Comp'+str(FILTER_COMPONENTS)+ '    hist_components'+'.svg')
plt.savefig(name+', Comp'+str(FILTER_COMPONENTS)+ '    hist_components'+'.png')

plt.show()
