

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
plt.title(name+'\n, Comp'+str(FILTER_COMPONENTS)+', Trafos'+str(trafos), wrap=True, fontsize=12)
plt.scatter(x,y,c='blue')
for i in range(len(x)):
    #plt.annotate(text[i], (x[i], y[i] + 0.2), fontsize=10)

    plt.xscale('log',base=10) 
    plt.yscale('log',base=10) 
plt.xlabel(x_name)
plt.ylabel(y_name)
plt.savefig(name+'_'+str(FILTER_COMPONENTS)+'_'+str(trafos)+'    degreeDistri'+'.png')
plt.show()



