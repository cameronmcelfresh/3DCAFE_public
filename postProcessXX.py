import damask
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


result_file = 'grainsID_tensionX.hdf5' #load file
result = damask.Result(result_file)

############################################
####Exract the GRAIN-averaged stress/strain data
############################################

# grid   = damask.Grid.load('grainsID.vti')
# grains = list(np.arange(grid.N_materials))

# data = {g:pd.DataFrame() for g in grains}
# for inc in result.get(['F','P']).values():
#     P = inc['P']
#     F = inc['F']
#     for g in grains:
#         points = grid.material.flatten(order='F')==g
#         P_11 = P[points,0,0].flatten()
#         F_11 = np.broadcast_to(np.average(F[:,0,0]),P_11.shape)
#         x = pd.DataFrame({'F_11':F_11,'P_11':P_11})
#         data[g] = pd.concat((data[g],x),ignore_index=True)

# for g in grains:
#    plot = sns.lineplot(y='P_11',x='F_11',data=data[g])
# fig = plot.get_figure()
# plt.show()


##############################################################
## Integrate over all points to extract the stress/strain data! 
##############################################################

data1 = pd.DataFrame()
P_11List = []
F_11List = []
strainList=[]

#Loop through all increments during which output was printed
for incIndex in result._incs:
    #print("Inc "+ str(incIndex) + "\n")
    resultValue = result.view(increments=incIndex)

    F_extract = resultValue.get(["F"])
    P_extract = resultValue.get(["P"])

    F_13 = F_extract[:,0,2] #extract the 13 component
    F_11 = F_extract[:,0,0] #extract the 11 component
    F_12 = F_extract[:,0,1] #extract the 12 component
    P_11 = P_extract[:,0,0] #extract the xx component

    #Integrate over the whole simulation volume to find the average
    P_11List.append(np.average(P_11)) 
    #F_11List.append(np.average(F_11)-1)

    strainList.append( (np.average(np.multiply(F_13,F_13)+  np.multiply(F_12,F_12)+np.multiply(F_11,F_11))-1)/2 )

data1 = pd.DataFrame({'strain':strainList,'stress':P_11List})

#Plotting the data
#plot1 = sns.lineplot(y='stress',x='strain',data=data1)
#fig1 = plot1.get_figure()
#plt.show()

data1.to_csv('stress_strainXX.csv', index=False) #save data

