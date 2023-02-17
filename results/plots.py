import matplotlib.pyplot as plt
import numpy as np
import pickle as pkl
import os
from matplotlib import rc,rcParams
rc('font', weight='bold',size=12)
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

path = os.path.realpath(os.path.dirname(__file__)) 
data_files = sorted(os.listdir(path))
files = []
for i in data_files:
    if '.' not in i:
        files.append(i)
 
def mean_std(data):
    mu = np.sum(data,axis=1)/data.shape[1]
    std = []
    for i in range(len(mu)):    
        std.append((np.sum((data[i]-mu[i])**2)/data.shape[1])**(0.5))
    return mu,np.array(std)

all_data = []
for f in files:
    file = open(path+'/'+f,'rb')
    data = pkl.load(file)
    num_robots = np.array(data)[:,0,0]
    search_times = np.array(data)[:,:,1]
    mu_t, std_t = mean_std(search_times)
    all_data.append(np.array([num_robots,mu_t,std_t]))
# all_data = np.array(all_data)
bar,ax = plt.subplots()
colors =['blue','orange','pink','green']
bench = [16,43,60] # 81, 216,301
sfc = [1,136,221]
gsfc = [1,136,221]
datas = [bench,sfc,gsfc]
ticks = []
for k in range(3):
    counter = 0
    for i in range(len(all_data)):
        w = 10
        if i != 0 and i%2 == 0:
            counter += 1
        ticks.append(all_data[i][0][datas[counter][k]])
        if k == 0 and (i+1)%2 != 0:
            ax.bar([all_data[i][0,datas[counter][k]]+i*11],[all_data[i][1,datas[counter][k]]],w,yerr=[all_data[i][2,datas[counter][k]]], color=colors[counter],ecolor='black',alpha=0.5,capsize=5,label=files[i])
        elif k == 0 and (i+1)%2 == 0:
            ax.bar([all_data[i][0,datas[counter][k]]+i*11],[all_data[i][1,datas[counter][k]]],w,yerr=[all_data[i][2,datas[counter][k]]], color=colors[counter],ecolor='black',alpha=0.5,capsize=5,label=files[i],hatch = '..')    
        elif k != 0 and (i+1)%2 == 0:
            ax.bar([all_data[i][0,datas[counter][k]]+i*11],[all_data[i][1,datas[counter][k]]],w,yerr=[all_data[i][2,datas[counter][k]]], color=colors[counter],ecolor='black',alpha=0.5,capsize=5,hatch = '..')    
        else:
            ax.bar([all_data[i][0,datas[counter][k]]+i*11],[all_data[i][1,datas[counter][k]]],w,yerr=[all_data[i][2,datas[counter][k]]], color=colors[counter],ecolor='black',alpha=0.5,capsize=5)

plt.xlabel('Number of Robots')
plt.ylabel('Search Time')
plt.xticks(ticks)
plt.legend()
plt.show()
