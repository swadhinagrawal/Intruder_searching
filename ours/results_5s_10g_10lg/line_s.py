import matplotlib.pyplot as plt
import numpy as np
import pickle as pkl
import os
from matplotlib import rc,rcParams
rc('font', weight='bold',size=20)
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

path = os.path.realpath(os.path.dirname(__file__)) 
data_files = sorted(os.listdir(path))
files = []
for i in data_files:
    if '.' not in i and 'sfcg' not in i:
        files.append(i)
 
def mean_std(data):
    mu = np.sum(data,axis=1)/data.shape[1]
    std = []
    for i in range(len(mu)):    
        std.append((np.sum((data[i]-mu[i])**2)/data.shape[1])**(0.5))
    return mu,np.array(std)

all_data = []
raw_data = []
for f in files:
    file = open(path+'/'+f,'rb')
    data = pkl.load(file)
    raw_data.append(np.array(data))
    num_robots = np.array(data)[:,0,0]
    search_times = np.array(data)[:,:,1]
    mu_t, std_t = mean_std(search_times)
    all_data.append(np.array([num_robots,mu_t,std_t]))

colors =['blue','orange','#FE6100','green']
labels = ['Benchmark','RS','CRS','SFC']
itype = [', M-I',', S-I']

line,ax1 = plt.subplots()
counter = 0
for i in range(len(all_data)):
    w = 0.5
    if i != 0 and i%2 == 0:
        counter += 1
    if (i+1)%2 == 0:
        ax1.plot(all_data[i][0],all_data[i][1], color=colors[counter],label=labels[int(i/2)],linewidth = 3)

ax1.set_xticks(range(0,39,6))
ax1.set_yticks(range(0,70,10))
ax1.set_title('Static Intruder',fontdict=dict(weight='bold',size=20))
ax1.legend(prop=dict(weight='bold',size=20),frameon=False)
plt.xlabel('Number of Robots',fontdict=dict(weight='bold',size=20))
plt.ylabel('Search steps',fontdict=dict(weight='bold',size=20))
line.savefig(path+'/line_s.pdf',format = "pdf",bbox_inches="tight",pad_inches=0.2)
plt.show()

