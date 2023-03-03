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
    if '.' not in i:
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

bar,ax = plt.subplots()
colors =['blue','orange','#FE6100','green','#9A03A2']
labels = ['Benchmark','RS','CRS','SFC','SFC-G']
itype = [', M-I',', S-I']

bench = [15,22,31,40,49] # 46, 67, 94, 121, 148
sfc = [36,57,84,111,138]
gsfc = [2,23,50,77,104]
random = [15,22,31,40,49]
crs_min = [15,22,31,40,49]
datas = [bench,random,crs_min,sfc,gsfc]
ticks = []
boxes = []
for k in range(5):
    counter = 0
    for i in range(len(all_data)):
        w = 2
        if i != 0 and i%2 == 0:
            counter += 1
        ticks.append(all_data[i][0,datas[counter][k]])
        if k == 0 and (i+1)%2 != 0:
            box1 = ax.boxplot([np.array(raw_data[i])[datas[counter][k],:,1]],medianprops=dict(color="red",linewidth=1),manage_ticks=False,widths=[w], positions=[all_data[i][0,datas[counter][k]]+i*1.5],patch_artist=True,showfliers=False)#,labels=[labels[int(i/2)]])#+ itype[0])
            box1['boxes'][0].set_facecolor(colors[counter])
            box1['boxes'][0].set_edgecolor(colors[counter])
            box1['boxes'][0].set_alpha(0.5)
            boxes.append(box1)
            
        elif (i+1)%2 != 0:
            box1 = ax.boxplot(np.array(raw_data[i])[datas[counter][k],:,1],medianprops=dict(color="red",linewidth=1),manage_ticks=False,widths=[w],positions=[all_data[i][0,datas[counter][k]]+i*1.5],patch_artist=True,showfliers=False)#+ itype[0])
            box1['boxes'][0].set_facecolor(colors[counter])
            box1['boxes'][0].set_edgecolor(colors[counter])
            box1['boxes'][0].set_alpha(0.5)

ax.set_ylim(-0.2)
ax.set_yticks(range(0,21,4))
ax.set_title('Moving Intruder',fontdict=dict(weight='bold',size=20))
print([labels[int(i/2)] for i in range(0,len(files)) if (i+1)%2 != 0])
ax.legend([b["boxes"][0] for b in boxes], [labels[int(i/2)] for i in range(0,len(files)) if (i+1)%2 != 0],prop=dict(weight='bold',size=18),frameon=False,bbox_to_anchor=(1, 1))
plt.xlabel('Number of Robots',fontdict=dict(weight='bold',size=20))
plt.ylabel('Search steps',fontdict=dict(weight='bold',size=20))
# plt.xticks(rotation=45)

ax.set_xticks(np.unique(ticks))
bar.savefig(path+'/bar_m.pdf',format = "pdf",bbox_inches="tight",pad_inches=0.2)
plt.show()
