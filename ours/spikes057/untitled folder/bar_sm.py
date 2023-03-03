import matplotlib.pyplot as plt
import numpy as np
import pickle as pkl
import os
from matplotlib import rc,rcParams
rc('font', size=14)
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

path = os.path.realpath(os.path.dirname(__file__)) 
data_files = sorted(os.listdir(path))
files = []
for i in data_files:
    if '.' not in i and 'untitled' not in i:
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
# all_data = np.array(all_data)
bar,ax = plt.subplots()
colors =['blue','orange','#FE6100','green']
labels = ['Baseline','RS','CRS','SFC']
itype = [', '+r'$\beta=0$',', '+r'$\beta=5$',', '+r'$\beta=7$']
bench = [3,15,22,31,40,49] # 10, 46, 67, 94, 121, 148
sfc = [0,36,57,84,111,138]
random = [3,15,22,31,40,49]
crs_min = [3,15,22,31,40,49]

datas = [bench,random,crs_min,sfc]#,gsfc,random,crs_min]

sfc0 = [9,36,57,84,111,138]
datas1 = [bench,random,crs_min,sfc0]#,gsfc,random,crs_min]

ticks = []
boxes_m = []
boxes_s = []
boxes = []
for k in range(1):
    counter = 0
    for i in range(len(all_data)):
        w = 3
        if i != 0 and i%3 == 0:
            counter += 1
        
        if k == 0 and (i)%3 == 1:
            box1 = ax.boxplot([np.array(raw_data[i])[datas[counter][k],:,1]],medianprops=dict(color="red",linewidth=1),manage_ticks=False,widths=[w], positions=[all_data[i][0,datas[counter][k]]+i*4],patch_artist=True,showfliers=False)#,labels=[labels[int(i/2)]])#+ itype[0])
            box1['boxes'][0].set_facecolor(colors[counter])
            box1['boxes'][0].set_edgecolor(colors[counter])
            box1['boxes'][0].set_alpha(0.5)
            boxes.append(box1)
            boxes_m.append(box1)
            ticks.append(all_data[i][0,datas[counter][k]])
            print('5:',all_data[i][0,datas[counter][k]])
            
        elif (i)%3 == 1:
            box1 = ax.boxplot(np.array(raw_data[i])[datas[counter][k],:,1],medianprops=dict(color="red",linewidth=1),manage_ticks=False,widths=[w],positions=[all_data[i][0,datas[counter][k]]+i*4],patch_artist=True,showfliers=False)#+ itype[0])
            box1['boxes'][0].set_facecolor(colors[counter])
            box1['boxes'][0].set_edgecolor(colors[counter])
            box1['boxes'][0].set_alpha(0.5)
            ticks.append(all_data[i][0,datas[counter][k]])

        if k == 0 and (i)%3 == 0:
            box1 = ax.boxplot([np.array(raw_data[i])[datas1[counter][k],:,1]],medianprops=dict(color="red",linewidth=1),manage_ticks=False,widths=[w], positions=[all_data[i][0,datas1[counter][k]]+i*4],patch_artist=True,showfliers=False)#,labels=[labels[int(i/2)]])#+ itype[0])
            box1['boxes'][0].set_facecolor('white')
            box1['boxes'][0].set_edgecolor(colors[counter])
            box1['boxes'][0].set_linewidth(2)
            box1['boxes'][0].set_alpha(0.5)
            # box1['boxes'][0].set(hatch='////')
            boxes.append(box1)
            boxes_s.append(box1)
            print('0:',all_data[i][0,datas1[counter][k]])
            ticks.append(all_data[i][0,datas1[counter][k]])
        elif k != 0 and (i)%3 == 0:
            box1 = ax.boxplot(np.array(raw_data[i])[datas1[counter][k],:,1],medianprops=dict(color="red",linewidth=1),manage_ticks=False,widths=[w],positions=[all_data[i][0,datas1[counter][k]]+i*4],patch_artist=True,showfliers=False)#+ itype[0])
            box1['boxes'][0].set_facecolor('white')
            box1['boxes'][0].set_edgecolor(colors[counter])
            box1['boxes'][0].set_linewidth(2)
            # box1['boxes'][0].set(hatch='////')
            box1['boxes'][0].set_alpha(0.5)
            ticks.append(all_data[i][0,datas1[counter][k]])
            
        if k == 0 and (i)%3 == 2:
            box1 = ax.boxplot([np.array(raw_data[i])[datas[counter][k],:,1]],medianprops=dict(color="red",linewidth=1),manage_ticks=False,widths=[w],positions=[all_data[i][0,datas[counter][k]]+i*4],patch_artist=True,showfliers=False)#,labels=[labels[int(i/2)]])#+ itype[0])
            box1['boxes'][0].set_facecolor('white')
            box1['boxes'][0].set_edgecolor(colors[counter])
            box1['boxes'][0].set_linewidth(2)
            box1['boxes'][0].set_alpha(0.5)
            box1['boxes'][0].set(hatch='////')
            boxes.append(box1)
            boxes_s.append(box1)
            ticks.append(all_data[i][0,datas[counter][k]])
            print('7:',all_data[i][0,datas[counter][k]])
        elif k != 0 and (i)%3 == 2:
            box1 = ax.boxplot(np.array(raw_data[i])[datas[counter][k],:,1],medianprops=dict(color="red",linewidth=1),manage_ticks=False,widths=[w],positions=[all_data[i][0,datas[counter][k]]+i*4],patch_artist=True,showfliers=False)#+ itype[0])
            box1['boxes'][0].set_facecolor('white')
            box1['boxes'][0].set_edgecolor(colors[counter])
            box1['boxes'][0].set_linewidth(2)
            box1['boxes'][0].set(hatch='////')
            box1['boxes'][0].set_alpha(0.5)
            ticks.append(all_data[i][0,datas[counter][k]])
ax.set_ylim(-0.2)
ax.set_yticks(range(0,160,30))
ax.set_title('Static Intruder',fontdict=dict(size=14))
# print([labels[int(i/2)] for i in range(0,len(files)) if (i+1)%2 != 0])
# ax.legend([b["boxes"][0] for b in boxes_m]+[b["boxes"][0] for b in boxes_s], [labels[int(i/2)]+itype[0] for i in range(0,len(files)) if (i+1)%2 != 0]+[labels[int(i/2)]+itype[1] for i in range(0,len(files)) if (i+1)%2 != 0],prop=dict(weight='bold',size=18),frameon=False,bbox_to_anchor=(1, 1))
plt.xlabel('Number of Robots',fontdict=dict(size=14))
plt.ylabel('Search steps',fontdict=dict(size=14))
# plt.xticks(rotation=45)

ax.set_xticks(np.unique(ticks))
bar.savefig(path+'/bar_s.pdf',format = "pdf",bbox_inches="tight",pad_inches=0.2)
plt.show()

