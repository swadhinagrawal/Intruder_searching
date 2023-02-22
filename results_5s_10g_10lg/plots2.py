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
for f in files:
    file = open(path+'/'+f,'rb')
    data = pkl.load(file)
    num_robots = np.array(data)[:,0,0]
    search_times = np.array(data)[:,:,1]
    mu_t, std_t = mean_std(search_times)
    all_data.append(np.array([num_robots,mu_t,std_t]))
# all_data = np.array(all_data)
bar,ax = plt.subplots()
colors =['blue','orange','slateblue','green','brown']
labels = ['Benchmark','RS','CRS','SFC','SFC-G']
itype = [', M-I',', S-I']
bench = [23,31,37] # 24, 32, 38
sfc = [14,22,28]
gsfc = [0,8,14]
random = [23,31,37]
crs = [22,30,36]
datas = [bench,random,crs,sfc,gsfc]
ticks = []
for k in range(1):
    counter = 0
    for i in range(len(all_data)):
        w = 0.5
        if i != 0 and i%2 == 0:
            counter += 1
        ticks.append(all_data[i][0][datas[counter][k]])
        if k == 0 and (i+1)%2 != 0:
            y_errormin = [all_data[i][1,datas[counter][k]]]
            y_errormax = [all_data[i][2,datas[counter][k]]]
            y_error = [y_errormin, y_errormax]
            ax.bar([all_data[i][0,datas[counter][k]]+i*0.5],[all_data[i][1,datas[counter][k]]],w, color=colors[counter],ecolor='black',alpha=0.5,label=labels[int(i/2)] + itype[0])
            ax.errorbar([all_data[i][0,datas[counter][k]]+i*0.5],[all_data[i][1,datas[counter][k]]],yerr=y_error,color="black",capsize=5)
        elif (i+1)%2 != 0:
            y_errormin = [all_data[i][1,datas[counter][k]]]
            y_errormax = [all_data[i][2,datas[counter][k]]]
            y_error = [y_errormin, y_errormax]
            ax.bar([all_data[i][0,datas[counter][k]]+i*0.5],[all_data[i][1,datas[counter][k]]],w, color=colors[counter],ecolor='black',alpha=0.5)
            ax.errorbar([all_data[i][0,datas[counter][k]]+i*0.5],[all_data[i][1,datas[counter][k]]],yerr=y_error,color="black",capsize=5)
        if k == 0 and (i+1)%2 == 0:
            y_errormin = [all_data[i][1,datas[counter][k]]]
            y_errormax = [all_data[i][2,datas[counter][k]]]
            y_error = [y_errormin, y_errormax]
            ax.bar([all_data[i][0,datas[counter][k]]+i*0.5],[all_data[i][1,datas[counter][k]]],w,color=colors[counter],ecolor='black',alpha=0.5,label=labels[int(i/2)]+ itype[1],hatch = '////')    
            ax.errorbar([all_data[i][0,datas[counter][k]]+i*0.5],[all_data[i][1,datas[counter][k]]],yerr=y_error,color="black",capsize=5)
        elif k != 0 and (i+1)%2 == 0:
            y_errormin = [all_data[i][1,datas[counter][k]]]
            y_errormax = [all_data[i][2,datas[counter][k]]]
            y_error = [y_errormin, y_errormax]
            ax.bar([all_data[i][0,datas[counter][k]]+i*0.5],[all_data[i][1,datas[counter][k]]],w, color=colors[counter],ecolor='black',alpha=0.5,hatch = '////')    
            ax.errorbar([all_data[i][0,datas[counter][k]]+i*0.5],[all_data[i][1,datas[counter][k]]],yerr=y_error,color="black",capsize=5)
        
ax.set_ylim(-0.2)
ax.set_title('Static vs Moving Intruder',fontdict=dict(weight='bold',size=20))
ax.legend(prop=dict(weight='bold',size=18),frameon=False,bbox_to_anchor=(1, 1))
plt.xlabel('Number of Robots',fontdict=dict(weight='bold',size=20))
plt.ylabel('Search steps',fontdict=dict(weight='bold',size=20))
ax.set_xticks(np.unique(ticks))
bar.savefig(path+'/bar_sm.pdf',format = "pdf",bbox_inches="tight",pad_inches=0.2)
plt.show()
# line,ax1 = plt.subplots()


# counter = 0
# for i in range(len(all_data)):
#     w = 0.5
#     if i != 0 and i%2 == 0:
#         counter += 1
#     ticks.append(all_data[i][0][datas[counter][k]])
#     if (i+1)%2 != 0:
#         ax1.plot(all_data[i][0],all_data[i][1], color=colors[counter],label=labels[int(i/2)],linestyle='-', linewidth=3)#,linewidth = counter+1)
#     # if (i+1)%2 == 0:
#     #     ax1.plot(all_data[i][0],all_data[i][1], color=colors[counter],label=labels[int(i/2)],linestyle='--',linewidth = 3)
# ax1.set_title('Moving Intruder',fontdict=dict(weight='bold',size=20))
# ax1.legend(prop=dict(weight='bold',size=20),frameon=False,bbox_to_anchor=(1, 1))
# plt.xlabel('Number of Robots',fontdict=dict(weight='bold',size=20))
# plt.ylabel('Search steps',fontdict=dict(weight='bold',size=20))
# line.savefig(path+'/line_m.pdf',format = "pdf",bbox_inches="tight",pad_inches=0.2)
# plt.show()

