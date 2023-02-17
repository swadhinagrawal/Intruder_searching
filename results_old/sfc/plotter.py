import numpy as np
import matplotlib.pyplot as plt
import copy as cp
import pickle as pkl
import os
np.random.seed(748657)


# Results plotter
run_vs_search_time = 0
if run_vs_search_time:
    path = os.path.realpath(os.path.dirname(__file__))
    B1 = open(path+'/results/benchmark/MRIS_s_bench', 'rb')
    B2 = open(path+'/results/benchmark/MRIS_d_bench', 'rb')
    M1 = open(path+'/results/sfc/MRIS_s_sfc', 'rb')# Not considering until k_min searcher case
    M2 = open(path+'/results/sfc/MRIS_d_sfc', 'rb')
    M3 = open(path+'/results/sfc/MRIS_s_gsfc', 'rb')
    M4 = open(path+'/results/sfc/MRIS_d_gsfc', 'rb')
    M5 = open(path+'/results/random/MRIS_s_random', 'rb')
    M6 = open(path+'/results/random/MRIS_d_random', 'rb')
    M7 = open(path+'/results/randomTA/MRIS_s_rta', 'rb')
    M8 = open(path+'/results/randomTA/MRIS_d_rta', 'rb')
    
    obj = pkl.load(file)

    plt.ioff()
    fig,ax = plt.subplots()
    
    ax.scatter(range(100),obj,linewidth=4,label='Data')
    ax.plot(range(100),[np.sum(obj)/len(obj) for i in range(100)],linewidth=4,label='Mean: '+str(np.sum(obj)/len(obj)),color='orange')
    t_font = {'weight': 'bold',
        'size': 15}
    ax_font = {'weight': 'bold',
        'size': 20}
    leg_font = {'weight': 'bold',
        'size': 12}
    plt.title('1RDIS: Random Search',fontdict=t_font)
    plt.xlabel('Runs',fontdict=ax_font)
    plt.ylabel('Search Time',fontdict=ax_font)
    plt.ticklabel_format(style='sci',scilimits=(0,3))
    plt.xticks(fontsize=15,fontweight='bold')
    plt.yticks(fontsize=15,fontweight='bold')
    plt.tight_layout()
    plt.legend(prop=leg_font)
    plt.savefig(path+'/results/1RIS_d_random_plot.pdf',format = "pdf",bbox_inches="tight",pad_inches=0)
    plt.show()

swarch_no_vs_search_time = 0
if swarch_no_vs_search_time:
    path = os.getcwd()
    file = open(path+'/results/MRIS_s_random', 'rb')
    obj = pkl.load(file)

    plt.ioff()
    fig,ax = plt.subplots()
    
    ax.plot(obj,range(1,10),linewidth=4,label='Static Intruder')
    # ax.plot(range(10),obj,linewidth=4)
    # ax.plot(range(9),[np.sum(obj)/len(obj) for i in range(10)],linewidth=4)
    t_font = {'weight': 'bold',
        'size': 15}
    ax_font = {'weight': 'bold',
        'size': 20}
    leg_font = {'weight': 'bold',
        'size': 12}
    plt.title('MRIS: Random Search',fontdict=t_font)
    plt.xlabel('Average search time',fontdict=ax_font)
    plt.ylabel('Search Number',fontdict=ax_font)
    plt.xticks(fontsize=15,fontweight='bold')
    plt.yticks(fontsize=15,fontweight='bold')
    plt.tight_layout()
    file.close()
    # plt.savefig(path+'/results/MRIS_random_s_plot.png')
    # plt.show()
    file = open(path+'/results/MRIS_d_random', 'rb')
    obj = pkl.load(file)

    plt.ioff()
    # fig,ax = plt.subplots()
    
    ax.plot(obj,range(1,10),linewidth=4,label='Dynamic Intruder')
    # ax.plot(range(10),obj,linewidth=4)
    # ax.plot(range(9),[np.sum(obj)/len(obj) for i in range(10)],linewidth=4)
    t_font = {'weight': 'bold',
        'size': 15}
    ax_font = {'weight': 'bold',
        'size': 20}
    leg_font = {'weight': 'bold',
        'size': 12}
    plt.title('MRIS: Random Search',fontdict=t_font)
    plt.xlabel('Average search time',fontdict=ax_font)
    plt.ylabel('Search Number',fontdict=ax_font)
    plt.xticks(fontsize=15,fontweight='bold')
    plt.yticks(fontsize=15,fontweight='bold')
    plt.tight_layout()
    plt.legend(prop=leg_font)
    plt.savefig(path+'/results/MRIS_random_plot.pdf',format = "pdf",bbox_inches="tight",pad_inches=0)
    plt.show()
