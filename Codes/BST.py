import numpy as np
import matplotlib.pyplot as plt

#calculate extrapolated scaling exponents
def Tau_ext(tau,omega):
    L_lst,q_lst,n_lst,omega_lst = Parameters()
    nl = np.size(L_lst)
    nq = np.size(q_lst)
    tau_ext = np.zeros(nq)
    D = np.zeros(nq)
    h = 1/L_lst
    T = np.zeros((nl-2,nl)) #table of extrapolants, with indices corresponding to m and L, respectively
    for q in range(nq):
        T[1] = tau[:,q]
        if q_lst[q] == 0: #avoid q=0 limit
            tau_ext[q] = -1
        elif q_lst[q] == 1: #avoid q=1 limit
            tau_ext[q] = 0
        else:
            for i in range(2,nl-2): #column of BST table
                for j in range(i-1,nl-i+1): #row of BST table
                    T[i,j] = T[i-1,j+1]+(T[i-1,j+1]-T[i-1,j])/((h[j]/h[i+j-1])**omega*(1-(T[i-1,j+1]-T[i-1,j])/(T[i-1,j+1]-T[i-2,j+1]))-1)
            tau_ext[q] = T[3,int(nl/2)] #extrapolant taken as one of the elements in final column
        D[q] = 2*(T[3,3]-T[3,2]) #difference between a pair of adjacent elements
    return tau_ext,D

#plot difference between adjacent elements in BST table
def Difference(tau):
    L_lst,q_lst,n_lst,omega_lst = Parameters()
    ng = 20
    nq = np.size(q_lst)
    omega_lst = np.linspace(3,6,ng)
    D = np.zeros((ng,nq))
    for i in range(ng):
        D[i] = Tau_ext(tau,omega_lst[i])[1]
    for i in range(1,nq,2):
        plt.plot(omega_lst,D[:,i],'-x',label=r'$q={}$'.format(round(q_lst[i],1)))
    plt.xlabel(r'$\omega$')
    plt.ylabel(r'$\epsilon(\omega)$')
    plt.ylim([-1,1])
    plt.legend(ncol=2)
    plt.show()

def Main():
    omega_lst = Parameters()[3]
    tau_ext_1 = Tau_ext(tau_1,omega_lst[0])[0]
    tau_ext_2 = Tau_ext(tau_2,omega_lst[1])[0]
    tau_ext_3 = Tau_ext(tau_3,omega_lst[2])[0]
    tau_ext_4 = Tau_ext(tau_4,omega_lst[3])[0]
    Plot(tau_ext_1,tau_ext_2,tau_ext_3,tau_ext_4)

def Plot(tau_ext_1,tau_ext_2,tau_ext_3,tau_ext_4):
    L_lst,q_lst,n_lst,omega_lst = Parameters()
    nq = np.size(q_lst)
    label = np.array(["a","b","c","d"])
    tau_ext = np.zeros((4,nq))
    tau_ext[0] = tau_ext_1
    tau_ext[1] = tau_ext_2
    tau_ext[2] = tau_ext_3
    tau_ext[3] = tau_ext_4
    plt.figure(figsize=(8,8),dpi=80)
    for i in range(4):
        plt.subplot(2,2,i+1)
        plt.plot(q_lst,tau_ext[i],'-x')
        plt.axhline(y=0,linestyle='--',xmin=0.35,color='k')
        plt.annotate(r'$\gamma T_0/2\pi={}$'.format(n_lst[i]),(0.05,0.92),xycoords='axes fraction')
        plt.annotate('CUE',(0.88,0.86),xycoords='axes fraction')
        plt.annotate('Loc.',(0.9,0.28),xycoords='axes fraction')
        plt.annotate('({})'.format(label[i]),(0.92,0.03),xycoords='axes fraction')
        plt.plot(q_lst,q_lst-1,'k--')
        plt.ylim([-1,2])
        ax = plt.gca()
        if i < 2:
            ax.set_xticklabels([])
        else:
            plt.xlabel(r'$q$')
        if i == 0 or i == 2:
            plt.ylabel(r'$\tau_q$')
        else:
            ax.set_yticklabels([])
        ax.set_aspect(1)
    plt.tight_layout()
    plt.show()

def Parameters():
    L_lst = np.array([6,8,10,12,14,16]) #system size
    q_lst = np.linspace(0,3,16) #fractal dimension parameter
    n_lst = np.array([0.5,1.0,1.5,2.0]) #resonance parameter, equal to γT_0/2π
    omega_lst = np.array([10,1.8,10,4.1]) #free parameter for BST algorithm
    return L_lst,q_lst,n_lst,omega_lst

tau_1 = np.load('Stark_Files_2/Tau/Data/tau_0.5_weak.npy')
tau_2 = np.load('Stark_Files_2/Tau/Data/tau_1.0_weak.npy')
tau_3 = np.load('Stark_Files_2/Tau/Data/tau_1.5_weak.npy')
tau_4 = np.load('Stark_Files_2/Tau/Data/tau_2.0_weak.npy')

Main()