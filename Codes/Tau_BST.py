import numpy as np
import matplotlib.pyplot as plt

#calculate extrapolated scaling exponents
def Tau_ext(tau,omega):
    L_lst,q_lst,n_lst = Parameters()
    nl = np.size(L_lst)
    nq = np.size(q_lst)
    tau_ext = np.zeros(nq)
    epsilon = np.zeros(nq)
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
        epsilon[q] = 2*(T[2,4]-T[2,3]) #difference between pairs of adjacent elements
    return tau_ext,epsilon

#check where solutions are for omega
def Check(tau,n):
    L_lst,q_lst,n_lst = Parameters()
    ng = 20 #omega increments
    omega_min = 1 #minimum omega value
    omega_max = 10 #maximum omega value
    omega_lst = np.linspace(omega_min,omega_max,ng)
    nq = np.size(q_lst) #q increments
    tau_ext = np.zeros((ng,nq))
    epsilon = np.zeros((ng,nq))
    for i in range(ng):
        tau_ext[i],epsilon[i] = Tau_ext(tau,omega_lst[i]) #obtain extrapolated scaling exponents for all omega and q
    plt.figure(figsize=(7,7),dpi=80)
    plt.rcParams.update({'font.size': 18})
    plt.axhline([0],linestyle='dashed',color='r',label='Loc.')
    a = np.where(np.abs(epsilon[:,n])==np.min(np.abs(epsilon[:,n])))[0][0]
    plt.axvline([omega_lst[a]],linestyle='dashed',color='k')
    plt.annotate(r'$\omega = {}$'.format(round(omega_lst[a],2)),(0.16,0.1),xycoords='axes fraction')
    #plt.axhline([q_lst[n]-1],linestyle='dashed',color='b',label='CUE')
    for j in range(1,nq,2):
        plt.plot(omega_lst,epsilon[:,j],'-x',label=r'$q={}$'.format({round(q_lst[j],1)}))
    #plt.plot(omega_lst,tau_ext[:,n],label=r'$\tau_{}$'.format({round(q_lst[n],1)}))
    #plt.annotate(r'$\tau_{} = {}$'.format(({round(q_lst[n],2)}),round(tau_ext[a,n],1)),(0.5,0.05),xycoords='axes fraction')
    plt.xlabel(r'$\omega$')
    plt.ylabel(r'$\epsilon_1^{(2)}(\omega)$')
    plt.ylim([-0.5,0.5])
    plt.legend(loc='lower right',ncol=2,fontsize=14)
    ax = plt.gca()
    ax.set_aspect(7)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    plt.show()

#check how omega changes with q
def Drift(tau):
    L_lst,q_lst,n_lst = Parameters()
    ng = 100 #omega increments
    omega_min = 0.1 #minimum omega value
    omega_max = 10 #maximum omega value
    omega_lst = np.linspace(omega_min,omega_max,ng)
    nq = np.size(q_lst) #q increments
    tau_ext = np.zeros((ng,nq))
    epsilon = np.zeros((ng,nq))
    a = np.zeros(nq)
    for i in range(ng):
        tau_ext[i],epsilon[i] = Tau_ext(tau,omega_lst[i]) #obtain extrapolated scaling exponents for all omega and q
    for i in range(nq):
        a[i] = omega_lst[np.where(epsilon[:,i]==np.min(epsilon[:,i]))[0][0]]
    plt.plot(q_lst,a,'x')
    plt.xlabel(r'$q$')
    plt.ylabel(r'$\omega(q)$')
    plt.show()

def Main():
    tau_ext_1 = Tau_ext(tau_1,100)[0]
    tau_ext_2 = Tau_ext(tau_2,1.8)[0]
    tau_ext_3 = Tau_ext(tau_3,100)[0]
    tau_ext_4 = Tau_ext(tau_4,4.1)[0]
    Plot(tau_ext_1,tau_ext_2,tau_ext_3,tau_ext_4)

def Plot(tau_ext_1,tau_ext_2,tau_ext_3,tau_ext_4):
    L_lst,q_lst,n_lst = Parameters()
    nq = np.size(q_lst)
    label = np.array(["a","b","c","d"])
    tau_ext = np.zeros((6,nq))
    tau_ext[0] = tau_ext_1
    tau_ext[1] = tau_ext_2
    tau_ext[2] = tau_ext_3
    tau_ext[3] = tau_ext_4
    plt.figure(figsize=(7,7),dpi=80)
    plt.rcParams.update({'font.size': 14})
    for i in range(4):
        plt.subplot(2,2,i+1)
        plt.plot(q_lst,tau_ext[i],'-x')
        plt.axhline(y=0,linestyle='--',xmin=0.35,color='k')
        plt.annotate(r'$\gamma T_0/2\pi={}$'.format(n_lst[i]),(0.05,0.90),xycoords='axes fraction')
        plt.annotate('CUE',(0.85,0.82),xycoords='axes fraction')
        plt.annotate('Loc.',(0.85,0.26),xycoords='axes fraction')
        plt.annotate('({})'.format(label[i]),(0.88,0.03),xycoords='axes fraction')
        plt.plot(q_lst,q_lst-1,'k--')
        plt.ylim([-1,2])
        if i == 0 or i == 2:
            plt.ylabel(r'$\tau_q$')
        else:
            ax = plt.gca()
            ax.set_yticklabels([])
        if i > 1:
            plt.xlabel(r'$q$')
        else:
             ax = plt.gca()
             ax.set_xticklabels([])
        ax = plt.gca()
        ax.set_aspect(1)
    plt.tight_layout()
    plt.subplots_adjust(hspace=-0.1)
    plt.show()

def Parameters():
    L_lst = np.array([6,8,10,12,14,16]) #system size
    q_lst = np.linspace(0,3,16) #fractal dimension parameter
    n_lst = np.array([0.5,1.0,1.5,2.0]) #resonance parameter, equal to γT_0/2π
    return L_lst,q_lst,n_lst

tau_1 = np.load('Stark_Files_2/Tau/Data/tau_0.5_weak.npy')
tau_2 = np.load('Stark_Files_2/Tau/Data/tau_1.0_weak.npy')
tau_3 = np.load('Stark_Files_2/Tau/Data/tau_1.5_weak.npy')
tau_4 = np.load('Stark_Files_2/Tau/Data/tau_2.0_weak.npy')

Main()