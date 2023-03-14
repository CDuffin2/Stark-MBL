import numpy as np
import matplotlib.pyplot as plt
import time

def Binom(n,k):
    ret = 1
    n = int(n)
    k = int(k)
    for i in range(min(k,n-k)):
        ret *= n-i
        ret /= i+1
    return int(ret)

def Load(a,b):
    L_lst = Parameters()[3]
    q_lst = Parameters()[4]
    nl = np.size(L_lst)
    nq = np.size(q_lst)
    tau = np.zeros((nl,nq))
    tau_8 = np.load('Stark_Files_2/Tau/Data_2/Tau_L=8.npy')
    tau_10 = np.load('Stark_Files_2/Tau/Data_2/Tau_L=10.npy')
    tau_12 = np.load('Stark_Files_2/Tau/Data_2/Tau_L=12.npy')
    tau_14 = np.load('Stark_Files_2/Tau/Data_2/Tau_L=14.npy')
    tau_16 = np.load('Stark_Files_2/Tau/Data_2/Tau_L=16.npy')
    tau[0] = tau_8[a,b]
    tau[1] = tau_10[a,b]
    tau[2] = tau_12[a,b]
    tau[3] = tau_14[a,b]
    tau[4] = tau_16[a,b]
    return tau

#calculate extrapolated scaling exponents
def Tau_ext(tau,omega,q):
    L_lst = Parameters()[3]
    nl = np.size(L_lst)
    h = 1/L_lst
    T = np.zeros((nl+1,nl)) #table of extrapolants, with indices corresponding to m and L, respectively
    k = int(0.5*((nl-2)**2+(nl-2))) #number of epsilons
    epsilon = np.zeros(k)
    T[1] = tau #input finite-size data
    if q == 0: #avoid q = 0 limit
        tau_ext = -1
    elif q == 1: #avoid q = 1 limit
        tau_ext = 0
    else:
        for i in range(2,nl+1): #row of BST table
            for j in range(nl-i+1): #column of BST table
                T[i,j] = T[i-1,j+1]+(T[i-1,j+1]-T[i-1,j])/((h[j]/h[i+j-1])**omega*(1-(T[i-1,j+1]-T[i-1,j])/(T[i-1,j+1]-T[i-2,j+1]))-1)
        c = 0
        for i in range(2,nl):
            for j in range(nl-i):
                epsilon[c+j] = 2*(T[i,j+1]-T[i,j]) #differences between adjacent T elements
            c += j+1
        tau_ext = T[nl,0] #extrapolant taken as element in final row
    return tau_ext,epsilon

#scan across omega
def Scan(tau,q_p):
    L_lst = Parameters()[3]
    nl = np.size(L_lst)
    omega_min,omega_max,ng = Omega_Parameters()
    omega_lst = np.linspace(omega_min,omega_max,ng) #omega values
    tau_ext = np.zeros(ng)
    k = int(0.5*((nl-2)**2+(nl-2))) #number of epsilons
    epsilon = np.zeros((ng,k))
    for i in range(ng):
        tau_ext[i],epsilon[i] = Tau_ext(tau,omega_lst[i],q_p) #obtain extrapolated scaling exponents for all omega and q
    return tau_ext,epsilon

#find optimum omega by checking solution for each epsilon 
def Solution(tau_ext,epsilon,q_p):
    L_lst = Parameters()[3]
    ng = Omega_Parameters()[2]
    nl = np.size(L_lst)
    k = int(0.5*((nl-2)**2+(nl-2))) #number of epsilons
    ind = np.zeros(k)
    diff = np.zeros(k)
    for i in range(k):
        ind[i] = np.where(np.abs(epsilon[:,i])==np.min(np.abs(epsilon[:,i])))[0][0] #select solution for each epsilon
        index = int(ind[i])
        diff[i] = np.abs((q_p-1)-tau_ext[index]) #find deviation of each solution from q-1 line
    index = int(ind[np.where(diff==np.min(diff))[0][0]]) #select best solution
    sum = np.abs(np.sum(epsilon,axis=1))
    index_sum = np.where(sum==np.min(sum))[0][0] #solution for sum over epsilons
    if np.abs((q_p-1)-tau_ext[index]) > np.abs((q_p-1)-tau_ext[index_sum]): #use sum solution if optimal
        index = index_sum
    if tau_ext[index] < tau_ext[ng-1]: #criterion for MBL: if specific solutions are below large omega solution then assume MBL
        index = ng-1 #select large omega solution
    return index

#check where solutions are for omega for given T_0, T_1 and q, and plot against all epsilon
def Check(p_1,p_2):
    gamma,n_1,n_2,L_lst,q_lst,T_0_lst,T_1_lst,q_p = Parameters()
    omega_min,omega_max,ng = Omega_Parameters()
    omega_lst = np.linspace(omega_min,omega_max,ng)
    tau = Load(p_1,p_2) #load finite size data
    tau_ext = np.zeros(ng)
    q_ind = np.where(np.round(q_lst,1)==q_p)[0][0]
    tau_ext,epsilon = Scan(tau[:,q_ind],q_p) #obtain extrapolated taus and epsilons for every omega
    a = Solution(tau_ext,epsilon,q_p) #obtain optimum solution
    b = np.where(np.abs(q_lst[q_ind]-1-tau_ext)==np.min(np.abs(q_lst[q_ind]-1-tau_ext)))[0][0] #point where CUE line is crossed
    plt.axhline([0],linestyle='dashed',color='r',label='Loc.')
    plt.axhline([q_lst[q_ind]-1],linestyle='dashed',color='b',label='CUE')
    plt.plot(omega_lst,epsilon)
    plt.plot(omega_lst,np.sum(epsilon,axis=1),'g--',label=r'$\sum_{i,m}\epsilon_m^{(i)}$')
    plt.plot(omega_lst,tau_ext,'k',label=r'$\tau_q$')
    plt.annotate(r'solution = {}'.format(round(omega_lst[a],2)),(0.5,0.1),xycoords='axes fraction')
    plt.annotate(r'$\tau_q = {}$'.format(round(tau_ext[a],2)),(0.5,0.05),xycoords='axes fraction')
    plt.annotate(r'$\gamma T_0/2\pi = {}$'.format(np.round(gamma*T_0_lst[p_1]/(2*np.pi),1)),(0.8,0.1),xycoords='axes fraction')
    plt.annotate(r'$T_1 = {}$'.format(np.round(T_1_lst[p_2],2)),(0.8,0.05),xycoords='axes fraction')
    plt.annotate(r'deloc. $\omega = {}$'.format(round(omega_lst[b],2)),(0.2,0.1),xycoords='axes fraction')
    plt.xlabel(r'$\omega$')
    plt.ylabel(r'$\epsilon_m^{(i)}(\omega)$')
    plt.ylim([-3,3])
    plt.legend(loc='upper right',ncol=2,fontsize=8)
    plt.show()

#plot extrapolated scaling exponents for all q at given T_0 and T_1
def Tau_Plot(a,b):
    gamma,n_1,n_2,L_lst,q_lst,T_0_lst,T_1_lst,q_p = Parameters()
    nq = np.size(q_lst)
    tau = np.zeros((6,nq))
    tau_ext = np.zeros(nq)
    tau = Load(a,b) #load finite size data
    for i in range(nq):
        taus,epsilon = Scan(tau[:,i],q_lst[i]) #obtain extrapolated taus for every omega
        index = Solution(taus,epsilon,q_lst[i]) #obtain optimum solution
        tau_ext[i] = taus[index]
    plt.plot(q_lst,tau_ext,'-x')
    plt.plot(q_lst,q_lst-1,'k--')
    plt.axhline(y=0,linestyle='--',xmin=0.35,color='k')
    plt.annotate(r'$\gamma T_0/2\pi={}$'.format(np.round(gamma*T_0_lst[a]/(2*np.pi),2)),(0.05,0.86),xycoords='axes fraction')
    plt.annotate(r'$T_1={}$'.format(np.round(T_1_lst[b],2)),(0.05,0.81),xycoords='axes fraction')
    plt.show()

#plot phase diagram using the quantity 1-diff, which is the fractal dimension at given q
def Map():
    start = time.time()
    gamma,n_1,n_2,L_lst,q_lst,T_0_lst,T_1_lst,q_p = Parameters()
    nq = np.size(q_lst)
    tau = np.zeros((6,nq))
    tau_ext = np.zeros((n_1,n_2))
    diff = np.zeros((n_1,n_2))
    q_ind = np.where(np.round(q_lst,1)==q_p)[0][0] #find index of selected q
    for i in range(n_1):
        for j in range(n_2):
            tau = Load(i,j) #load finite size data
            taus,epsilon = Scan(tau[:,q_ind],q_p) #obtain extrapolated taus for every omega
            index = Solution(taus,epsilon,q_p) #obtain optimum solution
            tau_ext[i,j] = taus[index]
            diff[i,j] = (q_p-1)-(tau_ext[i,j]) #obtain deviation from q-1 line at q_p
    if q_p > 1:
        diff[np.where(diff<0)] = 0 #truncate negative differences
    else:
        diff[np.where(diff>0)] = 0 #truncate positive differences
    diff /= (q_p-1) #scale difference relative to q: 1 = loc, 0 = deloc
    end = time.time()
    print("Runtime:", (end-start), "seconds")
    Plot(1-diff)

def Plot(diff):
    gamma,n_1,n_2,L_lst,q_lst,T_0_lst,T_1_lst,q_p = Parameters()
    c = gamma/(2*np.pi)
    plt.figure(figsize=(7,7),dpi=80)
    plt.rcParams.update({'font.size':16})
    plt.imshow(diff.T,interpolation='gaussian',cmap='PuOr',origin='lower',extent=[c*T_0_lst[0],c*T_0_lst[n_1-1],c*T_1_lst[0],c*T_1_lst[n_2-1]],aspect='auto')
    plt.xlabel(r'$\gamma T_0/2\pi$')
    plt.ylabel(r'$\gamma T_1/2\pi$')
    plt.colorbar(fraction=0.0378)
    ax = plt.gca()
    ax.set_aspect(1.3)
    plt.clim(0,1)
    #plt.savefig('Stark_Files_2/Tau/BST/Tau_Map_2.0.pdf',bbox_inches='tight')
    plt.show()

def Parameters():
    gamma = 5 #tilt
    n_1 = 11 #T_0 increments
    n_2 = n_1 #T_1 increments
    L_lst = np.array([8,10,12,14,16]) #system size
    q_lst = np.linspace(0,3,16) #fractal dimension parameter
    n_lst = np.linspace(0.25,2.75,n_1)
    T_0_lst = 2*np.pi*n_lst/gamma #unkicked period
    T_1_lst = np.linspace(0.2,2.2,n_2) #kicked period
    q_p = 2.0 #selected q point
    return gamma,n_1,n_2,L_lst,q_lst,T_0_lst,T_1_lst,q_p

def Omega_Parameters():
    omega_min = 1.5 #minimum possible value for free parameter
    omega_max = 20 #maximum possible value for free parameter
    ng = 1000 #omega increments
    return omega_min,omega_max,ng

#Map(): generate full phase diagram
#Check(a,b): plot epsilon as a function of omega, with a,b corresponding to index of T_0 and T_1 data points respectively
#Tau_Plot(a,b): plot extrapolated tau_q as a function of q

Map()