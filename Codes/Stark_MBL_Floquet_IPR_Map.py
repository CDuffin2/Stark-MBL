import numpy as np
import scipy
import scipy.sparse as spr
import scipy.sparse.linalg as splng
from numba import jit
import matplotlib.pyplot as plt
import time

#IPRs averaged over evecs and disorder calculated as a function
#of T_0 and T_1

#calculate dimensionality for half filling
@jit
def Binom(n,k):
    ret = 1
    for i in range(min(k,n-k)):
        ret *= n-i
        ret /= i+1
    return ret

#state to label
@jit
def S2L(L,l_vec):
    N = L/2
    l_counter = 1
    n_counter = N
    for i in range(L):
        if l_vec[i] == 0:
            l_counter += Binom(L-(i+1),n_counter-1)
        if l_vec[i] == 1:
            n_counter -= 1
        if n_counter == 0:
            return l_counter

#label to state
@jit
def L2S(L,l):
    N = L/2
    n_counter = N
    l_counter = l
    l_vec = np.zeros(L)
    for i in range(L):
        binom = Binom(L-(i+1),n_counter-1)
        if n_counter == 0:
            l_vec[i] = 0
            continue
        if l_counter > binom:
            l_vec[i] = 0
            l_counter -= binom
            continue
        if l_counter <= binom:
            l_vec[i] = 1
            n_counter -= 1
    return l_vec

#generate states via hopping terms
def Hop(L,input):
    input = np.ndarray.astype(input,int) #input state
    output_lst = [] #appended list of output states
    for i in range(L-1):
        if input[i] == 1 and input[i+1] == 0: #hopping to the right
            output = np.copy(input)
            output[i] = 0
            output[i+1] = 1
            output_lst.append(output)
        if input[i] == 0 and input[i+1] == 1: #hopping to the left
            output = np.copy(input)
            output[i] = 1
            output[i+1] = 0
            output_lst.append(output)
    ns = np.shape(output_lst)[0] #number of output states
    return output_lst,ns

#diagonal elements of static Hamiltonian
def Diag(L,Dim,gamma,w,V):
    row = np.empty(Dim*L)
    col = np.empty(Dim*L)
    data = np.empty(Dim*L)
    nnz = 0
    for i in range(Dim):
        state = L2S(L,i+1) #select vector
        v = 0 #nearest-neighbour interaction counter
        g = 0 #field interaction counter
        h = 0 #disorder interaction counter
        for k in range(L):
            g = g-gamma*k*state[k] #count field interactions
            h = h+np.random.normal(0,w)*state[k] #count disorder interactions
        for k in range(L-1):
            v = v+V*state[k]*state[k+1] #count nearest-neighbour interactions
        row[nnz] = i #H[i,i] = g+h+v
        col[nnz] = i
        data[nnz] = g+h+v
        nnz += 1
    return spr.csc_matrix((data[:nnz],(row[:nnz],col[:nnz])),shape=(Dim,Dim))

#off-diagonal elements of static Hamiltonian
def Offdiag(L,Dim,J):
    row = np.empty(Dim*L)
    col = np.empty(Dim*L)
    data = np.empty(Dim*L)
    nnz = 0
    for i in range(Dim):
        input = L2S(L,i+1) #generate Fock state
        output_lst,ns = Hop(L,input) #generate new states via hopping
        for k in range(ns): #loop through all states generated by hopping
            output = output_lst[k] #select output state
            j = int(S2L(L,output)-1) #convert to label
            row[nnz] = i #H[i,j] = J
            col[nnz] = j
            data[nnz] = J
            nnz += 1
    return spr.csc_matrix((data[:nnz],(row[:nnz],col[:nnz])),shape=(Dim,Dim))

#calculate IPRs
def IPR(T_0,T_1,H_0,H_1):
    A = np.diag(np.exp(np.diag(-1j*H_0.todense()*T_0))) #exponentiate diagonal matrix
    B = scipy.linalg.expm(-1j*H_1.todense()*T_1) #exponentiate off-diagonal matrix
    F = A@B #Floquet operator
    U = np.linalg.eig(F)[1] #calculate evecs
    U = np.array(U) #convert from matrix object to array
    I = np.mean(np.sum(np.abs(U)**4,axis=1)) #average IPR over evecs
    return np.mean(I) #average IPR over disorder

#loop through system sizes and unkicked periods
def Main():
    L,J,V,gamma,w,n_1,n_2,T_0_lst,T_1_lst,n_lst = Parameters() #system parameters
    I = np.zeros((n_1,n_2)) #size of IPR array
    Dim = int(Binom(L,L/2)) #dimensionality
    for i in range(n_1):
        T_0 = T_0_lst[i] #update T_0
        H_1 = Offdiag(L,Dim,J) #off-diagonal part of Hamiltonian
        for j in range(n_2):
            T_1 = T_1_lst[j] #update T_1
            H_0 = Diag(L,Dim,gamma,w,V) #diagonal part of Hamiltonian
            I[i,j] = IPR(T_0,T_1,H_0,H_1) #calculate IPR
    return I

#plot IPRs
def Plot(I):
    L,J,V,gamma,w,n_1,n_2,T_0_lst,T_1_lst,n_lst = Parameters()
    plt.imshow(I.T,interpolation='gaussian',cmap='inferno',origin='lower',extent=[n_lst[0],n_lst[n_1-1],T_1_lst[0],T_1_lst[n_2-1]],aspect='auto')
    plt.xlabel(r'$\gamma T_0/2\pi$')
    plt.ylabel(r'$T_1$')
    plt.colorbar()
    plt.title(r'$L={}$'.format(L))
    plt.show()

#system parameters
def Parameters():
    L = 14 #system size
    J = 0.5 #hopping amplitude
    V = 0.5 #interaction strength
    gamma = 5 #potential strength
    w = 0.05 #disorder strength
    n_1 = 29 #T_0 iterations
    n_2 = n_1 #T_1 iterations
    n_lst = np.linspace(0.1,2.9,n_1)
    T_0_lst = 2*np.pi*n_lst/gamma #unkicked period
    T_1_lst = np.linspace(0,1,n_2) #kicked period
    return L,J,V,gamma,w,n_1,n_2,T_0_lst,T_1_lst,n_lst

start = time.time()
I = Main()
end = time.time()
print("Runtime:", (end-start)/60, "minutes")
Plot(I)
