# Import libraries
import numpy as np
import random
import math
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------

# Input Parameters
Nx = 16
Ny = 8
L = 40
MaxL = 9000
beta = 2
alpha = 3
Nsites = Nx*Ny
Nbonds = Nsites-1
max_iter = 1000
iter = 0
n = 0
sp = np.zeros(MaxL,np.int16)
a = np.zeros(MaxL,np.int16)
b = np.zeros(MaxL,np.int16)
spins = [(random.randint(1,Nsites)-0.5)/abs(random.randint(1,Nsites)-0.5) for i in range(Nsites)]

initial_magnetization = np.sum(spins)

# ------------------------------------------------------------------------------

# Generate a list of labels for each bond in a 2D lattice
def generate_bonds(Nx,Ny):
    dim = Nx*Ny
    bonds = np.zeros([2,2*dim],np.int32)
    for i in range(Ny):
        for j in range(Nx):
            s = j+i*Nx
            x = (j+1)%Nx 
            y = i
            s_p = x+y*Nx
            bonds[0,s] = s
            bonds[1,s] = s_p
            x = j
            y = (i+1)%Ny
            s_p = x+y*Nx
            bonds[0,s+dim] = s
            bonds[1,s+dim] = s_p
    return bonds  

bonds = generate_bonds(Nx,Ny)

"""-----Diagonal update-----"""

# Diagonal update addition probability  
def prob_add(L,Nb,beta,n):
    if n == L:
        return 1
    else:
        return min((beta*Nb)/(2*(L-n)),1)

# Diagonal update removal probability
def prob_remove(L,Nb,beta,n):
    return min(((L-n+1))/((beta/2)*Nb),1)

def diag_update(L,Nb,beta,n,spin,a,b):
    for p in range(0,L):
        bnd = random.randint(1,2*Nb) 
        # Adjacent spins in the bond labelled by the number bnd
        #print(bonds[0,bnd],bonds[1,bnd])
        spin_i = spin[bonds[0,bnd]]
        spin_j = spin[bonds[1,bnd]]
        if a[p] == 0:
            if spin_i == spin_j:
                continue
            # Metropolis update
            if random.random() < prob_add(L,Nb,beta,n):
                # accept the operator on the bond 
                #print('added')
                b[p] = bnd
                a[p] = 1
                sp[p] = 2*b[p]
                n += 1
        elif a[p] == 1:
            # Metropolis update
            if random.random() < prob_remove(L,Nb,beta,n):
                # remove the operator on the bond
                b[p] = 0
                a[p] = 0
                sp[p] = 0
                n -= 1
        else:
            # Switch the spins
            spin[bonds[0,b[p]]] *= -1
            spin[bonds[1,b[p]]] *= -1
    return n

"""-----Generate a vertix list-----"""

def build_vertex_list(sp,Nsites,L):
    xv = -1*np.ones(L*4,np.int32)
    v_first = -1*np.ones(Nsites,np.int32)
    v_last = -1*np.ones(Nsites,np.int32)
    for i in range(L):
        if sp[i] == 0:
            continue
        v0 = 4*i 
        bnd = int(sp[i]/2)
        s1 = bonds[0,bnd]
        s2 = bonds[1,bnd]
        v1 = v_last[s1]
        v2 = v_last[s2]

        if v1 != -1:
            xv[v1] = v0
            xv[v0] = v1     
        else:
            v_first[s1] = v0
        if v2 != -1:
            xv[v2] = v0+1
            xv[v0+1] = v2              
        else:
            v_first[s2] = v0+1
        v_last[s1] = v0+2
        v_last[s2] = v0+3
        
    for i in range(Nsites):
        f = v_first[i]
        if f != -1:
            l = v_last[i]
            xv[f] = l
            xv[l] = f

    return v_first,v_last,xv

"""-----Loop Update-----"""

# Flip the bit of the number
def flip_bit(number, i):
    mask = 1 << i  # Create a mask with the i-th bit set to 1
    flipped_number = number ^ mask  # XOR the number with the mask to flip the i-th bit
    return flipped_number

# Traverses a vertex list and marks
def traverse_loop(v0,xv,flip,sp):
    v = v0
    while True:
        if flip:
            xv[v] = -2
            p = int(v/4)
            sp[p] = flip_bit(sp[p],0)
        else:
            xv[v] = -1
        v_prime = flip_bit(v,0)
        v = xv[v_prime]
        if flip:
            xv[v_prime] = -2
        else:
            xv[v_prime] = -1
        if v == v0:
            break

# Updates the operator types in b according to loops given by xv
def loop_update(L,xv,a,sp,vfirst,Nsites,spins):
    for v0 in range(0,4*L,2):
        if xv[v0] < 0:
            continue
        if random.uniform(0,1) < 1/2:
            flip = False
            traverse_loop(v0,xv,flip,sp)
        else:
            flip = True
            traverse_loop(v0,xv,flip,sp)
    for i in range(Nsites):
        v = vfirst[i]
        if v == -1:
            if random.uniform(0,1) < 1/2:
                spins[i] *= -1
        else:
            v += 1
            if xv[v] == -2:
                spins[i] *= -1
    for i in range(len(sp)):
        if sp[i] > 0:
            a[i] = (sp[i]%2)+1
     
counter_n = []
counter_L = []
counter_iter =[]
while iter < max_iter:
    
    n = diag_update(L,Nbonds,beta,n,spins,a,b)
    if n<=0:
        continue
    
    v_first,v_last,xv = build_vertex_list(sp,Nsites,L)
    loop_update(L,xv,a,sp,v_first,Nsites,spins)
    
    if L-n < n/alpha:
        L = math.trunc(n+(n/alpha))
    counter_n.append(n)
    counter_L.append(L)
    
    iter += 1
    counter_iter.append(iter)
    
plt.plot(counter_iter,counter_n,label='n')
plt.plot(counter_iter,counter_L,label='L')
plt.legend()
plt.show()

""" ------- Compute the observables of the system ------- """

# Compute the energy of the system
avg_n = np.sum(counter_n)/len(counter_n)
energy = -avg_n/beta

# Compute the heat capacity of the system
avg_n_sq = np.sum(np.multiply(counter_n,counter_n))/len(counter_n)
heat_capacity = avg_n_sq - avg_n**2 - avg_n

# Compute the magnetization of the system
final_magnetization = np.sum(spins)
print(np.multiply([1,2],[1,3]))
print('The total energy of the system is:', energy)
print('The total heat capacity of the system is:', heat_capacity/Nsites)
print('The total initial and final magnetization is:', initial_magnetization, final_magnetization)