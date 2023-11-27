import numpy as np
import math 
from matplotlib import pyplot as plt
import random
random.seed(1345677767)

ite = 100
lx = 16
ly = 16
nos = lx*ly
nb = 2 * nos
m = 20
nh = 0
beta = 20
energy = 0
s_heat = 0
msteps = 100


spin = np.zeros((nos+1),int)
fullconfig = np.zeros((nos+1,m*2),int)
opstring = np.zeros(m+1,int)
vertexlist = np.zeros(4*m+1,int)

bsites = np.zeros((2,nb+1),int)
for y1 in range(0,ly,1):
    for x1 in range(0,lx,1):
        s = 1+x1 +y1*lx
        x2 = (x1+1)%lx
        y2 = y1
        bsites[0,s] = s
        bsites[1,s] = 1+x2+y2*lx
        x2 = x1
        y2 = (y1+1)%ly
        bsites[0,s+nos] = s
        bsites[1,s+nos] = 1+x2+y2*lx
 
 
#-------------------------------------------------------------#
 
for i in range(nos):
    r = random.randint(0,10)
    if r<6:
        spin[i+1] = -1
    else:
        spin[i+1] = 1

def diagonalupdate(m,nh,opstring):
    for i in range(1,m+1,1):
    
        b = random.randint(1,nos) # choosing random site
        r = (random.randint(0,10000)/10000)
        s1,s2 = int(bsites[0,b]), int(bsites[1,b]) 
        sp1,sp2 = spin[s1],spin[s2]  # determining spins at that site and next to that site  
        addpb = (nb * (beta/2)) / (m-nh) #calculating probability for adding a diagonal operator 
        dropb = (m- nh + 1) / (nb * (beta/2))# for deleteing a diagonal operator 
        if opstring[i] == 0:  
            if r <= addpb:
                if sp1 != sp2:                         #adding a diagonal 
                    opstring[i] = 2*b 
                    nh =nh + 1
                    
        elif opstring[i]%2 ==0:
            if r<= dropb:
                opstring[i] = 0                        #deleteing a diagonal 
                nh = nh-1
                
        else:                                             #switching the spins after encountering off diagonal operators 
            b = int(opstring[i]/2)
            
            spin[s1],spin[s2] = -sp1,-sp2
    return m,nh,opstring

def vertexlinker(m,opstring,vertexlist):
    vertexlist[:] = 0
    fullconfig[:,:] = 0 
    for i in range(1,m+1,1):
        p = opstring[i]
        vert = int(p/2)
        s1,s2 = bsites[0,vert],bsites[1,vert]

        if p !=0:
            fullconfig[s1,2*i-2] = 1
            fullconfig[s2,2*i-2] = 2
            fullconfig[s2,2*i-1] = 4
            fullconfig[s1,2*i-1] = 3
    

    
    
    
    for i in range(0,2*m,2):
        op = opstring[int(i/2) +1]
    
        found = np.zeros(4,int)
        if op !=0:
            p = int(i/2) +1
            vert = int(op/2)
            s1,s2 = bsites[0,vert],bsites[1,vert]
            v1 = p*4-3
            v2 = p*4-2
            v3 = p*4-1
            v4 = p*4
            for j in range(i,0,-1):      # searching above the operator upto p=1 for the connection of leg =0 of pth op
                if fullconfig[s1,j] > 2:
                    if fullconfig[s1,j] == 4: # if 0th leg connects to l=3 of some op
                        pfound = int(j/2) + 1
                        vertexlist[v1] = pfound*4
                        vertexlist[pfound*4] = v1 
                        found[0] =1
                        break
                
                    if fullconfig[s1,j] == 3: #if 0th leg connects to l=2 of some op
                        pfound = int(j/2) + 1
                        vertexlist[v1] = pfound*4-1
                        vertexlist[pfound*4-1] = v1 
                        found[0] = 1
                        break
            
            for j in range(i,0,-1):   # searching above the operator upto p=1 for the connection of leg =1 of pth op
                if fullconfig[s2,j] > 2:
                    if fullconfig[s2,j] == 4: # if 0th leg connects to l=3 of some op
                        pfound = int(j/2) + 1
                        vertexlist[v2] = pfound*4
                        vertexlist[pfound*4] = v2 
                        found[1] =1
                        break
                    
                    if fullconfig[s2,j] == 3: #if 0th leg connects to l=2 of some op
                        pfound = int(j/2) + 1
                        vertexlist[v2] = pfound*4-1
                        vertexlist[pfound*4-1] = v2 
                        found[1] = 1
                        break

            for j in range(2*m-1,0,-1): #starting the search from p=m to p=0 
                if found[0]==1 and found[1] == 1: #breaking if found connection for both leg
                    break
                if found[0] == 0: #finding for 0th leg
                    if fullconfig[s1,j] == 4:
                        pfound = int(j/2) + 1
                        vertexlist[v1] = pfound*4
                        vertexlist[pfound*4] = v1 
                        found[0] =1
                
                    if fullconfig[s1,j] == 3:
                        pfound = int(j/2) + 1
                        vertexlist[v1] = pfound*4-1
                        vertexlist[pfound*4-1] = v1 
                        found[0] = 1
                if found[1] == 0: #finding for 1st leg
                    if fullconfig[s2,j] == 4:
                        pfound = int(j/2) + 1
                        vertexlist[v2] = pfound*4
                        vertexlist[pfound*4] = v2
                        found[1] =1
                
                    if fullconfig[s2,j] == 3: #finding for 2nd leg
                        pfound = int(j/2) + 1
                        vertexlist[v2] = pfound*4-1
                        vertexlist[pfound*4-1] = v2 
                        found[1] = 1
   #---------------------------------------------------------------------------     
    #starting the search for connection for lower legs(leg 2 and 3)
            for j in range(i+1,2*m,1): 
                if fullconfig[s1,j] < 3:
                    if fullconfig[s1,j] == 2:
                        pfound = int(j/2) + 1
                        vertexlist[v3] = pfound*4-2
                        vertexlist[pfound*4-2] = v3 
                        found[3] =1
                        break
                
                    if fullconfig[s1,j] == 1:
                        pfound = int(j/2) + 1
                        vertexlist[v3] = pfound*4-3
                        vertexlist[pfound*4-3] = v3
                        found[2] = 1
                        break
            
            for j in range(i,2*m,-1):
                if fullconfig[s2,j] < 3:
                    if fullconfig[s2,j] == 2:
                        pfound = int(j/2) + 1
                        vertexlist[v4] = pfound*4-2
                        vertexlist[pfound*4-2] = v4 
                        found[3] =1
                        break
                
                    if fullconfig[s2,j] == 1:
                        pfound = int(j/2) + 1
                        vertexlist[v4] = pfound*4-3
                        vertexlist[pfound*4-3] = v4 
                        found[2] = 1
                        break
    
    return vertexlist

def adjv(i):
    if i%2 == 0:
        v1= i
        v2 = v1^1 
        v2 = v2 -2
        return v2
    else :
        v1 = i
        v2 = v1^1
        v2 = v2 +2
        return v2

def loopupdate(m,vertexlist,opstring,spin):
    global nos,F
    for i in range(1,4*m+1,2):
        r = random.randint(0,1000)/1000
        if vertexlist[i]>0:
            v1 = i
            if r>0:
                for j in range(1,100):
                    if opstring[int((v1-1)/4)+1]%2 ==0:
                        opstring[int((v1-1)/4)+1] = opstring[int((v1-1)/4)+1]+1
                        
                    else:
                        opstring[int((v1-1)/4)+1] = opstring[int((v1-1)/4)+1] -1
                        
                    vertexlist[v1] = -1
                    v2 = adjv(v1)
                    v1 = vertexlist[v2]
                    vertexlist[v2] = -1
                    if v1 == i:
                        break
            else:
                for j in range(1000):
                    vertexlist[v1] = 0
                    v2 = adjv(v1)
                    v1 = vertexlist[v2]
                    vertexlist[v2] = 0
                    if v1 ==i:
                        break
    for i in range(1,4*m+1,1):
        if vertexlist[i] == -1:
            p = int((i-1)/4) + 1
            op = opstring[p]
            b = int(op/2)
            s1,s2 = int(bsites[0,b]), int(bsites[1,b])
            if i%4 == 1 or i%4 == 2:
                spin[s1] = -spin[s1]
            else:
                spin[s2] = -spin[s2]
            
    return opstring,vertexlist,spin

def adjustcutoff(m,nh):
    #for adjusting cutoff
    mnew = nh + int(nh/2)
    print(nh)
    if mnew > m: 
        temp = np.zeros(m+1,int)
        temp = opstring
        opstring = np.zeros(mnew+1,int)
        opstring[0:m+1] = temp
        mnew = m
        fullconfig = np.zeros((nos+1,m*2),int)
        vertexlist = np.zeros(4*m+1,int)
        

    return m

# MAIN FUNCTION
mdata = np.zeros(ite,int)
nhdata = np.zeros(ite,int)
p_ad = np.zeros(ite)
p_del = np.zeros(ite)
for k in range(ite):
    
    m,nh,opstring = diagonalupdate(m,nh,opstring)
    vertexlist = vertexlinker(m,opstring,vertexlist)
    opstring,vertexlist,spin = loopupdate(m,vertexlist,opstring,spin)
    mnew = nh + int(nh/2)
    p_ad[k] = (nb * (beta/2)) / (m-nh)
    p_del[k] = (m- nh + 1) / (nb * (beta/2))
    if mnew > m: 
        temp = np.zeros(m+1,int)
        temp = opstring
        opstring = np.zeros(mnew+1,int)
        opstring[0:m+1] = temp
        m = mnew
        fullconfig = np.zeros((nos+1,m*2),int)
        vertexlist = np.zeros(4*m+1,int)
    
    mdata[k] = m
    nhdata[k] =nh
       
print("nh:",nh,"m:",m)

x = np.linspace(0,ite,ite)
plt.plot(x,mdata)
plt.show()

plt.plot(x,nhdata)
plt.show()

plt.plot(x,p_ad)
plt.show()

plt.plot(x,p_del)
plt.show()