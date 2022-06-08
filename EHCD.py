"""
EHCD: Efficient Heisenberg Chain Decomposition
Copyright (C) 2022 Quinn van der Velden

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import numpy as np
import sympy as sym
import sympy.physics.wigner as wig
import concurrent.futures as fut
import time

#determination of all ways of coupling
def couplings(spins):
    
    #determination of all ways of coupling for two particles
    if spins.shape[0] == 2:
        s = spins[0] + spins[1]
        stop = abs(spins[0] - spins[1])
        coupling = np.empty((0,3))
        while s >= stop:
            coupling = np.vstack((coupling, np.array([spins[0],spins[1], s])))
            s += -1
        return coupling
    
    #determination of first way of coupling for more than two particles
    JM = np.array([spins[0],spins[1], spins[0] + spins[1]])
    for i in range(2, spins.size):
        JM = np.append(JM, [spins[i], JM[-1] + spins[i]])
        
    coupling = np.array([JM]) #array with ways of coupling
    
    #generating all ways of coupling
    while True:
        JM_old = coupling[-1]
        j = JM_old.shape[0] - 1
        while True:
            if JM_old[j] > abs(JM_old[j-2] - JM_old[j - 1]):
                JM_new = JM_old[0:j]
                break
            else:
                j = j - 2
        if j < 2:
            break
        JM_new = np.append(JM_new, [JM_old[j]-1])
        while j < JM_old.size - 1:
            JM_new = np.append(JM_new, [spins[int(j/2 + 1)], JM_new[-1] + spins[int(j/2 + 1)]])
            j = j + 2
        coupling = np.vstack((coupling, JM_new))
        if JM_old.shape[0] == 3:
            break
    
    #sorting of ways of coupling by resulting spin
    coupling = coupling[(-1*coupling[:, JM.size - 1]).argsort()]
    return coupling

#determination of coupling of tensors, make sure i < j
def incoup(i,j, N):
    if i == 0 and j == 1:
        k = [1,1,0]
    elif i == 1:
        k = [0,1,1]
    elif i == 0:
        k = [1,0,1]
    else:
        k = [0,0,0]
    for l in range(2, N):
        if i == l or j == l:
            if k[-1] == 0:
                k = np.append(k, [1,1])
            else:
                k = np.append(k, [1,0])
        else:
            if k[-1] == 0:
                k = np.append(k, [0,0])
            else:
                k = np.append(k, [0,1])
    return k

#determination of reduced matrix element of spin STO of rank 0
def red_s_0(s):
    return np.sqrt(2*s+1)

#determination of reduced matrix element of spin STO of rank 1
def red_s_1(s):
    return np.sqrt(s*(s+1)*(2*s+1))

#optimized determination of the wigner-9j symbol for the case j_3,j_6,j_9=0
def wig_9000(a,b,c,d,e,f):
    if a == b and c == d and e == f and a + c >= e >= np.abs(a - c):
        return (-1)**(a + b + 2*c + 2*e)/(np.sqrt((2*a+1)*(2*c + 1)*(2*e+1)))
    else:
        return 0

#optimized determination of the wigner-9j symbol for the case j_3,j_6 = 1 ,j_9=0
def wig_9110(a,b,c,d,e,f):
    if e == f:
        return (-1)**(1 + b + c + e)*sym.N(wig.wigner_6j(a, b, 1, d, c, e)/(np.sqrt(3*(2*e+1))))
    else:
        return 0
    
#optimized determination of the wigner-9j symbol for the case j_3=0 ,j_6,j_9=1
def wig_9011(a,b,c,d,e,f):
    if a == b:
        return (-1)**(1 + 2*a + b + 2*c + d + e + 2*f)*sym.N(wig.wigner_6j(e, f, 1, d, c, a)/(np.sqrt(3*(2*a+1))))
    else:
        return 0

#optimized determination of the wigner-9j symbol for the case j_3=1 ,j_6=0 ,j_9=1
def wig_9101(a,b,c,d,e,f):
    if c == d:
        return (-1)**(a + 2*b + 2*c + d + f + 2*e + 1)*sym.N(wig.wigner_6j(a, b, 1, f, e, c)/(np.sqrt(3*(2*c+1))))
    else:
        return 0

#determination of the matrixblock for some term
def mat_t(J, i, j, spins, coup):
    N = spins.shape[0]
    d = coup.shape[0] #dimension of the block
    coup_T = incoup(i,j,N) #determination of the coupling of the STOs
    
    
    #working out the prefactor of the block
    red_base = -1*np.sqrt(3/(2*coup[0,-1] + 1))*J
    for l in range(0, N):
        if l != i and  l != j:
            red_base = red_base*red_s_0(spins[l])
        else:
            red_base = red_base*red_s_1(spins[l])
    T = np.full((d,d), red_base)
    
    #Numbering of coupling sequences
    index = np.array([[0]])
    for l in range(1, d):
        index = np.vstack((index, [[l]]))
    coup = np.hstack((coup, index))
    
    #Computation of term
    for l in range(0, N - 1):
        
        #Permutation of coupling sequences to create blocks that recieve same coupling factor
        ind = np.lexsort((coup[:, 2*l+2], coup[:, 2*l+1], coup[:, 2*l]))
        coup = coup[ind]
        T = T[ind, :]
        T = T[:, ind]
        
        #Determining blocks
        c = (coup[0][2*l] + 1, 0, 0)
        blocks = []
        for p in range(0, d):
            if c != (coup[p][2*l], coup[p][2*l+1], coup[p][2*l+2]):
                blocks.append(p)
                c = (coup[p][2*l], coup[p][2*l+1], coup[p][2*l+2])
        blocks.append(d)
        
        #Determination and application of coupling factors for upper diagonal
        for q in range(0, len(blocks) - 1):
            for p in range(0, q + 1):
                
                if coup_T[2*l] == 0 and coup_T[2*l+1] == 0 and coup_T[2*l+2] == 0:
                    frac = np.sqrt((2*coup[blocks[p]][2*l+2]+1)*(2*coup[blocks[q]][2*l+2]+1)*(2*coup_T[2*l+2]+1))*wig_9000(coup[blocks[p]][2*l], coup[blocks[q]][2*l], coup[blocks[p]][2*l+1], coup[blocks[q]][2*l+1], coup[blocks[p]][2*l+2], coup[blocks[q]][2*l+2])
                elif coup_T[2*l] == 1 and coup_T[2*l+1] == 1 and coup_T[2*l+2] == 0:
                    frac = np.sqrt((2*coup[blocks[p]][2*l+2]+1)*(2*coup[blocks[q]][2*l+2]+1)*(2*coup_T[2*l+2]+1))*wig_9110(coup[blocks[p]][2*l], coup[blocks[q]][2*l], coup[blocks[p]][2*l+1], coup[blocks[q]][2*l+1], coup[blocks[p]][2*l+2], coup[blocks[q]][2*l+2])
                elif coup_T[2*l] == 0 and coup_T[2*l+1] == 1 and coup_T[2*l+2] == 1:
                    frac = np.sqrt((2*coup[blocks[p]][2*l+2]+1)*(2*coup[blocks[q]][2*l+2]+1)*(2*coup_T[2*l+2]+1))*wig_9011(coup[blocks[p]][2*l], coup[blocks[q]][2*l], coup[blocks[p]][2*l+1], coup[blocks[q]][2*l+1], coup[blocks[p]][2*l+2], coup[blocks[q]][2*l+2])
                else:
                        frac = np.sqrt((2*coup[blocks[p]][2*l+2]+1)*(2*coup[blocks[q]][2*l+2]+1)*(2*coup_T[2*l+2]+1))*wig_9101(coup[blocks[p]][2*l], coup[blocks[q]][2*l], coup[blocks[p]][2*l+1], coup[blocks[q]][2*l+1], coup[blocks[p]][2*l+2], coup[blocks[q]][2*l+2])
                T[blocks[p]:blocks[p+1],blocks[q]:blocks[q+1]] = T[blocks[p]:blocks[p+1],blocks[q]:blocks[q+1]] * frac
                
        #mirroring upper to lower diagonal
        for q in range(1,d):
            for p in range(0, q):
                T[q,p] = T[p,q]
        
        #Undoing permutation of coupling sequences
        ind = coup[:, -1].argsort()
        coup = coup[ind]
        T = T[ind,:]
        T = T[:,ind]
        
    return T



#computes the block corresponding to spin S
def block(spins, J, coupling, S):
    prog_coup = coupling[coupling[:,-1] == S]
    d = prog_coup.shape[0] #dimension of block
    H_S = np.zeros((d,d)) #block
    
    #Computation of terms, terms are distributed over cpu cores
    results = []
    runtime = time.time()
    count = 0
    max_count = 0
    with fut.ProcessPoolExecutor() as executor:
        for j in range(0,J.shape[0]):
            for i in range(0, j+1):
                if J[i,j] != 0:
                    results.append(executor.submit(mat_t, J[i,j], i,j, spins, prog_coup))
                    max_count += 1
        for T in fut.as_completed(results):
            H_S = H_S + T.result()
            count += 1
            if time.time() > runtime + 60:
                print(count, '/', max_count)
                runtime = time.time()
    return H_S

def uniline(N, J_s):
    J = np.zeros((N,N))
    for i in range(0,N):
        for j in range(0,N):
            if np.abs(i - j) == 1:
                J[i,j] = J_s
    return J

#returns coupling matrix for the Heisenberg Circle of N particles with coupling J_s
def unicircle(N, J_s):
    J = np.zeros((N,N))
    for i in range(0,N):
        for j in range(0,N):
            if np.abs(i - j) == 1:
                J[i,j] = J_s
    J[0,-1] = -1
    J[-1,0] = -1
    return J