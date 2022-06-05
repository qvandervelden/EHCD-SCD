import numpy as np
import numpy.linalg as lin
import time
import EHCD
import SCD as dc
import os
import matplotlib.pyplot as plt


if __name__ ==  '__main__':
    
    filename = 'spectrum'
    
    start = time.time()
    
    S = 1
    N = 9
    J_s = -1
    
    spins = np.array([1,1,1,1,1,1,1,1,1])
    J = np.array([[0,1,0,1,0,0,0,0,0],
                  [1,0,1,0,1,0,0,0,0],
                  [0,1,0,1,0,1,0,0,0],
                  [1,0,0,0,1,0,1,0,0],
                  [0,1,0,1,0,1,0,1,0],
                  [0,0,1,0,1,0,0,0,1],
                  [0,0,0,1,0,0,0,1,0],
                  [0,0,0,0,1,0,1,0,1],
                  [0,0,0,0,0,1,0,1,0]])
    
    print(J)
    
    coupling = EHCD.couplings(spins)
    
    E = np.empty((0,2))
    
    total_s = dc.decompose(spins)[:,1]
    
    for s in total_s:
        print(s, time.time() - start)
        H_s = EHCD.block(spins, J, coupling, s)
        E_s = lin.eigvalsh(H_s)
        E_s = np.vstack((E_s, np.full((1, np.shape(E_s)[0]), s)))
        E = np.vstack((E, np.transpose(E_s))) 
    
    #Plotting spectrum
    plt.figure()
    plt.scatter(E[:,1],E[:,0], marker='_', linewidths=0.5)
    plt.show()
        
    
    if not os.path.exists('./Results/' + filename + '.txt'):
        np.savetxt('./Results/' + filename + '.txt', E, header = 'Runtime: ' + str(time.time() - start))
    else:
        i = 1
        while True:
            if not os.path.exists('./Results/' + filename + str(S) + '(' + str(i) + ')' + '.txt'):
                np.savetxt('./Results/' + filename + '(' + str(i) + ')' + '.txt', E, header = 'Runtime: ' + str(time.time() - start))
                break
            i += 1
            
    print(time.time() - start)