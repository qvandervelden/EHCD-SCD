"""
HC-Circle: Heisenberg Chain Circle
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
import numpy.linalg as lin
import time
import EHCD
import SCD as dc
import os
import matplotlib.pyplot as plt


if __name__ ==  '__main__':
    
    start = time.time() #Starting the timer for the runtime
    
    S = 0.5 #The spin of the particles
    N = 14 #The number of particles
    J_s = -1 #The coupling between neighbouring particles
    
    spins = EHCD.unispin(S, N) #The set of spins
    J = EHCD.unicircle(N, J_s) #The coupling matrix
    
    print(J)
    
    coupling = EHCD.couplings(spins) #Determination of the coupling sequences
    
    E = np.empty((0,2)) #Array for eigenenergies
    
    #Determination of highest spin with a nontrivial spinsector
    total_s = dc.decompose(spins)[:,1]
    
    #Determination of the eigenenergies for every spin-sector
    for s in total_s:
        print(s, time.time() - start)
        H_s = EHCD.block(spins, J, coupling, s)
        E_s = lin.eigvalsh(H_s)
        E_s = np.vstack((E_s, np.full((1, np.shape(E_s)[0]), s)))
        E = np.vstack((E, np.transpose(E_s)))
    
    #plotting spectrum
    plt.figure()
    plt.scatter(E[:,1],E[:,0], marker='_', linewidths=0.5)
    plt.show()
        
    #Saving spectrum
    if not os.path.exists('./Results/Spectrum-' + str(N) + 'x' + str(S) + 'c.txt'):
        np.savetxt('./Results/Spectrum-' + str(N) + 'x' + str(S) + 'c.txt', E, header = 'Runtime: ' + str(time.time() - start))
    else:
        i = 1
        while True:
            if not os.path.exists('./Results/Spectrum-' + str(N) + 'x' + str(S) + '(' + str(i) + ')' + 'c.txt'):
                np.savetxt('./Results/Spectrum-' + str(N) + 'x' + str(S) + '(' + str(i) + ')' + 'c.txt', E, header = 'Runtime: ' + str(time.time() - start))
                break
            i += 1
            
    print(time.time() - start) #Printing runtime