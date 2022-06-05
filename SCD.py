import numpy as np

#Computation of the multiplicities in the direct sum decomposition of the tensor product between a a direct sum of spins and a spin s
def tensor(a, s):
    new_a = np.empty((0,2))
    S = a[0,1] + s
    while S >= 0:
        mult = 0
        for spin in a:
            if S >= abs(spin[1] - s) and S <= spin[1] + s:
                mult += spin[0]
        new_a = np.vstack((new_a, [mult, S]))
        S += -1
    
    return new_a

#Decomposition of a product of spins
def decompose(spins):
    a = np.array([[1,spins[0]]])
    for spin in spins[1:]:
        a = tensor(a, spin)
    return a

#Analyse decomposition of product of spins
def reduction(spins):
    reduced = decompose(spins)
    tot_hilbert = 1
    print(reduced)
    for spin in spins:
        tot_hilbert = tot_hilbert * (2*spin + 1)
    print(tot_hilbert)
    print('the matrix is recuduced into blocks of the following size')
    print('The biggest block has size', int(reduced.max()))
    print('this is a reduction of', (1 - reduced.max()/tot_hilbert)*100, 'percent with respect to the original')
    tot_mat = 1
    for spin in reduced:
        tot_mat += spin[0]**2
    print('There are', int(tot_mat), 'matrix elements to be determined')
    print('This a reduction of', (1-tot_mat/(float(tot_hilbert)**2))*100, 'percent')