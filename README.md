# EHCD (Efficient Heisenberg Chain Diagonalization) + SCD (Spin Coupling Decomposition)
Software for determining the spectra of Heisenberg chains

This folder contains the following file

  * EHCD.py: Python file containing functions to determine the blocks of the matrix of a Heisenberg chain disgonalized with respect to total spin.
  * SCD.py: Python file containing functions to study the decomposition of a system of coupled spins
  * HC-Circle.py: Python file for determining the spectra of the Heisenberg circle: a Heisenberg chain with particles in a line, coupled to their neighbours and the particles at both ends also coupled, all with the same coupling strenght.
  * HC-General.py: Python file for determining the spectra of a general Heisenberg chain.
  * HC-Line-bipartite.py: Heisenberg chain with particles in a line, coupled to their neighbours. The particles can be divided in two sets A en B with particles in the same set being coupled ferromagnetically with strength J_s and particles of different sets coupled antiferromagnetically with strength -J_s
  * HC-Line.py: Python file for determining the spectra of the Heisenberg line: a Heisenberg chain with particles in a line, coupled to their neighbours, all with the same coupling strenght.


All the files above above were written by Quinn van der Velden (q.vandervelden@students.uu.nl). The code is made available under the GPLv3 License. See the
file LICENSE for details.


EHCD.py uses the following packages:
* Numpy
* Sympy (including sympy.physics.wigner)
* Concurrent (including concurrent.futures)
SCD.py uses the following packages:
* Numpy

The other files (starting with "HC") are meants as examples of how one may use the functions from EHCD.py and SCD.py in determining spectra of different types of Heisenberg chains. They use, among others, the EHCD.py and SCD.py module.

EHCD.py contains the following functions which one may want to use in their own code
 * coupling(spins): for determining all the ways of coupling a set of spins
    * Takes as input:
      * spins: a 1D Numpy array of positive half-integers which represents spins
    * Returns an array with all ways of coupling the given spins
 * mat_t(J, i, j, spins, coup): for determining reduced matrix elements of a term S_i * S_j for i not equal to j
    * Takes as input:
      * J: a Numpy array which represents the coupling matrix
      * spins: a 1D Numpy array of positive half-integers which represents spins
      * coup: a 2D numpy array with all the ways of coupling the spins. This type of input is returned by the coupling functions
      * i: index of a spin
      * j: index of a spin
    * Returns a 2D Numpy array with the reduced matrix element of S_i * S_j for the interactions between the subsystems in the direct sum decomposition that are given by the coupling sequences from coup multiplied by sqrt((2S^tot_i + 1)).
 * block(spins, J, coupling, S): #computes the diagonal block of the Heisenberg Chain corresponding to spin S
   * Takes as input:
     * spins: a 1D Numpy array of positive half-integers which represents spins
     * coupling: a 2D numpy array with all the ways of coupling the spins. This type of input is returned by the coupling functions
     * J: a Numpy array which represents the coupling matrix
     * S: a half integer representating a spin
   * Returns a 2D numpy array with all the reduced matrix elements of the Heisenberg hamiltonian for some spin S
 * uniline(N, J_s): Return the coupling matrix corresponding to the heisenberg line of N particles with coupling between neighbours equal to J_s
   * Takes as input:
     * N: the number of spins
     * J_s: The coupling between neighbouring spins
   * Returns a 2D Numpy array representing the coupling matrix corresponding to the heisenberg line of N particles with coupling between neighbours equal to J_s
 * unicircle(N, J_s): Return the coupling matrix corresponding to the heisenberg circle of N particles with coupling between neighbours equal to J_s
   * Takes as input:
     * N: the number of spins
     * J_s: The coupling between neighbouring spins
   * Returns a 2D Numpy array representing the coupling matrix corresponding to the heisenberg circle of N particles with coupling between neighbours equal to J_s

SCD.py contains the following functions which one may want to use in their own code
 * decompose: Computes the multiplictity of each spin in the direct sum decomposition of a system of coupled spin
   * Takes as input:
     * 1D Numpy array with half integers representing the spins of particles
   * Returns a 2D array with two collumns. the first collumns represents the multiplicity and the second the spin
  * reduction: Prints information about the decomposition of of a system of coupled spins
    * Takes as input:
     * 1D Numpy array with half integers representing the spins of particles
    * Does not return anything
