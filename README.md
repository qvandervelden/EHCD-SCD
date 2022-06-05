# EHCD (Efficient Heisenberg Chain Diagonalization)
Software for determining the spectra of Heisenberg chains

This folder contains the following file

  • EHCD.py: Python file containing functions to determine the blocks of the matrix of a Heisenberg chain disgonalized with respect to total spin.
  
  • SCD.py: Python file containing functions to study the decomposition of a system of coupled spins
  
  • HC-Circle.py: Python file for determining the spectra of a Heisenberg chain with particles in a line, coupled to their neighbours and the particles at both ends also coupled, all with the same coupling strenght.
  
  • HC-General.py: Python file for determining the spectra of a general Heisenberg chain.
  
  • HC-Line-bipartite.py: Heisenberg chain with particles in a line, coupled to their neighbours. The particles can be divided in two sets A en B with particles in the same set being coupled ferromagnetically with strength J_s and particles of different sets coupled antiferromagnetically with strength -J_s
  
  • HC-Circle.py: Python file for determining the spectra of a Heisenberg chain with particles in a line, coupled to their neighbours, all with the same coupling strenght.


All the files above above were written by Quinn van der Velden (q.vandervelden@students.uu.nl). The code is made available under the GPLv3 License. See the
file COPYING for details.


EHCD.py uses the following packages:
• Numpy
• Sympy (including sympy.physics.wigner)
• Concurrent (including concurrent.futures)
SCD.py uses the following packages:
• Numpy

The other files (starting with "HC")
