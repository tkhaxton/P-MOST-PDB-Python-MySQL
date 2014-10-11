
INTRODUCTION
------------

P-MOST-PDB-Map-Python-MySQL is software that uses Protein Data Bank structure 
files to calculate optimal backmapping parameters for the Protein Model with 
Oriented SiTes (P-MOST) and compute the root-mean-square displacement for each 
atom type.

LICENSE
-------

P-MOST-PDB-Map-Python-MySQL is free software: you can redistribute it and/or 
modify it under the terms of the GNU General Public License as published by the 
Free Software Foundation, either version 3 of the License, or (at your option) 
any later version.

P-MOST-PDB-Map-Python-MySQL is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
more details.

You should have received a copy of the GNU General Public License along with
P-MOST-PDB-Map.  If not, see <http://www.gnu.org/licenses/>.

PUBLICATION
-----------

Details of the oriented coarse-grained protein model appear in my paper:

[1] T. K. Haxton, "High-resolution coarse-grained modeling using oriented 
coarse-grained sites," submitted (available online at arxiv.org/abs/1409.8658).

AUTHORS, CITATION, AND CONTACT INFORMATION
-------------------------------------------

P-MOST-PDB-Map-Python-MySQL was developed by Tom Haxton at the Molecular 
Foundry, Lawrence Berkeley National Laboratory.

I ask that users please cite my publication [1] in any publication presenting 
results produced with P-MOST-PDB-Map-Python-MySQL.

Contact: Tom Haxton (tomhaxton@gmail.com)

PLATFORM
--------

P-MOST-PDB-Map has been tested on OS X.

REQUIREMENTS
------------

1. MySQL Community Server
2. Python
3. MySQL-python
4. PDB structure files (available from www.rcsb.org)

See also the C version, P-MOST-PDB-Map, available at 
http://nanotheory.lbl.gov/peptoid_code.

INSTALLATION
------------

None required.

CONTENTS/USAGE
--------------

Start the MySQL server, and run "sudo mysql" to open the SQL shell as root.

In a separate (bash) shell, run "python roundtrip.py <directory>", where
<directory> is the top directory containing Protein Data Bank (PDB) files.

MySQL OUTPUT TABLES
-------------------

pdb_file_error_table: List of PDB entries excluded from the calculation, as well
   as the reason for exclusion

pdb: List of atoms in the original PDB configurations

cg: Table of coarse-grained sites, listing their positions and orthonormal 
   vectors describing their orientations

cov: List of information for each linear correction backmapping function,
   including the intercept (3-vector) and slope (3 x 3 matrix)

rmsd: List of component-wise and overall root-mean-square displacement (rmsds)
   between original and backmapped atomic configurations, separated by amino
   acid residue and atom type
