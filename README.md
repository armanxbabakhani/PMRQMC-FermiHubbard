# PMRQMC-FermiHubbard
This code converts a FermiHubbard input into a spin-1/2 input file

### Format of the input file for the Fermi-Hubbard QMC

t = 1.0  
U = 2.0  
node1: neighbor11 neighbor12 ...  
node2: neighbor21 neighbor22 ...  
.  
.  
.  

### Description of the input file
Essentially, the first two lines specify the t and U parameter, and the rest of the file specifies the adjacency list for the graph on which Fermi-Hubbard model is computed.

### Description of what this code does
This code converts the Fermi-Hubbard input file into a QMC-spin1/2 input file. This is done so by creating two copies of the original lattice. 
One for spin up and the other for spin down. In other words, however many sites there are for the original lattice, there are double the number of sites for the QMC code.
The permutation terms do not mix, and only the diagonal term has a mixing of spin up and down.

## Code Execution
After compiling the FHtoQMC.cpp code, one can run the program through the command ./FHtoQMC "FHinput.txt", where "FHinput.txt" is the Fermi-Hubbard input file, with the format
described above. Then, the code will perform the conversion to a QMC spin-1/2 input file and create a "FHinput_QMC.txt" file.
