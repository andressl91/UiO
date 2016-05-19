This program assumes that the sides of the matrices divided by the size of the process grid is a 
perfect square. 

This can very easily be fixed by changing from a 2D - partition to a 1D.

An Alternative is to pad the matrix with enough zeros so the size becomes nice enough for the  
code to work. One could work wit hthis extended matrix throughout, and only write back the relevant  
part, or one could let each process know the relevant sloce of its block. This last point would  
require a rewrite of the implementatin of Cannon's algorithm, by exchanging MPI_Isendrecv_replace  
with MPI_Isend, MPI_Irecv, MPI_probe and MPI_Get_count.   

