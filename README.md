# TopoProject
An algorithm to compute the dual (marked) Thurston unit ball given an element of the mapping class group and a curve on the surface of the 2-torus.

INSTRUCTIONS:

Firstly you will need to install the python modules (can be found by google-searching):

scipy

shapely

matplotlib

.
The generators of the mapping class group are called a, A, ..., d, D where a, b, c, d are the usual generators and the 
capitalisation indicates inversion. 

The curve on the surface is to be inputted as a word in x, X, y, Y, w, W, z, Z
where these generators are as in the literature you gave me, and again capitalisation indicates inversion.

To produce the dual unit ball, you will use the MCG_to_unitball function. This takes two arguments. The first is the surface
curve that you want to deform. Input this in the form [x, y, X, Y] where this would represent the word xyXY.
The second is the element of the MCG that you want to use to deform your curve. Input this in the form 
[A, c, D] where this represents the word AcD in the MCG.

To produce the MARKED dual unit ball, you use the DualWithMarked function. The inputs work the same way as before. The only 
reason I have left the MCG_to_unitball function in is that it runs much faster for some large examples I tested. 


