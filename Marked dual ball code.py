import operator
import scipy
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Polygon
from shapely.geometry import MultiPoint, Point

def WalkOnList(L):
    #produces the walk determined by adding points together succesively
    length = len(L)
    pts = [(0,0)]
    for k in range(length):
        pts.append(tuple(map(operator.add, pts[-1], L[k])))
    return(pts)  

def HullPlot(L):
    #plots the convex hull of the walk above
    pts = WalkOnList(L)
    newpts = np.asarray(pts)
    hull = ConvexHull(newpts)
    plt.plot(newpts[:, 0], newpts[:,1], 'o')
    for simplex in hull.simplices:
            plt.plot(newpts[simplex, 0], newpts[simplex, 1], 'k-')
            
from sympy import symbols, Symbol
x, y, z, w, X, Y, Z, W = symbols('x, y, z, w, X, Y, Z, W')

#These coming functions a, b, ... C, D represent the standard generators of the MCG and their inverses

def a(t):
    #input as a symbol
    if t == Symbol('x'):
        return (Symbol('Z'), Symbol('x'))
    if t == Symbol('X'):
        return (Symbol('X'), Symbol('z'))
    else:
        return(t)
        
def b(t):
    #input as a symbol
    if t == Symbol('z'):
        return (Symbol('x'), Symbol('z'))
    if t == Symbol('Z'):
        return (Symbol('Z'), Symbol('X'))
    else:
        return(t)
        
def c(t):
    #input as a symbol
    if t == Symbol('x'):
        return (Symbol('x'), Symbol('w'), Symbol('Z'))
    if t == Symbol('y'):
        return (Symbol('y'), Symbol('w'), Symbol('Z'))
    if t == Symbol('X'):
        return (Symbol('z'), Symbol('W'), Symbol('X'))
    if t == Symbol('Y'):
        return (Symbol('z'), Symbol('W'), Symbol('Y'))
    else:
        return(t)
    
def d(t):
    #input as a symbol
    if t == Symbol('w'):
        return (Symbol('Y'), Symbol('w'))
    if t == Symbol('W'):
        return (Symbol('W'), Symbol('y'))
    else:
        return(t)
        
def e(t):
    #input as a symbol
    if t == Symbol('y'):
        return (Symbol('w'), Symbol('y'))
    if t == Symbol('Y'):
        return (Symbol('Y'), Symbol('W'))
    else:
        return(t)
        
def A(t):
    #input as a symbol
    if t == Symbol('x'):
        return (Symbol('z'), Symbol('x'))
    if t == Symbol('X'):
        return (Symbol('X'), Symbol('Z'))
    else:
        return(t)
        
def B(t):
    #input as a symbol
    if t == Symbol('z'):
        return (Symbol('X'), Symbol('z'))
    if t == Symbol('Z'):
        return (Symbol('Z'), Symbol('x'))
    else:
        return(t)
        
def C(t):
    #input as a symbol
    if t == Symbol('x'):
        return (Symbol('x'), Symbol('z'), Symbol('W'))
    if t == Symbol('y'):
        return (Symbol('y'), Symbol('z'), Symbol('W'))
    if t == Symbol('X'):
        return (Symbol('w'), Symbol('Z'), Symbol('X'))
    if t == Symbol('Y'):
        return (Symbol('w'), Symbol('Z'), Symbol('Y'))
    else:
        return(t)

def D(t):
    #input as a symbol
    if t == Symbol('w'):
        return (Symbol('y'), Symbol('w'))
    if t == Symbol('W'):
        return (Symbol('W'), Symbol('Y'))
    else:
        return(t)
        
def E(t):
    #input as a symbol
    if t == Symbol('y'):
        return (Symbol('W'), Symbol('y'))
    if t == Symbol('Y'):
        return (Symbol('Y'), Symbol('w'))
    else:
        return(t)
    
def ListOp(f, L):
    l = [f(v) for v in L]
    return(l)

def SlowReduce(word):
    #reduces a word in the mapping class group only so far as killing elements and their inverses when next to one another
    x, y, z, w, X, Y, Z, W = symbols('x, y, z, w, X, Y, Z, W')
    kill_list = [(x, X), (y, Y), (z, Z), (w, W), (X, x), (Y, y), (W, w), (Z, z)]
    bad_indices = []
    for i in range(len(word) - 1):
        elt = word[i]
        next_elt = word[i+1]
        if (elt, next_elt) in kill_list:
            bad_indices.append(i)
            bad_indices.append(i+1)
            i+=2
    wanted_indices = [i for i in range(len(word)) if i not in bad_indices]
    reduced_word = [word[i] for i in wanted_indices]
    return(reduced_word)

def CyclicReduce(word):
    #continues on from SlowReduce to cyclically reduce the word
    x, y, z, w, X, Y, Z, W = symbols('x, y, z, w, X, Y, Z, W')
    kill_list = [(x, X), (y, Y), (z, Z), (w, W), (X, x), (Y, y), (W, w), (Z, z)]
    word = SlowReduce(word)
    length = len(word) - 1
    while (word[length], word[0]) in kill_list:
        word = [word[i] for i in range(1, length)]
        length -= 2
    return(word)

def FullClean(word):
    #gets rid of unnecessary symbols in a word in the MCG for the purposes below. We can disregard the letters z,w,Z,W by...
    #a result I proved
    x, y, z, w, X, Y, Z, W = symbols('x, y, z, w, X, Y, Z, W')
    word = [v for v in word if v == x or v == y or v == X or v == Y]

    new_word = CyclicReduce(word)
    return(new_word)

def map_to_vecs(word):
    #maps a word in the MCG to their corresponding vectors as used by Thurston
    x, y, z, w, X, Y, Z, W = symbols('x, y, z, w, X, Y, Z, W')
    word = FullClean(word)
    list_of_vecs = []
    for elt in word:
        if elt == x:
            list_of_vecs.append([1,0])
        if elt == y:
            list_of_vecs.append([0,1])
        if elt == X:
            list_of_vecs.append([-1,0])
        if elt == Y:
            list_of_vecs.append([0,-1])
    new_list = [np.asarray(v) for v in list_of_vecs]
    return(new_list)

def ReduceList(L):
   
    reduced_list = []
    for element in L:
        if isinstance(element, tuple):
            for val in element:
                reduced_list.append(val)
        else:
            reduced_list.append(element)
    return reduced_list

def DehnOnCurve(L, gamma):
    #gives reusult of action of a sequence of Dehn twists (given by a word in the MCG) on a curve on the surface
    for func in L[::-1]: 
        gamma = SlowReduce(ReduceList(ListOp(func, gamma)))
        print(gamma)
    return(map_to_vecs(gamma))

            
def MarkedList(L):
    #gives list of the marked vertices
    pts = WalkOnList(L)
    M = []
    for v in pts:
        check = 0
        for w in pts:
            if w == v:
                check += 1
        if check == 1:
            M.append(v)
    return(M)      
    

def DualUnitBall(L):
    #produces the dual unit ball for an element of the MCG
    pts = WalkOnList(L)
    newpts = np.asarray(pts)
    hull = ConvexHull(newpts)
    poly = MultiPoint(pts).convex_hull
    MD = []
    for i in hull.vertices:
        vertex = np.asarray(pts[i])
        neighbours = [vertex + (0,1), vertex + (1,1), vertex + (1,0), vertex + (1,-1), vertex + (0,-1), vertex + (-1,-1), vertex + (-1,0), vertex + (-1,1)]
        a1 = Point(tuple(neighbours[0]))
        a2 = Point(tuple(neighbours[1]))
        a3 = Point(tuple(neighbours[2]))
        a4 = Point(tuple(neighbours[3]))
        a5 = Point(tuple(neighbours[4]))
        a6 = Point(tuple(neighbours[5]))
        a7 = Point(tuple(neighbours[6]))
        a8 = Point(tuple(neighbours[7]))
        if poly.intersects(a1):
            if poly.intersects(a2) and poly.intersects(a3):
                MD.append(vertex + (0.5, 0.5))
            if poly.intersects(a7) and poly.intersects(a8):
                MD.append(vertex + (-0.5, 0.5))
        if poly.intersects(a5):
            if poly.intersects(a3) and poly.intersects(a4):
                MD.append(vertex + (0.5, -0.5))
            if poly.intersects(a6) and poly.intersects(a7):
                MD.append(vertex + (-0.5, -0.5))
    return(MD)
    
def DualMarkedVertices(L):
    #produces marked vertices for the DUAL unit ball for an element of the MCG
    DM = []
    marked = MarkedList(L)
    for v in DualUnitBall(L):
        nbs = [v + (0.5, 0.5), v + (0.5, -0.5), v + (-0.5, -0.5), v + (-0.5, 0.5)]
        for i in range(3):
            if tuple(nbs[i]) in marked:
                DM.append(v)
                break
    return(DM)

def MCG_to_unitball(L, gamma):
    #produces the dual unit ball for the result of a sequence of Dehn twists on a surface curve
    return(DualUnitBall(DehnOnCurve(L, gamma)))

def DualWithMarked(L, gamma):
    #plots the marked dual unit ball for the result of a sequence of Dehn twists on a surface curve
    verts = np.array(DehnOnCurve(L, gamma))
    marked = np.array(DualMarkedVertices(verts))
    
    hull2 = ConvexHull(verts)
    plt.plot(verts[:, 0], verts[:,1], 'o')
    for simplex in hull2.simplices:
            plt.plot(verts[simplex, 0], verts[simplex, 1], 'k-')
    plt.scatter(marked[:, 0], marked[:, 1], 'o')
    plt.show()
