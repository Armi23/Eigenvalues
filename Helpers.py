# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 11:33:57 2013

@author: Armi, Nick, Hassan, Ibrahim
"""

from __future__ import division

import numpy as np
from random import randint
import math as m


#####################
# MAKING NEW THINGS #
#####################    
#Makes an empty matrix
def make_empty_matrix(): 
    return np.array([[]])
    
#Make a new matrix of size (rows,cols), allow default specifization
def new_matrix((rows, cols), default=0):
    if default == 0:
        return np.zeros((rows,cols))
    elif default == 1:
        return np.ones((rows,cols))
    else:
        print "Not allowed, giving you zeros"
        new_matrix(rows,cols,default=0)  
        
#makes matrix with random elements
def random_matrix(size,lower=0,upper=9): 
    matrix = new_matrix((size,size))
    for i in xrange(0,len(matrix)):
        for j in xrange(0,len(matrix[0])):
            update_el(matrix,i,j,randint(lower,upper))
    return matrix
    
#Make a new identity matrix
def new_identity(size):
    return np.identity(size)    

def new_vector(vals, default=0):
    if default == 0:
        return np.zeros(vals)
    elif default == 1:
        return np.ones(vals)
    else:
        print "Not allowed, giving you zeros"  
        
#makes matrix with random elements
def random_vector(size,lower=0,upper=9): 
    vector = new_vector(size)
    for i in xrange(0,len(vector)):
        update_el_v(vector,i,randint(lower,upper))
    return vector
    
#Makes a unit row of given size in given dimension
def unit_maker(size,dim=0): 
    result = make_empty_row()
    for i in xrange(0,size):
        if i == dim:
            result = add_el_to_r(result,1)
        else:
            
            result = add_el_to_r(result,0)
    return result
    
#####################
# UPDATING A MATRIX #
#####################
#Update rows and elements
def update_el(A,row,col,newval):
    A[row,col] = newval
def update_row(A,row,newvec):
    A[row] = newvec

#Adds a row to a matrix
def add_row_to_matrix(matrix, vector): 
    if (matrix.shape == (1,0)) :
        return np.array([vector])
    else:    
        return np.append(matrix, [vector],axis=0)

#Returns a matrix with chosen row removed
def remove_row(matrix, number): 
    result = make_empty_matrix()
    for i in xrange (0, len(matrix)):
        if i != number:
            result = add_row_to_matrix(result,get_row(matrix,i))
    return result

#Adds a column
def add_col_to_matrix(matrix,vector):  
    if (matrix.shape == (1,0)) :
        for i in xrange (0,len(vector)):
            el = make_empty_row()
            el = add_el_to_r(el,vector[i])
            matrix = add_row_to_matrix(matrix,el)
        return matrix
    else:    
        answer = make_empty_matrix()
        for i in xrange (0,len(vector)):
            row = add_el_to_r(get_row(matrix,i),vector[i])
            answer = add_row_to_matrix(answer,row)
        return answer

#Return a matrix without chosen column
def remove_col(matrix, number): 
    result = make_empty_matrix()
    for i in xrange (0, len(matrix[0])):
        if i != number:
            result = add_col_to_matrix(result,get_col(matrix,i))
    return result    

############################
# RETRIEVING FROM A MATRIX #
############################
#Returns i element of row
def get_el_r(row, i): 
    return row[i]

#Returns ij element of matrix
def get_el(matrix, i, j): 
    return matrix[i][j]

#Returns row of matrix
def get_row(matrix, row): 
    return matrix[row]
        
#Returns column (as a vector/row)
def get_col(matrix, col): 
    return matrix[:,col]
    
#########################
# VECTOR/ROW OPERATIONS #    
#########################
#Vectors are rows, to simplify and adhere to abstractions of functions
#Makes an empty row
def make_empty_row(): 
    return np.array([])
    
#use to expand row
def add_el_to_r(row, element): 
    return np.append(row,element)
    
#updates vector element
def update_el_v(v,i,newval):
    v[i] = newval

#retrieve from vector  
def get_el_v(row, i): 
    return row[i]
    
#changes length and keeps abstraction
def lengthen(vector, scalar):  
    def mapmult(element):
        return element * scalar
    npmult = np.vectorize(mapmult)
    return npmult(vector)     

#finds magnitude of vector
def magnitude(vector): 
    return (dot(vector,vector)**(0.5))

#Returns finds dot product of vectors (rows)
def dot(v1, v2): 
    Lenv1 = len(v1)
    Lenv2 = len(v2)
    if (Lenv1 != Lenv2):
        raise Exception("Cannot take dot product")
    else:     
        answer = 0
        for i in xrange(0, len(v1)):
            answer += v1[i] * v2[i]
        return answer
    
#returns unit vector of input vector
def normalize(vector): 
    return lengthen(vector,1/(magnitude(vector)))   

#adds 2 vectors and returns a vector
def add(v1,v2): 
    return v1 + v2

#projects v2 onto v1
def projection(v1,v2): 
    return lengthen(v1,(dot(v1,v2)/magnitude(v1)))    
 
#appends vector to another
def row_combine(v1,v2): 
    return np.append(v1,v2)

# compairs if 2 vectors are the same
def compare_vec(v1,v2): 
    return (v1 == v2).all()     
    
#Map over a matrix
def vecmap(f,A):
    elements, = shape(A)
    for i in xrange(elements):
        update_el_v(A,i,f(get_el_v(A,i)))
    return A

#Round each element of a matrix in place
def vecround(A,decplace):
    vecmap(lambda x: round(x,decplace),A)
    
#####################
# MATRIX OPERATIONS #   
#####################
#returns matrix with orthonormal columns (and rows)
def orthogonalize(matrix): 
    result = make_empty_matrix()
    for i in xrange(0,len(matrix)):
        col = get_col(matrix,i)
        rcol = col
        for j in xrange(0,len(get_row(result,0))):
            rcol = add(rcol,lengthen(projection(get_col(result,j),col),-1))
        result = add_col_to_matrix(result,normalize(rcol))
    return result
    
#matrix to a power
def power(matrix, num): # matrix to a power
    result = matrix
    for i in xrange(0,num-1):
        result = multiply(result,matrix)
    return result  
        
#Temporarily multiplies 2 matrices
def multiply(m1,m2): 
    return np.dot(m1,m2)

# finds the transpose of matrix depending on its shape`
def transpose_matrix(matrix):    
    if len(matrix) == 0:
        print("Not a matrix \n")
    elif len(matrix) == 1:
        if (len(matrix[0]) == 0):
            return make_empty_matrix()
        else:
            matrix_list = matrix.tolist()
            return np.array(map(lambda x : [x], matrix_list[0]))
    else:
        matrix_list = matrix.tolist()
        transpose_tup = zip(*matrix_list)
        return np.array(map(list, transpose_tup))

#Check if a matrix is balanced (n x m), assumes numpy array input
def check_balance(matrix):
    try:
        matrix.shape[1]
    except:
        raise Exception("Matrix not valid")

#Check if a matrix is square, assumes numpy array input
def check_square(matrix):
    #Try/catch gets unbalanced matrices
    check_balance(matrix)
    if matrix.shape[0] != matrix.shape[1]:
        raise Exception("Expected square matrix")

# Gets the shape of A under our abstraction, returns a tuple (rows, cols)
def shape(A):
    return A.shape

#Invert matrix A under our abstraction
def invert(A):
    return np.linalg.inv(A)
    
# multiplies all the entries of the matrix by scalar
def scale_mat(mat, scalar):
    return scalar*mat

# converts matrix into a list of lists
def matrix_to_list(mat):
    return mat.tolist()
    
#Map over a matrix
def matmap(f,A):
    rows,cols = shape(A)
    for i in xrange(rows):
        for j in xrange(cols):
            update_el(A,i,j,f(get_el(A,i,j)))
    return A

#Round each element of a matrix in place
def matround(A,decplace):
    matmap(lambda x: round(x,decplace),A)
    
##################
# MISC FUNCTIONS #  
##################
#Makes a matrix from user input
def make_into_matrix(inp):
    A = np.array(inp)
    check_balance(A)
    return A
    
def make_into_vec(inp):    
    A = np.array(inp)
    return A
    
#Solves polynomials
def roots(coeffs): 
    return np.roots(coeffs)

# Operates matrix A on vector b,
# assumes both are numpy arrays  
def operate(A, b):
    n = A.shape[0]  
    output = np.zeros(n)
    for i in xrange(n):
        output[i] = dot(A[i],b)
    return output

def nextPower (n):
    return 2**int(m.ceil(m.log(n,2)))
    
if __name__ == "__main__":   
    print "Testing add row to matrix"
    A = make_into_matrix([[1,2,3],[4,5,6],[7,8,9]])
    B = new_vector(3,1)
    Rowtest = make_into_matrix([[1,2,3],[4,5,6],[7,8,9],[1,1,1]])
    assert (add_row_to_matrix(A,B) == Rowtest).all()
    
    print "Testing add column to matrix"
    Coltest = make_into_matrix([[1,2,3,1],[4,5,6,1],[7,8,9,1]])
    assert (add_col_to_matrix(A,B) == Coltest).all()
    
    print "Testing remove row"
    RemoveRtest = make_into_matrix([[1,2,3],[7,8,9]])
    assert (remove_row(A,1) == RemoveRtest).all()
    
    print "Testing remove col"
    RemoveCtest = make_into_matrix([[1,3],[4,6],[7,9]])
    assert (remove_col(A,1) == RemoveCtest).all()
    
    print "Testing dot product"
    C = random_vector(3)
    D = random_vector(3)
    assert dot(C,D) == np.dot(C,D)
    
    print "Testing dot product"
    C = random_vector(3)
    D = random_vector(3)
    assert dot(C,D) == np.dot(C,D)
    
    print "Testing lengthen"
    assert (lengthen(C,3) == np.dot(C,3)).all()
    
    print "Testing magnitude"
    assert (magnitude(C) == np.linalg.norm(C))
    
    print "Testing normalize"
    Norm = normalize(C)
    vecround(Norm,2)
    check = C/np.linalg.norm(C)
    vecround(check,2)
    assert (Norm == check).all()
    
    print "Testing projection"
    assert (projection(C,D) == np.dot(C,(np.dot(C,D)/np.linalg.norm(C)))).all()
    
    print "Testing orthogonalize (makes absolute value because of direction)"
    Testmatrix = random_matrix(5)
    Qtest = abs(orthogonalize(Testmatrix))
    matround(Qtest,6)
    q,r = np.linalg.qr(Testmatrix)
    q = abs(q)
    matround(q,6)
    assert (Qtest == q).all()
    
    print "Testing power"
    p = power(Testmatrix,3)
    pcheck = np.linalg.matrix_power(Testmatrix,3)
    assert (p == pcheck).all()
    
    print "Testing unit_maker"
    row1 = unit_maker(3,0)
    row2 = unit_maker(3,1)
    test = make_empty_matrix()
    test = add_row_to_matrix(test,row1)
    test = add_row_to_matrix(test,row2)
    check = make_into_matrix([[1,0,0],[0,1,0]])
    assert (test == check).all()
