# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 11:58:24 2013

@author: Armi
"""

from __future__ import division
from decimal import *
getcontext().prec = 100

import Helpers as H
import scipy.linalg as sc
import numpy as np

################
# Trace Section
################
#adds up diagonal elements
def trace(matrix): 
    H.check_square(matrix)
    answer = 0    
    for i in xrange(0, len(matrix)):
        answer += H.get_el(matrix,i,i)
    return answer
    
#########################################
# Determinant Section (simple functions)
#########################################
#returns determinant of 2x2 matrix
def _simple_determinant(matrix): 
    a = H.get_el(matrix,0,0)
    b = H.get_el(matrix,0,1)
    c = H.get_el(matrix,1,0)
    d = H.get_el(matrix,1,1)
    return a*d - b*c

#Regular Laplace co-factors
def _R_determinant(matrix):    
    if len(matrix) == 2:
        return _simple_determinant(matrix)
    else:
        answer = 0
        for i in xrange(0, len(matrix)):
            el = H.get_el(matrix,0,i)            
            smaller = H.remove_col(matrix,i)
            smaller = H.remove_row(smaller,0)
            answer += el*(-1.0)**(i+1)*_R_determinant(smaller)
    return answer
    
#Swtiches signs based on length of matrix
def _sign_determinant(matrix): 
    return (-1.0)**(len(matrix))*_R_determinant(matrix)
    
#Determinant of triangular matrix
def _T_determinant(matrix): 
    answer = 1.0
    for i in xrange(0, len(matrix)):
        answer *= H.get_el(matrix,i,i)
    return answer 

########################################
# Determinant Section (QR Decomposition)
########################################   
#makes matrix based on v*vT (explanations are below)
def _vvT(vector): 
    matrix = H.make_empty_matrix()
    for i in xrange(0,len(vector)):
        row = H.make_empty_row()
        for j in xrange(0,len(vector)):
            el = H.get_el_r(vector,i)*(H.get_el_r(vector,j))
            row = H.add_el_to_r(row,el)
        matrix = H.add_row_to_matrix(matrix,row)
    return matrix
    
# pads smaller matrices for QR determinant
def _padding(matrix,size): 
    result = H.make_empty_matrix()
    pad = size - len(matrix)
    for i in xrange(0,pad):
        unit = H.unit_maker(size,i)
        result = H.add_row_to_matrix(result,unit)
    pad_row = H.new_vector(pad)
    for i in xrange (0,len(matrix)):
        row = H.row_combine(pad_row,H.get_row(matrix,i))
        result = H.add_row_to_matrix(result,row)
    return result
    
#makes Qi based on I-2vvT
def _make_Q(vector,size): 
    I = H.new_identity(len(vector))
    V = _vvT(vector)
    Q = I - 2*V
    QP = _padding(Q,size)
    return (Q,QP)
    
# sets up Q maker and returns Q to QR determinant
def _Q_setup(matrix,size): 
    e = H.unit_maker(len(matrix),0)
    col = H.get_col(matrix,0)
    a = H.magnitude(col)
    if (H.get_el(matrix,0,0) > 0):
        a *= -1.0
    ae = H.lengthen(e,a)
    u = H.add(col,ae)
    if (u == H.new_vector(len(matrix))).all():
        V = H.normalize(col)
        (Q,QP) = _make_Q(V,size)
        return (Q,QP)      
    else:   
        V = H.normalize(u)
        (Q,QP) = _make_Q(V,size)
    return (Q,QP)
            
'''            
This is the O(n^3) method with sign. Diagnolizes each column by reducing 
to one value and zeros and then moving on to the next minor matrix. Read 
more about it in this link. Many variables usually share names between the 
ones used in the link and the ones in the program.
http://en.wikipedia.org/wiki/QR_decomposition#Using_Householder_reflections
'''            
def _QR_determinant(matrix):    
    Qlist = []
    A = matrix
    size = len(matrix)
    for i in xrange(0,len(matrix)-1):
        (Q,QP) = _Q_setup(A,size)        
        Qlist.append(QP)
        A = H.multiply(Q,A)
        A = H.remove_col(A,0)
        A = H.remove_row(A,0)
    result = matrix
    for i in xrange(0,len(Qlist)):
        result = H.multiply(Qlist[i],result)
    return (-1.0)**(len(matrix)+1)*round(_T_determinant(result),4)
    
'''
This is absolute value QR decomposition. It runs much quicker than QR 
decomposition because it uses another method. You can read about it here, 
but the names of the functions explain it rather well. 
http://en.wikipedia.org/wiki/QR_decomposition#Using_the_Gram.E2.80.93Schmidt_process
'''
def _AQR_determinant(matrix): 
    Q = H.orthogonalize(matrix)
    QT = H.transpose_matrix(Q)
    R = H.multiply(QT,matrix)
    Rdet =_T_determinant(R)
    return round(Rdet,4)
    
#finds determinant by calling best method
def determinant(matrix): 
    H.check_square(matrix)
    if len(matrix) < 6:
        return _sign_determinant(matrix)
    else:
        return _QR_determinant(matrix)
        
#finds absolute value of the determinant by calling best method
def absdeterminant(matrix): 
    H.check_square(matrix)
    if len(matrix) < 6:
        return abs(_R_determinant(matrix))
    else:
        return _AQR_determinant(matrix)  
    
#############################################
# Newton's Identities section (first attempt)
# Built each coefficent based on previous
############################################# 
'''
This function find coefficients of characteristic polynomials by using
Newton's Identities. You can read about them here:
http://en.wikipedia.org/wiki/Newton%27s_identities#Formulation_in_terms_of_symmetric_polynomials
However, they are a numerically unstable way to calculate coefficients. 
They are refered to as Leverrier's Method in the following paper which explains
why they are numerically unstable:
http://www4.ncsu.edu/~ipsen/ps/charpoly3.pdf
'''
def _N_identity(matrix): 
    elist = [Decimal(1)]
    mlist = matrix
    plist = [Decimal(trace(matrix))]
    elist.append(plist[0])
    for i in xrange(1,len(matrix)-1):
        mlist = H.multiply(matrix,mlist)
        plist.append(Decimal(trace(mlist)))
        exp = Decimal(0)
        for j in xrange(len(elist)):
            sign = Decimal((-1)**j)
            e = elist[len(elist)-j-1]
            p = plist[j]
            exp += sign*e*p
        elist.append(exp/Decimal(i+1))
    return elist

################################################
# Newton's Identities section (second attempt)
# Makes system of equations to solve for coeffs
################################################
#makes list of traces
def _tlist_maker(matrix): 
    tlist = H.make_empty_row()
    tlist = H.add_el_to_r(tlist,trace(matrix))
    powerA = matrix
    for i in xrange(len(matrix)-1):
        powerA = H.multiply(matrix,powerA)
        tlist =  H.add_el_to_r(tlist,trace(powerA))
    return tlist

#makes matrix that we need to solve
def _system_maker(tlist): 
    system = H.make_empty_matrix()
    for i in xrange(len(tlist)):
        col = H.new_vector(i)
        col = H.add_el_to_r(col,i+1)
        j = 0
        while len(col) < len(tlist):
            col = H.add_el_to_r(col,tlist[j])
            j += 1
        system = H.add_col_to_matrix(system,col)
    return system

'''   
This is the second attempt to make coeffs using traces. This also failed,
but you can read about how it works on page 7 of the following pdf.
http://scipp.ucsc.edu/~haber/ph116A/charpoly_11.pdf
'''
def _N_system(matrix): 
    tlist = _tlist_maker(matrix)
    system = _system_maker(tlist)
    coeffs = np.linalg.solve(system,tlist)
    return coeffs                

#########################################################
# Root (solutions) list manipulation section
#########################################################
#rounds a list of numbers
def _rounds(lst,num): 
    def rounding(element):
        return round(element,num)
    return map(rounding,lst)
    
#rounds a list of imaginary numbers
def _irounds(lst,num,inum): 
    result = []
    def rounding(element):
        real = round(element.real,num)
        imag = round(element.imag,inum)
        result.append((real,imag))
    map(rounding,lst)
    return result
    
#seperates imaginary and real solutions to desired accuracy(acc and iacc)
def _splits(lst,acc,iacc): 
    real = []
    imag = []
    for i in xrange(0,len(lst)):
        if round(lst[i].imag,acc) == 0:
            real.append(lst[i].real)
        else:
            imag.append(lst[i])
    real = _rounds(real,acc)
    imag = _irounds(imag,acc,iacc)
    return (real,imag)
    
#########################################################
# Coefficient making section. It calls _N_identity rather
# than _N_system because system fails quicker
#########################################################
#makes coefficients for characteristic polynomial 
def _coeffs(matrix):    
    newton = _N_identity(matrix)
    coeffs = []
    for i in xrange(0,len(newton)):
        coeffs.append(Decimal((-1)**(i+1))*newton[len(newton)-i-1])
    det = Decimal(determinant(matrix))
    coeffs.reverse()
    coeffs.append(det)
    coeffs = _rounds(coeffs,6)
    return coeffs

####################################################################
# Characteristic and Eigenvalue both call coeffs if the matrix is
# small enough to use my methods, otherwise it finds them through
# professional packages
####################################################################
#Useful to put +/- for equation
def _string_num(num): 
    result = ""    
    if num >= 0:
        result += " + "
    else: 
        result += " - "    
    result += str(abs(num))
    return result

#builds coefficients
def _build(coeffs): 
    equation = ""
    for i in xrange(0,len(coeffs)):
        if i == 0:
            equation += str(coeffs[i]) + "x^" + str(len(coeffs)-i-1)
        elif i == len(coeffs)-1:
            equation += _string_num(coeffs[i])  
        else:
            equation += _string_num(coeffs[i]) + "x^"
            equation += str(len(coeffs)-i-1)
    return equation

#returns the characteristic equation, as string
def characteristic(matrix): 
    H.check_square(matrix)
    if len(matrix) < 15:    
        coeffs = _coeffs(matrix) 
    else: 
        coeffs = np.poly(matrix)
    return _build(coeffs)
    
#returns list of real and complex eigenvalues separately
def eigenvalue(matrix,acc=4,iacc=4): 
    H.check_square(matrix)
    if len(matrix) < 15:
        coeffs = _coeffs(matrix)
        solutions = H.roots(coeffs)
    else:
        solutions = sc.eigvals(matrix)
    (real,imag) = _splits(solutions,acc,iacc)
    return real,imag
    
if __name__ == "__main__":
    print "Testing trace"    
    test1 = H.random_matrix(5)
    assert round(trace(test1),0) ==  round(np.trace(test1),0)   
    
    print "Testing Laplace Co-factors determinant"    
    assert round(determinant(test1),0) ==  round(np.linalg.det(test1),0)
    
    print "Testing QR decomposition determinant"    
    test2 = H.random_matrix(10)
    assert round(determinant(test2),0) == round(np.linalg.det(test2),0)
    
    print "Testing eigvalues (sometimes fails due to instability and order)"    
    solutions = sc.eigvals(test2)
    a = _splits(solutions,2,2)
    print eigenvalue(test2,2,2) == a
    
    print "Testing characteristic (sometimes fails due to instability)"    
    coeffs = np.poly(test2)
    a = _build(coeffs)
    print characteristic(test2) == a