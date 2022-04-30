from __future__ import division
import numpy as np
import scipy.linalg as slg
import math
import sys
from numpy import genfromtxt
import csv
from numpy import linalg as LA
from numpy.linalg import inv

def str2float(s):
    return float(s.replace('D','e'))

def CutLTMatOut(fname, head_line, index=0):
    cur_index = 0
    f = open(fname)
    status = None
    r = 0
    indices = []
    mat = []
    for line in f:
        if head_line in line:
            status = "Column Indices"
            continue
        if status == "Rows":
            splitted = line.split()
            if all([x.isnumeric() for x in splitted]):
                status = "Column Indices"
            else:
                try:
                    r = int(splitted[0])
                    floats = [str2float(s) for s in splitted[1:]]
                    if len(mat) < r:
                        mat.append(floats)
                    else:
                        mat[r-1] += floats
                except ValueError as e:
                    break
                continue
        if status == "Column Indices":
            status = "Rows"
            indices = [int(I) for I in line.split()]
            continue
    f.close()
    return mat

def toNPMatrix(mat_array, symmetry = 'symmetric'):
    dim = len(mat_array)
    mat = np.zeros((dim,dim))
    for i in range(dim):
        mat[i,:i+1] = mat_array[i]
    if symmetry == 'symmetric':
        diag = np.diag(np.diag(mat))
        return mat + mat.T - diag
    elif symmetry == 'anti-symmetric':
        return mat - mat.T


def ReadMatrix(fname, head_line, index=0):
    return np.array(CutLTMatOut(fname, head_line, index))

def ReadLowerTri(fname, head_line, index=0, symmetry = 'symmetric'):
    return toNPMatrix(CutLTMatOut(fname, head_line, index), symmetry)

def outMatrix(A):
    if np.iscomplexobj(A):
        print("Real part:")
        outMatrix(A.real)
        print()
        print("Imag part:")
        outMatrix(A.imag)
        return
    print("--------------------------------------------------------------------------------")
    curColStart, curColEnd = 0, 5
    while curColStart < A.shape[1]:
        curColEnd = min(curColStart + 5, A.shape[1])
        print("     "+"               ".join(list(map(str, range(curColStart + 1, curColEnd + 1)))))
        for i in range(A.shape[0]):
            print("%3d   "%(i+1,) + " ".join(list(map(lambda x : "%15.8e"%(x if abs(x) > 1e-8 else 0,) , A[i, curColStart : curColEnd]))))
        curColStart += 5
        print()

    print("--------------------------------------------------------------------------------")

fname = "./hcn_ccpvdz_pb4d.log"
head_line = '*** NEO Overlap ***'

smat = ReadLowerTri(fname, head_line)

print("Overlap Matrix: ")
outMatrix(smat)

# Diagonalize smat, and obtain the eigenvalues and eigenvectors
Sigma_S, L_S  = LA.eig(smat)

print("Eigenvalues of smat:")
print(Sigma_S)

print('Determinant of smat', np.linalg.det(smat))



fname = "./hcn_ccpvdz_pb4.log"
head_line = '*** NEO Overlap ***'

smat_nd = ReadLowerTri(fname, head_line)

print("Overlap Matrix (no D): ")
outMatrix(smat_nd)

# Diagonalize smat, and obtain the eigenvalues and eigenvectors
Sigma_S_nd, L_S_nd  = LA.eig(smat_nd)

print("Eigenvalues of smat (no D):")
print(Sigma_S_nd)

print('Determinant of smat (no D):', np.linalg.det(smat_nd))





#Sigma_S_minusOneHalf = (np.diag(Sigma_S**(-0.5)) )
#S_minusOneHalf = np.dot(L_S, np.dot(Sigma_S_minusOneHalf, np.transpose(L_S)))
#
#print("V: ")
#outMatrix(S_minusOneHalf)
#
#iden = np.dot(np.dot(S_minusOneHalf.transpose(), gs), S_minusOneHalf)
#
#print("Should be identity: ")
#outMatrix(iden)
#
#print("V+SV")
#
#iden2 = np.dot(np.dot(gv.transpose(), gs), gv)
#outMatrix(iden2)
