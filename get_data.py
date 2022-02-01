
import numpy as np
import h5py
from numpy.fft import fft, fftfreq, ifft
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

def prettyPrint(matrix):
    print('\n'.join(['\t'.join([str(cell) for cell in row]) for row in matrix]))

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

## Locate the data stored in the HDF5 groups, and put them into arrays
#with h5py.File('coh2_sto3g_et.bin', 'r') as f_x:
#    base_items = list(f_x.items())
#    print('Items in the SCF directory: \n', base_items)
#    SCF = f_x.get('SCF')
#    SCF_items = list(SCF.items())
#    print('\n Items in SCF: ', SCF_items)
#    fock = np.array(SCF.get('FOCK_SCALAR'))
#    #time =  np.array(RT.get('TIME'))
#
# Locate the data stored in the HDF5 groups, and put them into arrays
with h5py.File('coh2_sto3g_et.bin', 'r') as f_x:
    base_items = list(f_x.items())
#    print('Items in the base directory: \n', base_items)
    PINTS = f_x.get('PINTS')
    INTS = f_x.get('INTS')
    PINTS_items = list(PINTS.items())
    INTS_items = list(INTS.items())
#    print('\n Items in PINTS: ', PINTS_items)
    ps = np.array(PINTS.get('OVERLAP'))
    pt = np.array(PINTS.get('KINETIC'))
    pv = np.array(PINTS.get('POTENTIAL'))
    es = np.array(INTS.get('OVERLAP'))
    et = np.array(INTS.get('KINETIC'))
    ev = np.array(INTS.get('POTENTIAL'))
    eh = et+ev
    #time =  np.array(RT.get('TIME'))


print("NEO Overlap: ")
outMatrix(ps)
print("NEO Kinetic: ")
outMatrix(pt)
print("NEO Potential: ")
outMatrix(pv)
print("Electronic Overlap: ")
outMatrix(es)
print("Electronic Kinetic: ")
outMatrix(et)
print("Electronic Potential: ")
outMatrix(ev)
print("Electronic H")
outMatrix(eh)
