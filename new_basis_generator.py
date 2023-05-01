import dGr_orbitals
import sys
import get_basis as gb
import numpy as np
import scipy

old_basis=sys.argv[1]
orbs_file=sys.argv[2]

def GetSpacedElements(array,deg,initial):
    out=[]
    number=int(len(array)/deg)
    for i in range(number):
        out.append('{0:.6e}'.format(float(array[initial+i*deg])))
    return out

def Test_coef(array,threshold=10**-5):
    test_val=0
    for i in array:
        if abs(float(i)) > threshold:
            test_val=1
            break
    return test_val

def Norm(array):
    S = 0
    for i in range(len(array)):
        S+= float(array[i])**2
    #print('Ini S ',S)
    S = np.sqrt(S)
    for i in range(len(array)):
        #print(array[i],S,float(array[i])/S)
        array[i] = '{:.5e}'.format(float(array[i])/S)
        #print(array[i])
	
    

old_basis,exp_size,contr_size,sym_list,title=gb.get_basis(old_basis)
Co_orb = dGr_orbitals.Molecular_Orbitals.from_file(orbs_file)

#print(old_basis)

# Merge each de-contraction basis matrix in a block diagonal one of dimmension PxN.
for i in range(len(sym_list)):
     if i == 0:
         open_basis=old_basis[0]
     else:
         open_basis=scipy.linalg.block_diag(open_basis,old_basis[i])

Co_primt=np.matmul(open_basis,Co_orb[0])  #Decontracted the opt orbitals only for the most symmetrycal sym 
#Co_primt=np.matmul(open_basis,Co_orb[3]) #It can be easylly done to all sym

#print(Co_primt)

#for j in range(len(Co_primt[0][:])):
#    for i in range(len(Co_primt)):
#        print(i,j,Co_primt[i][j])
#exit()
primt_size=[]
primt_tot=0
# Print only the coef of the same sym in the molpro format
for i in range(len(sym_list)):
    print(title[i])
    primt_size.append(exp_size[i]*(2*i+1))
    Co_irr=Co_primt[primt_tot:primt_tot+primt_size[i],:]
    #print(Co_irr)
    for j in range(Co_primt.shape[1]):
        for k in range(2*sym_list[i]+1):
            #print(exp_size[i],j,k)
            Co_irr_no_deg=np.array(GetSpacedElements(Co_irr[:,j],(2*sym_list[i]+1),k),dtype=str)
            #print(len(Co_irr_no_deg))
            #print(Co_irr_no_deg)
            if Test_coef(Co_irr_no_deg) == 1:
#                Norm(Co_irr_no_deg)
                print('c, 1.'+str(exp_size[i])+', '+', '.join(Co_irr_no_deg))
                #pass
            #exit()
            #print(Co_irr[:,j])
        
    primt_tot+=primt_size[i]
#print(size)
