import numpy as np

def orbital_sym(l_or_s,value):
    '''Given the quantum number l (l) or the symmetry (s) return the other
    '''
    orbital_list={ 0 : 's',
                   1 : 'p',
                   2 : 'd',
                   3 : 'f',
                   4 : 'g',
                   5 : 'h',
                   6 : 'i',
                   7 : 'j'}
    if l_or_s == 'l':
        try:
            out=orbital_list[value]
        except:
            return None
        return out
    elif l_or_s == 's':
        try:
            out=list(orbital_list.keys())[list(orbital_list.values()).index(value)]
        except:
            return None
        return out
    raise ValueError('Error in get the orbital symmetry or number')


def get_basis(file_name):
    '''recieve a basis set molpro input file and return a matrix with
       the contracted coefficients for each contracted basis function
       
       return: a set of arrays described below
       
       - 'basis' conatining len(sym_list) matrix are generated
       one for eache symmtry
       
       Matrix structure for sym n:
       lines: P_n*(2*n+1)
         where P_n is the number of primitive with sym n and (2*n+1)
         the degenerescence factor;
       
       columns: N_n*(2*n+1),
         where N_n is the number of contracted with sym n and (2*n+1)
         the degenerescence factor.
       
       - 'exp' number of primitive functions for each sym; 
       - 'contr' number of contracted functions for each sym;
       - 'sym_list' order of the sym in for the other arrays;
       - 'sym_title' the molpro basis expoents for each sym.

    '''
    f=open(file_name,'r')
    temp=[]
    basis=[]
    exp=[]
    contr=[]
    sym_list=[]
    sym_title=[]
    sym=-1
    for line in f:
        line=line.rstrip('\n')
        test=orbital_sym('s',line[0])
        if test is not None:
            if sym >= 0:
                sym_list.append(sym)
                contr_size=len(temp)
                basis.append(temp_matrix_converter(temp,exp_size,contr_size,sym,sym_title[-1]))
                contr.append(contr_size)
                temp=[]
            exp_size=len(line.split(','))-2
            exp.append(exp_size)
            sym_title.append(line)
            sym=test
        elif line[0] == 'c':
             temp.append(line) 
    sym_list.append(sym)
    contr_size=len(temp)
    contr.append(contr_size)
    basis.append(temp_matrix_converter(temp,exp_size,contr_size,sym,sym_title[-1]))
    return basis,exp,contr,sym_list,sym_title


def temp_matrix_converter(temp,exp_size,contr_size,s,exp_title):
    '''Generate a matrix for the sym n:
       Matrix structure for the p case:
       contracted:     first      second    . . .     N_p-th
                    px  py  pz  px  py  pz  . . . px   py   pz  
       g(e_1^p)    |c11  0   0  c12  0   0  . . . c1Np  0    0   |
       g(e_1^p)    | 0  c11  0   0  c12  0  . . .  0   c1Np  0   |
       g(e_1^p)    | 0   0  c11  0   0  c12 . . .  0    0   c1Np |
       g(e_2^p)    |c21  0   0  c22  0   0  . . . c2Np  0    0   |
       g(e_2^p)    | 0  c21  0   0  c22  0  . . .  0   c2Np  0   |
       g(e_2^p)    | 0   0  c21  0   0  c22 . . .  0    0   c2Np |
       .           | .   .   .   .   .   .  . . .  .    .    .   |
       .           | .   .   .   .   .   .  . . .  .    .    .   |
       .           | .   .   .   .   .   .  . . .  .    .    .   |
       g(e_Pp^p)   |cPp1 0   0  cPp2 0   0  . . . cPpNp 0    0   |
       g(e_Pp^p)   | 0  cPp1 0   0  cPp2 0  . . .  0   cPpNp 0   |
       g(e_Pp^p)   | 0   0  cPp1 0   0  cPp2. . .  0    0   cPpNp|
    '''
    #print(exp_size,contr_size)
    #print(exp_title)
    #print(s)
    Matrix=np.zeros((exp_size*(2*s+1),contr_size*(2*s+1)))
    for i in range(contr_size):
        temp_actv=temp[i].split(',')
        #print(temp_actv)
        S=contraction_normalizer(temp_actv,exp_title,s)
        #print(S)
        first_orb=int(temp_actv[1].split('.')[0])-1
        #print(first_orb)
        for j in range(len(temp_actv)-2):
            #print(i,j,j+2)
            for k in range(2*s+1):
                Matrix[(j+first_orb)*(2*s+1)+k][i*(2*s+1)+k]=float(temp_actv[j+2])/np.sqrt(S)
    return Matrix

def contraction_normalizer(coef,exp,sym):
    ''' Normalize the contracted functions'''
    exp=exp.split(',')
    first_orb=int(coef[1].split('.')[0])-1
    #print(exp)
    #print(coef)
    S=0
    for i in range(2,len(coef)):
        #print('i',coef[i],exp[i+first_orb])
        for j in range(2,len(coef)):
            #print('j',coef[j],exp[j])
            S+=float(coef[i])*float(coef[j])*((4*float(exp[i+first_orb])*float(exp[j+first_orb])/(float(exp[i+first_orb])+float(exp[j+first_orb]))**2)**((2*sym+3)/4))
    #print(S)
    return S
        
