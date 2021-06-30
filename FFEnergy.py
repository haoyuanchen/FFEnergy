#!/Users/Haoyuan/anaconda3/bin/python

'''
calculates binding energy in host-guest complexes using force fields

Haoyuan Chen

example coordinate file (X/Y/Z in A, charge in e) (no empty lines in the end):
C            8.8379    4.3797   23.6832  -0.10
C            9.8205    5.0770   22.9870  -0.10
H           10.2176    6.0204   23.3654   0.20
Cu          13.1715    7.4586   20.6275   1.00
Cu          13.1715    5.7128   18.8869   1.00
O           11.7809    8.3626   19.5668  -0.50
O           14.5620    6.4018   21.5362  -0.50
O           14.5620    8.3626   19.5668  -0.50
O           11.7809    6.4018   21.5362  -0.50
O           11.7809    4.8088   19.9476  -0.50

example force field file (epsilon in K, sigma in A) (no empty lines in the end):
H           22.14    2.57
C           52.83    3.43
O           30.19    3.12
Cu          2.52     3.11
Kr          166.4    3.636
Xe          221.0    4.10

example pair file (C12, C6, C4 (all in units derived from K and A)) (no empty lines in the end) (order in pair doesn't matter):
Mg    Ow    36605101  66677  90833
'''

import math

class atom(object):
    def __init__(self, crd_line, crd_idx):
        self.idx = crd_idx  #1-start
        self.element = crd_line.strip().split()[0]
        self.x = float(crd_line.strip().split()[1])
        self.y = float(crd_line.strip().split()[2])
        self.z = float(crd_line.strip().split()[3])
        self.charge = float(crd_line.strip().split()[4])
        self.epsilon = 0.0
        self.sigma = 0.0
    def getff(self, ff_dict):
        self.epsilon = ff_dict[self.element][0]
        self.sigma = ff_dict[self.element][1]

class pair(object):
    def __init__(self, atom1, atom2):
        self.element1 = atom1.element
        self.element2 = atom2.element
        self.x1 = atom1.x
        self.x2 = atom2.x
        self.y1 = atom1.y
        self.y2 = atom2.y
        self.z1 = atom1.z
        self.z2 = atom2.z
        self.r = math.sqrt((self.x1-self.x2)**2+(self.y1-self.y2)**2+(self.z1-self.z2)**2)
        self.c12 = 0.0
        self.c6 = 0.0
        self.c4 = 0.0
    def getff(self, pair_dict):
        self.c12 = pair_dict[self.element1,self.element2][0]
        self.c6 = pair_dict[self.element1,self.element2][1]
        self.c4 = pair_dict[self.element1,self.element2][2]

def crd_parse(crd_file,guests):
    #guests: a list of the indices of guest atoms, 1-start
    f = open(crd_file, 'r')
    d = f.readlines()
    f.close()
    host_atoms = []
    guest_atoms = []
    for i in range(len(d)):
       if (i+1) in guests:
           guest_atoms.append(atom(d[i],i+1))
       else:
           host_atoms.append(atom(d[i],i+1))
    return host_atoms,guest_atoms
 
def ff_parse(ff_file):
    f = open(ff_file, 'r')
    d = f.readlines()
    f.close()
    ff_dict = {}
    for l in d:
        ll = l.strip().split()
        element = ll[0]
        epsilon = float(ll[1])
        sigma = float(ll[2])
        ff_dict[element] = [epsilon, sigma]
    return ff_dict

def pair_parse(pair_file):
    f = open(pair_file, 'r')
    d = f.readlines()
    f.close()
    pair_dict = {}
    for l in d:
        ll = l.strip().split()
        element1 = ll[0]
        element2 = ll[1]
        c12 = float(ll[2])
        c6 = float(ll[3])
        c4 = float(ll[4])
        pair_dict[(element1,element2)] = [c12, c6, c4]
        pair_dict[(element2,element1)] = [c12, c6, c4]
    return pair_dict

def lj(atom1,atom2,cut=12.8,scheme='T',mix='LB',debug=False):
    r = math.sqrt((atom1.x-atom2.x)**2+(atom1.y-atom2.y)**2+(atom1.z-atom2.z)**2)
    if mix == 'LB':  #Lorentz-Berthelot
        epsilon = math.sqrt(atom1.epsilon*atom2.epsilon)
        sigma = (atom1.sigma+atom2.sigma)*0.5
    elif mix == 'WH':  #Waldman-Hagler
        epsilon = 2.0*math.sqrt(atom1.epsilon*atom2.epsilon)*(atom1.sigma*atom2.sigma)**3/((atom1.sigma)**6+(atom2.sigma)**6)
        sigma = (((atom1.sigma)**6+(atom2.sigma)**6)/2.0)**(1/6)
    #elif mix == 'TT':  #Tang-Toennies
        #epsilon = 
        #sigma =
    e = 0.0
    if scheme == 'T':  #truncated directly at cutoff
        if r >= cut:
            e = 0.0
        else:
            sixth = (sigma/r)**6
            e = 4*epsilon*(sixth**2-sixth)
    elif scheme == 'S':  #shifted vertically to have 0 value at cutoff
        if r >= cut:
            e = 0.0
        else:
            sixth_end = (sigma/cut)**6
            e_end = 4*epsilon*(sixth_end**2-sixth_end)
            sixth = (sigma/r)**6
            e = 4*epsilon*(sixth**2-sixth) - e_end
    #elif scheme == 'TTC':  #truncated with tail corrections
    if debug:
        print(atom1.element,atom2.element,r,'A,',e*0.0083,'kJ/mol,LJ')    
    return e*0.0083    #convert to kJ/mol

def lj_pair(pair1,cut=12.8,scheme='T',debug=False):  #specifically-defined pairwise interactions (like 12-6-4) that will override regular LJ
    r = pair1.r
    c12 = pair1.c12
    c6 = pair1.c6
    c4 = pair1.c4
    e = c12/(r**12) - c6/(r**6) - c4/(r**4)
    if debug:
        print(pair1.element1,pair1.element2,r,'A,',e*0.0083,'kJ/mol,LJ_Pair')    
    return e*0.0083    #convert to kJ/mol

def coulomb(atom1,atom2,cut=12.8,debug=False):
    r = math.sqrt((atom1.x-atom2.x)**2+(atom1.y-atom2.y)**2+(atom1.z-atom2.z)**2)
    e = 0.0
    if r >= cut:
        e = 0.0
    else:
        e = atom1.charge*atom2.charge/r
    if debug:
        print(atom1.element,atom2.element,r,'A,',e*1389,'kJ/mol,Coulomb')   
    return e*1389    #convert to kJ/mol

def energy(crd_file,ff_file,pair_file,pairs,guests,lj_cut=12.8,lj_scheme='T',lj_mix='LB',coulomb_cut=12.8,debug=False):
    host_atoms,guest_atoms = crd_parse(crd_file,guests)
    ff_dict = ff_parse(ff_file)
    #print(ff_dict)    #for testing
    if pair_file:
        pair_dict = pair_parse(pair_file)
    #print(pair_dict)    #for testing
    for h in host_atoms:
        h.getff(ff_dict)
    for g in guest_atoms:
        g.getff(ff_dict)
    e_lj = 0.0
    e_coulomb = 0.0
    for h in host_atoms:
        for g in guest_atoms:
            if ([h.idx,g.idx] in pairs) or ([g.idx,h.idx] in pairs):  #not all pairs of atoms with defined specific pairwise interactions will be treated that way, only the ones in 'pairs' will
                p = pair(h,g)
                p.getff(pair_dict)
                e_lj += lj_pair(p,lj_cut,lj_scheme,debug)
            else:
                e_lj += lj(h,g,lj_cut,lj_scheme,lj_mix,debug)
            e_coulomb += coulomb(h,g,coulomb_cut,debug)
    if debug:
        print('E_LJ = %16.4f kJ/mol\nE_Coulomb = %16.4f kJ/mol\nE_Total = %16.4f kJ/mol'%(e_lj,e_coulomb,e_lj+e_coulomb))
    return e_lj, e_coulomb, e_lj+e_coulomb

#single run example
#energy('SPC_dimer.crd','UFF_And_SPC.ff','',[],[1,2,3],debug=True)
#energy('MOF841_Node_Water_2_Re.crd','UFF_And_SPC.ff','',[],[71,72,73],debug=True)
#energy('MOF841_Node_Water_2_Re_HtoHw.crd','UFF_And_SPC.ff','',[],[71,72,73],debug=True)
energy('MgComplex.crd','UFF_And_SPC.ff','LJ1264.pair',[[1,13]],[13,14,15],debug=True)
