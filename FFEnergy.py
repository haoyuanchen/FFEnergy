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

example force field file (epsilon in K, sigma in A) (no empty lines in the end) (can have comments after the numbers):
H           22.14    2.57
C           52.83    3.43
O           30.19    3.12
Cu          2.52     3.11

example LJ1264 pair file (C12, C6, C4 (all in units derived from K and A)) (no empty lines in the end) (can have comments after the numbers) (order in pair doesn't matter):
Mg    OwTIP4PEW    36605101  66677  90833  #http://ambermd.org/tutorials/advanced/tutorial20/12_6_4.htm#ref3

example Morse pair file (De, re, a (all in units derived from K and A)) (no empty lines in the end) (can have comments after the numbers) (order in pair doesn't matter):
Mg    OwTIP4PEW    1000  2.0  1.0  #made up
'''

from Classes import *

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
 
def ff_parse(ff_files):
    d = []
    for ff_file in ff_files:
        f = open(ff_file, 'r')
        dd = f.readlines()
        d += dd
        f.close()
    ff_dict = {}
    for l in d:
        ll = l.strip().split()
        element = ll[0]
        epsilon = float(ll[1])
        sigma = float(ll[2])
        ff_dict[element] = [epsilon, sigma]
    return ff_dict

def pair_parse(pair_file,pair_style):
    f = open(pair_file, 'r')
    d = f.readlines()
    f.close()
    pair_dict = {}
    for l in d:
        ll = l.strip().split()
        element1 = ll[0]
        element2 = ll[1]
        if pair_style == 'LJ1264':
            c12 = float(ll[2])
            c6 = float(ll[3])
            c4 = float(ll[4])
            pair_dict[(element1,element2)] = [c12, c6, c4]
            pair_dict[(element2,element1)] = [c12, c6, c4]
        elif pair_style == 'Morse':
            de = float(ll[2])
            re = float(ll[3])
            a = float(ll[4])
            pair_dict[(element1,element2)] = [de, re, a]
            pair_dict[(element2,element1)] = [de, re, a]
    return pair_dict

def lj(atom1,atom2,cut=12.8,scheme='T',mix='LB',debug=False):
    r = pair(atom1,atom2).r
    if mix == 'LB':  #Lorentz-Berthelot
        epsilon = np.sqrt(atom1.epsilon*atom2.epsilon)
        sigma = (atom1.sigma+atom2.sigma)*0.5
    elif mix == 'WH':  #Waldman-Hagler
        epsilon = 2.0*np.sqrt(atom1.epsilon*atom2.epsilon)*(atom1.sigma*atom2.sigma)**3/((atom1.sigma)**6+(atom2.sigma)**6)
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

def pair_lj1264(pair1,cut=12.8,scheme='T',debug=False):  
    r = pair1.r
    c12 = pair1.c12
    c6 = pair1.c6
    c4 = pair1.c4
    e = 0.0
    if scheme == 'T':  #truncated directly at cutoff
        if r >= cut:
            e = 0.0
        else:
            e = c12/(r**12) - c6/(r**6) - c4/(r**4)
    elif scheme == 'S':  #shifted vertically to have 0 value at cutoff
        if r >= cut:
            e = 0.0
        else:
            e_end = c12/(cut**12) - c6/(cut**6) - c4/(cut**4)
            e = c12/(r**12) - c6/(r**6) - c4/(r**4) - e_end
    #elif scheme == 'TTC':  #truncated with tail corrections
    if debug:
        print(pair1.element1,pair1.element2,r,'A,',e*0.0083,'kJ/mol,LJ1264')    
    return e*0.0083    #convert to kJ/mol

def pair_morse(pair1,cut=12.8,scheme='T',debug=False):
    r = pair1.r
    de = pair1.de
    re = pair1.re
    a = pair1.a
    e = 0.0
    if scheme == 'T':  #truncated directly at cutoff
        if r >= cut:
            e = 0.0
        else:
            e = de*(1-np.exp(a*(re-r)))**2
    elif scheme == 'S':  #shifted vertically to have 0 value at cutoff
        if r >= cut:
            e = 0.0
        else:
            e_end = de*(1-np.exp(a*(re-cut)))**2
            e = de*(1-np.exp(a*(re-r)))**2 - e_end
    #elif scheme == 'TTC':  #truncated with tail corrections
    if debug:
        print(pair1.element1,pair1.element2,r,'A,',e*0.0083,'kJ/mol,Morse')
    return e*0.0083    #convert to kJ/mol

def coulomb(atom1,atom2,cut=12.8,debug=False):
    r = pair(atom1,atom2).r
    e = 0.0
    if r >= cut:
        e = 0.0
    else:
        e = atom1.charge*atom2.charge/r
    if debug:
        print(atom1.element,atom2.element,r,'A,',e*1389,'kJ/mol,Coulomb')   
    return e*1389    #convert to kJ/mol

def energy(crd_file,ff_files,pair_file,pair_style,pairs,guests,lj_cut=12.8,lj_scheme='T',lj_mix='LB',coulomb_cut=12.8,debug=False):
    host_atoms,guest_atoms = crd_parse(crd_file,guests)
    ff_dict = ff_parse(ff_files)
    #print(ff_dict)    #for testing
    if pair_file and pair_style and pairs:
        pair_dict = pair_parse(pair_file,pair_style)
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
                p.getff(pair_dict,pair_style)
                if pair_style == 'LJ1264':
                    e_lj += pair_lj1264(p,lj_cut,lj_scheme,debug)
                elif pair_style == 'Morse':
                    e_lj += pair_morse(p,lj_cut,lj_scheme,debug)
            else:
                e_lj += lj(h,g,lj_cut,lj_scheme,lj_mix,debug)
            e_coulomb += coulomb(h,g,coulomb_cut,debug)
    if debug:
        print('E_LJ = %16.4f kJ/mol\nE_Coulomb = %16.4f kJ/mol\nE_Total = %16.4f kJ/mol'%(e_lj,e_coulomb,e_lj+e_coulomb))
    return e_lj, e_coulomb, e_lj+e_coulomb

#single run example
#energy('MOF841_Node_Water_2_Re_SPC.crd',['UFF.ff','Water.ff'],'','',[],[71,72,73],debug=True)
#energy('MOF841_Node_Water_2_Re_HtoHw_SPC.crd',['UFF.ff','Water.ff'],'','',[],[71,72,73],debug=True)
#energy('MOF841_Node_Water_2_Re_HtoHw_TIP4PEW.crd',['UFF.ff','Water.ff'],'','',[],[71,72,73,74],debug=True)
#energy('MgComplex_TIP4PEW.crd',['UFF.ff','Water.ff'],'LJ1264.pair','LJ1264',[[1,13]],[13,14,15,16],debug=True)
#energy('MgComplex_TIP4PEW.crd',['UFF.ff','Water.ff'],'Morse.pair','Morse',[[1,13]],[13,14,15,16],debug=True)
