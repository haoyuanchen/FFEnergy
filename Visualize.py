#!/Users/Haoyuan/anaconda3/bin/python

'''
convert .crd to .xyz or .gjf (Gaussian input) so it can be visualized (or calculated too)

Haoyuan Chen
'''

from Classes import *

def crd_parse(crd_file):
    f = open(crd_file, 'r')
    d = f.readlines()
    f.close()
    atoms = []
    for i in range(len(d)):
        atoms.append(atom(d[i],i+1))
    return atoms

def xyz_write(crd_file):
    atoms = crd_parse(crd_file)
    f = open(crd_file.replace('.crd','.xyz'), 'w')
    n = len(atoms)
    f.write(str(n)+'\n\n')
    for a in atoms:
        if a.realelement != 'M':
            f.write(a.xyzwrite()+'\n')
    f.close()

def gjf_write(crd_file,chg=0,multi=1,theory='B3LYP',basis='def2SVP',ecp=['',''],dispersion='gd3bj',opt=False,freq=False,scrf='',solvent='',pop='',nproc=1,mem=1,chk=True,additional=''):
    #just basic and most common cases that could be useful in FF-related stuff, not meant to be comprehensive at all
    atoms = crd_parse(crd_file)
    f = open(crd_file.replace('.crd','.gjf'), 'w')
    n = len(atoms)
    f.write('%'+'nproc=%d\n'%nproc)
    f.write('%'+'mem=%dGB\n'%mem)
    if chk:
        f.write('%'+'chk=%s\n'%crd_file.replace('.crd','.chk'))
    mainline = '#p ' + theory + '/'
    if ''.join(ecp):
        mainline += 'gen'
    else:
        mainline += basis
    if dispersion:
        mainline += ' empiricaldispersion='+dispersion
    if opt:
        mainline += ' opt '
    if freq:
        mainline += ' freq '
    if scrf:
        mainline += ' scrf('+scrf+',solvent='+solvent+') '
    if pop:
        mainline += ' pop='+pop
    mainline += ' '+additional
    f.write(mainline+'\n')
    f.write('\nFFEnergy\n\n')
    f.write('%d %d\n'%(chg,multi))
    for a in atoms:
        if a.realelement != 'M':
            f.write(a.xyzwrite()+'\n')
    f.write('\n')
    elemlist = []
    for a in atoms:
        ae = a.realelement
        if ae not in elemlist:
            elemlist.append(ae)
    if ''.join(ecp):
        ecpelems = ecp[0].split()
        regelems = list(set(elemlist)-set(ecpelems))
        if 'M' in regelems:
            regelems.remove('M')
        l1 = ''
        for e in regelems:
            l1 += e
            l1 += ' '
        l1 += ' 0\n'
        f.write(l1)
        f.write(basis+'\n')
        f.write('****\n')
        l2 = ''
        for e in ecpelems:
            l2 += e
            l2 += ' '
        l2 += ' 0\n'
        f.write(l2)
        f.write(ecp[1]+'\n')
        f.write('****\n\n')
        f.write(l2)
        f.write(ecp[1]+'\n')
    f.write('\n')
    f.close()


#single run example
#xyz_write('MOF841_Node_Water_2_Re_HtoHw.crd')
#gjf_write('HKUST1_Node_Methane.crd',chg=0,multi=1,theory='M06L',basis='def2SVP',ecp=['Cu','SDD'],dispersion='gd3',opt=True,freq=False,scrf='',solvent='',pop='',nproc=28,mem=16,chk=True,additional='SCF=XQC')
#gjf_write('MOF841_Node_Water_2_Re_HtoHw_TIP4PEW.crd',chg=0,multi=1,theory='M06L',basis='def2SVP',ecp=['Zr','SDD'],dispersion='gd3',opt=False,freq=False,scrf='',solvent='',pop='',nproc=28,mem=16,chk=True,additional='')
#gjf_write('MgComplex_TIP4PEW.crd',chg=0,multi=1,theory='B3LYP',basis='6-31G*',ecp=['',''],dispersion='gd3bj',opt=True,freq=False,scrf='PCM',solvent='water',pop='',nproc=28,mem=16,chk=True,additional='')
#gjf_write('MgComplex_TIP4PEW_RIGID.crd',chg=0,multi=1,theory='B3LYP',basis='6-31G*',ecp=['',''],dispersion='gd3bj',opt=True,freq=False,scrf='PCM',solvent='water',pop='',nproc=28,mem=16,chk=True,additional='')
