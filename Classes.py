#!/Users/Haoyuan/anaconda3/bin/python

'''
useful classes

Haoyuan Chen
'''

import math

class atom(object):
    def __init__(self, crd_line, crd_idx=1):
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
    def crdwrite(self):
        line = '%s %12.6f %12.6f %12.6f %12.6f'%(self.element,self.x,self.y,self.z,self.charge)
        return line

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

