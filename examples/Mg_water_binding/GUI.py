#!/Users/Haoyuan/anaconda3/bin/python

'''
GUI for FFEnergy program

Haoyuan Chen
'''

from tkinter import *
from tkinter.filedialog import *
from tkinter.ttk import *
from ase.visualize import view
from ase import Atoms
from Classes import *
import FFEnergy, WaterModelAssign

class MainPage(object):

    def __init__(self, master, info='GUI for FFEnergy program'):
        self.master = master
        self.crd_file = ''
        self.ff_files = ['','','']
        self.pair_file = ''
        self.guests = []
        self.host_atoms = []
        self.guest_atoms = []
        self.all_atoms = []
        self.ff_dict = {}
        self.watermodel = ''
        self.waters = []
        self.watergeom = 0
        self.watergeomval = IntVar(value=self.watergeom)
        self.pairs = []
        self.pairstyle = ''
        self.elj = 0.0
        self.ec = 0.0
        self.et = 0.0
        self.lj_cut = 12.8
        self.lj_scheme = 'T'
        self.lj_mix = 'LB'
        self.coulomb_cut = 12.8
        self.debug = 0
        self.debugval = IntVar(value=self.debug)
        self.row = 0

        self.section1_label = Label(master, text='1. Load Coordinates and Force Field Files', justify=CENTER, foreground='blue', font=('TkDefaultFont',16))
        self.section1_label.grid(row=self.row, columnspan=12)
        self.row += 1 # next row

        self.crd_label = Label(master, text='Coordinates file:')
        self.crd_label.grid(row=self.row, sticky=W)
        self.crd_entry = Entry(master)
        self.crd_entry.grid(row=self.row, column=2, columnspan=8)
        self.crd_entry.insert(0,'.crd')
        self.crd_button = Button(master, text='Choose', command=self.choose_crd_file)
        self.crd_button.grid(row=self.row, column=11)
        self.row += 1 # next row

        self.ff1_label = Label(master, text='Force field file 1 (general, required):')
        self.ff1_label.grid(row=self.row, sticky=W)
        self.ff1_entry = Entry(master)
        self.ff1_entry.grid(row=self.row, column=2, columnspan=8)
        self.ff1_entry.insert(0,'.ff')
        self.ff1_button = Button(master, text='Choose', command=self.choose_ff1_file)
        self.ff1_button.grid(row=self.row, column=11)
        self.row += 1 # next row

        self.ff2_label = Label(master, text='Force field file 2 (specific, optional):')
        self.ff2_label.grid(row=self.row, sticky=W)
        self.ff2_entry = Entry(master)
        self.ff2_entry.grid(row=self.row, column=2, columnspan=8)
        self.ff2_entry.insert(0,'.ff')
        self.ff2_button = Button(master, text='Choose', command=self.choose_ff2_file)
        self.ff2_button.grid(row=self.row, column=11)
        self.row += 1 # next row

        self.ff3_label = Label(master, text='Force field file 3 (specific, optional):')
        self.ff3_label.grid(row=self.row, sticky=W)
        self.ff3_entry = Entry(master)
        self.ff3_entry.grid(row=self.row, column=2, columnspan=8)
        self.ff3_entry.insert(0,'.ff')
        self.ff3_button = Button(master, text='Choose', command=self.choose_ff3_file)
        self.ff3_button.grid(row=self.row, column=11)
        self.row += 1 # next row

        self.load_button = Button(master, text='Load Files', command=self.load_files)
        self.load_button.grid(row=self.row, sticky=W)
        self.view_button = Button(master, text='Visualize Structure', command=self.view_crd)
        self.view_button.grid(row=self.row, column=2)
        self.row += 1 # next row

        self.section2_label = Label(master, text='(Optional) 2. Define Water Molecules and Assign Water Model', justify=CENTER, foreground='green', font=('TkDefaultFont',16))
        self.section2_label.grid(row=self.row, columnspan=12)
        self.row += 1 # next row

        self.water_label = Label(master, text='Water atoms (1-start, separate with , within water and ; between waters, OHH order):')
        self.water_label.grid(row=self.row, sticky=W)
        self.water_entry = Entry(master)
        self.water_entry.grid(row=self.row, column=2, columnspan=8)
        self.water_entry.insert(0,'')
        self.row += 1 # next row

        self.wm_label = Label(master, text='Water model:')
        self.wm_label.grid(row=self.row, sticky=W)
        self.wm_combo = Combobox(master, values=('SPC','SPCE','TIP3P','TIP4P','TIP4PEW'), state="readonly")
        self.wm_combo.grid(row=self.row, column=2)
        self.wm_combo.current(0)
        self.row += 1 # next row

        self.adjust_check = Checkbutton(master, text='Adjust the geometry of waters?', variable=self.watergeomval)
        self.adjust_check.grid(row=self.row, sticky=W)
        self.water_button = Button(master, text='Set Waters', command=self.set_waters)
        self.water_button.grid(row=self.row, column=11)
        self.row += 1 # next row

        self.section3_label = Label(master, text='3. Select Guest Atoms', justify=CENTER, foreground='blue', font=('TkDefaultFont',16))
        self.section3_label.grid(row=self.row, columnspan=12)
        self.row += 1 # next row

        self.guest_label = Label(master, text='Guest atoms (1-start, separate with , , do not forget the M site in 4-pt water model):')
        self.guest_label.grid(row=self.row, sticky=W)
        self.guest_entry = Entry(master)
        self.guest_entry.grid(row=self.row, column=2, columnspan=8)
        self.guest_entry.insert(0,'')
        self.guest_button = Button(master, text='Set Guests', command=self.set_guests)
        self.guest_button.grid(row=self.row, column=11)
        self.row += 1 # next row

        self.section4_label = Label(master, text='(Optional) 4. Add Special Pair Interactions', justify=CENTER, foreground='green', font=('TkDefaultFont',16))
        self.section4_label.grid(row=self.row, columnspan=12)
        self.row += 1 # next row

        self.pair_label = Label(master, text='Atom pairs with special interactions (1-start, separate with , within pair and ; between pairs):')
        self.pair_label.grid(row=self.row, sticky=W)
        self.pair_entry = Entry(master)
        self.pair_entry.grid(row=self.row, column=2, columnspan=8)
        self.pair_entry.insert(0,'')
        self.row += 1 # next row

        self.pf_label = Label(master, text='Pair file:')
        self.pf_label.grid(row=self.row, sticky=W)
        self.pf_entry = Entry(master)
        self.pf_entry.grid(row=self.row, column=2, columnspan=8)
        self.pf_entry.insert(0,'.pair')
        self.pf_button = Button(master, text='Choose', command=self.choose_pair_file)
        self.pf_button.grid(row=self.row, column=11)
        self.row += 1 # next row

        self.ps_label = Label(master, text='Special pair type:')
        self.ps_label.grid(row=self.row, sticky=W)
        self.ps_combo = Combobox(master, values=('LJ1264','Morse'), state="readonly")
        self.ps_combo.grid(row=self.row, column=2)
        self.ps_combo.current(0)
        self.pair_button = Button(master, text='Set Pairs', command=self.set_pairs)
        self.pair_button.grid(row=self.row, column=11)
        self.row += 1 # next row

        self.section5_label = Label(master, text='5. Define Parameters for the Interaction Potentials', justify=CENTER, foreground='blue', font=('TkDefaultFont',16))
        self.section5_label.grid(row=self.row, columnspan=12)
        self.row += 1 # next row

        self.ljmixmethod_label = Label(master, text='LJ mixing scheme:')
        self.ljmixmethod_label.grid(row=self.row, sticky=W)
        self.ljmixmethod_combo = Combobox(master, values=('Lorentz-Berthelot','Waldman-Hagler'), state="readonly")
        self.ljmixmethod_combo.grid(row=self.row, column=2)
        self.ljmixmethod_combo.current(0)
        self.row += 1 # next row

        self.ljcutmethod_label = Label(master, text='LJ cut-off scheme:')
        self.ljcutmethod_label.grid(row=self.row, sticky=W)
        self.ljcutmethod_combo = Combobox(master, values=('Truncated','Shifted'), state="readonly")
        self.ljcutmethod_combo.grid(row=self.row, column=2)
        self.ljcutmethod_combo.current(0)
        self.row += 1 # next row

        self.ljcut_label = Label(master, text='LJ cut-off (Angstrom):')
        self.ljcut_label.grid(row=self.row, sticky=W)
        self.ljcut_entry = Entry(master)
        self.ljcut_entry.grid(row=self.row, column=2, columnspan=8)
        self.ljcut_entry.insert(0,'12.8')
        self.row += 1 # next row

        self.ccut_label = Label(master, text='Coulomb cut-off (Angstrom):')
        self.ccut_label.grid(row=self.row, sticky=W)
        self.ccut_entry = Entry(master)
        self.ccut_entry.grid(row=self.row, column=2, columnspan=8)
        self.ccut_entry.insert(0,'12.8')
        self.param_button = Button(master, text='Set Parameters', command=self.set_params)
        self.param_button.grid(row=self.row, column=11)
        self.row += 1 # next row

        self.section6_label = Label(master, text='6. Calculate the Binding Energy!', justify=CENTER, foreground='red', font=('TkDefaultFont',16))
        self.section6_label.grid(row=self.row, columnspan=12)
        self.row += 1 # next row

        self.run_button = Button(master, text='Get Binding Energy', command=self.run)
        self.run_button.grid(row=self.row, sticky=W)
        self.debug_check = Checkbutton(master, text='Debug (will print details on command line)', variable=self.debugval)
        self.debug_check.grid(row=self.row, column=2)
        self.row += 1 # next row

        self.elj_label = Label(master, text='VDW Energy (kJ/mol):')
        self.elj_label.grid(row=self.row, sticky=W)
        self.elj_entry = Entry(master)
        self.elj_entry.grid(row=self.row, column=2, columnspan=8)
        self.elj_entry.insert(0,'0.0')
        self.row += 1 # next row

        self.ec_label = Label(master, text='Coulomb Energy (kJ/mol):')
        self.ec_label.grid(row=self.row, sticky=W)
        self.ec_entry = Entry(master)
        self.ec_entry.grid(row=self.row, column=2, columnspan=8)
        self.ec_entry.insert(0,'0.0')
        self.row += 1 # next row

        self.et_label = Label(master, text='Total Energy (kJ/mol):')
        self.et_label.grid(row=self.row, sticky=W)
        self.et_entry = Entry(master)
        self.et_entry.grid(row=self.row, column=2, columnspan=8)
        self.et_entry.insert(0,'0.0')
        self.row += 1 # next row

        self.reset_button = Button(master, text='Reset', command=self.reset)
        self.reset_button.grid(row=self.row, sticky=W)
        self.row += 1 # next row

    def choose_crd_file(self):
        self.crd_file = askopenfilename()
        self.crd_entry.delete(0,END)
        self.crd_entry.insert(0,self.crd_file)

    def choose_ff1_file(self):
        self.ff_files[0] = askopenfilename()
        self.ff1_entry.delete(0,END)
        self.ff1_entry.insert(0,self.ff_files[0])

    def choose_ff2_file(self):
        self.ff_files[1] = askopenfilename()
        self.ff2_entry.delete(0,END)
        self.ff2_entry.insert(0,self.ff_files[1])

    def choose_ff3_file(self):
        self.ff_files[2] = askopenfilename()
        self.ff3_entry.delete(0,END)
        self.ff3_entry.insert(0,self.ff_files[2])

    def choose_pair_file(self):
        self.pair_file = askopenfilename()
        self.pf_entry.delete(0,END)
        self.pf_entry.insert(0,self.pair_file)

    def set_guests(self):
        self.guests = []
        self.guest_atoms = []
        self.host_atoms = []
        gs = self.guest_entry.get().strip().split(',')
        for g in gs:
            self.guests.append(int(g))
        f = open(self.crd_file, 'r')
        d = f.readlines()
        f.close()
        for i in range(len(d)):
            if (i+1) in self.guests:
                self.guest_atoms.append(atom(d[i],i+1))
            else:
                self.host_atoms.append(atom(d[i],i+1))

    def set_waters(self):
        self.waters = []
        self.watermodel = self.wm_combo.get()
        ws = self.water_entry.get().strip().split(';')
        for w in ws:
            ww = w.split(',')
            w1 = int(ww[0])
            w2 = int(ww[1])
            w3 = int(ww[2])
            self.waters.append([w1,w2,w3])
        self.watergeom = self.watergeomval.get()

    def set_pairs(self):
        self.pairs = []
        self.pairstyle = self.ps_combo.get()
        ps = self.pair_entry.get().strip().split(';')
        for p in ps:
            pp = p.split(',')
            p1 = int(pp[0])
            p2 = int(pp[1])
            self.pairs.append([p1,p2])

    def set_params(self):
        ljmix = self.ljmixmethod_combo.get()
        if ljmix == 'Lorentz-Berthelot':
            self.lj_mix = 'LB'
        elif ljmix == 'Waldman-Hagler':
            self.lj_mix = 'WH'
        ljsch = self.ljcutmethod_combo.get()
        if ljsch == 'Truncated':
            self.lj_scheme = 'T'
        elif ljsch == 'Shifted':
            self.lj_scheme = 'S'
        self.lj_cut = float(self.ljcut_entry.get())
        self.coulomb_cut = float(self.ccut_entry.get())


    def load_files(self):
        self.crd_parse()
        self.ff_parse()

    def view_crd(self):
        symbol = ''
        pos = []
        for a in self.all_atoms:
            symbol += a.realelement
            posline = (a.x,a.y,a.z)
            pos.append(posline)
        atoms = Atoms(symbol,pos)
        view(atoms)

    def crd_parse(self):
        self.all_atoms = []
        f = open(self.crd_file, 'r')
        d = f.readlines()
        f.close()
        for i in range(len(d)):
            self.all_atoms.append(atom(d[i],i+1))

    def ff_parse(self):
        self.ff_dict = {}
        d = []
        for ff_file in self.ff_files:
            if len(ff_file):
                f = open(ff_file, 'r')
                dd = f.readlines()
                d += dd
                f.close()
        for l in d:
            ll = l.strip().split()
            element = ll[0]
            epsilon = float(ll[1])
            sigma = float(ll[2])
            self.ff_dict[element] = [epsilon, sigma]

    def run(self):
        self.debug = self.debugval.get()
        if self.waters:
            crd = WaterModelAssign.water(self.crd_file,self.watermodel,self.waters,self.watergeom)    
        else:
            crd = self.crd_file
        (self.elj, self.ec, self.et) = FFEnergy.energy(crd,self.ff_files,self.pair_file,self.pairstyle,self.pairs,self.guests,self.lj_cut,self.lj_scheme,self.lj_mix,self.coulomb_cut,self.debug)
        self.elj_entry.delete(0,END)
        self.elj_entry.insert(0,str(self.elj))
        self.ec_entry.delete(0,END)
        self.ec_entry.insert(0,str(self.ec))
        self.et_entry.delete(0,END)
        self.et_entry.insert(0,str(self.et))

    def reset(self):
        self.crd_file = ''
        self.ff_files = ['','','']
        self.pair_file = ''
        self.guests = []
        self.host_atoms = []
        self.guest_atoms = []
        self.all_atoms = []
        self.ff_dict = {}
        self.watermodel = ''
        self.waters = []
        self.watergeomval = IntVar()
        self.watergeom = 0
        self.pairs = []
        self.pairstyle = ''
        self.elj = 0.0
        self.ec = 0.0
        self.et = 0.0
        self.elj_entry.delete(0,END)
        self.ec_entry.delete(0,END)
        self.et_entry.delete(0,END)
        self.crd_entry.delete(0,END)
        self.ff1_entry.delete(0,END)
        self.ff2_entry.delete(0,END)
        self.guest_entry.delete(0,END)
        self.water_entry.delete(0,END)
        self.pair_entry.delete(0,END)
        self.pf_entry.delete(0,END)

if __name__ =='__main__':

    root = Tk()
    root.title("GUI for FFEnergy program")
    gui = MainPage(root)
    root.mainloop()

