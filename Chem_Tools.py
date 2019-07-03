import re
import numpy as np

# These are almost becoming universal functions I need to create full standalone scripts for.
#Currently Unused for Zero Point Energy Correction.
def get_ZPE(filename):
    try:
        with open(filename, 'r') as file:
            file = file.read()
            ZPE = re.findall(r'Non-thermal .ZPE. correction              (.*) Eh', file) #Pull ZPE values in Hartree
            if not ZPE:
                return 0
            else:
                ZPE = float(ZPE[-1])
                return ZPE*27.212 # Convent ZPE to eV
    except:
        print ("Can not open file for ZPE")

def get_element_symbol(element_number):
    symbol = ('H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl',
              'Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As',
              'Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',
              'Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb',
              'Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl',
              'Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk',
              'Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh',
              'Fl','Mc','Lv','Ts','Og')
    if type(element_number) == int:
        return str(symbol[element_number-1])

    symbols =[]
    for num in element_number:
        symbols.append(str(symbol[num-1]))          
    return symbols




def cart_to_sphere ( coord ):
    x,y,z = coord[0], coord[1], coord[2]
    
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arctan(y/x)
    phi = np.arctan( np.sqrt(x**2 + y**2) / z )
    return [r,theta,phi]

def freq_in(freq, mass1, mass2):
    # freq must be in wavenumbers
     
    dynes_per_cm = ((freq/4.12)**2)*(mass1*mass2/(mass1+mass2))
    
    mdynes_per_ang = dynes_per_cm*1000/(100000000)
    return mdynes_per_ang