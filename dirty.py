# SPE energies      eV
# freqencies        1/cm
# displacements     delta angstrom


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



#######################################################
#    Main Code                                        #
#######################################################
import cclib
import glob
import re
import os
import Molecule
import DBShotgun
import itertools
import time

import numpy as np
import pandas as pd




# TODO generalize this function
def minimize_coords( coords ):
    # centering should be done first, hence the subtraction
    # convert all to spherical, centering on Oxygen
    Mn_sc = cart_to_sphere( np.subtract(coords[0] , coords[1] ))
    O_sc  = [0,0,0]
    H1_sc = cart_to_sphere( np.subtract(coords[2] , coords[1] ))
    H2_sc = cart_to_sphere( np.subtract(coords[3] , coords[1] ))
    
    # convert all to minimized
    # Mn Tail only
    Mn_sc = np.subtract( Mn_sc,[0, Mn_sc[1], Mn_sc[2]] )
    H1_sc = np.subtract( H1_sc,[0, Mn_sc[1], Mn_sc[2]] )
    H2_sc = np.subtract( H2_sc,[0, Mn_sc[1], Mn_sc[2]] )
    #flatten H1 onto plane
    H1_sc = np.subtract( H1_sc,[0,0, H1_sc[2]] )
    H2_sc = np.subtract( H2_sc,[0,0, H1_sc[2]] )
    
    Mnr,H1r,H1theta,H2r,H2theta,H2phi = Mn_sc[0], H1_sc[0], H1_sc[1], H2_sc[0], H2_sc[1], H2_sc[2]

    return Mnr,H1r,H1theta,H2r,H2theta,H2phi




#Put points in queue
def add_points_to_queue(stationary_point, target_points):
    print("\n Currently accesing: {} \n".format(stationary_point) )
    data = cclib.io.ccopen(stationary_point ).parse()

    # Get the well's eq geom to apply the displacements to
    corrected_coords = minimize_coords(data.atomcoords[-1] )
    print("\nEquilibrium Geometry: \n{}".format( corrected_coords ) )

    orthoganal_modes = []
    for normal_freq, normal_disp in zip(data.vibfreqs, data.vibdisps):
        # this puts the displacements in terms of internal coordinates
        corrected_disp = minimize_coords( normal_disp )     
        print("-----------Mode-----------")
        print(normal_freq)
        print("Displacements:\n{}".format( corrected_disp ) )
        orthoganal_modes.append(corrected_disp)

    # Generate Jobs
    current_range = np.linspace(0,1,2) # change this later to be more clever

    for p in itertools.permutations(current_range, r=len(orthoganal_modes)):
        # combine orthogonal modes with eachother and a magnitude factor

        displacement = np.zeros_like(orthoganal_modes[0])
        for number,mode in enumerate(orthoganal_modes):
            displacement += np.multiply(orthoganal_modes[number],p[number])

        displaced_coords = np.add( displacement, corrected_coords)

        print("Adding: \n Displacement Vector: \n{}\n Geometry: \n {}\n".format(p,displaced_coords)) 


        this_point = pd.DataFrame({   
                                'Mnr'       :   displaced_coords[0],
                               
                               'H1r'        :   displaced_coords[1],
                               'H1theta'    :   displaced_coords[2],
                               
                                'H2r'       :   displaced_coords[3],
                                'H2theta'   :   displaced_coords[4],
                                'H2phi'     :   displaced_coords[5]
                              }, index=[0])


        target_points = target_points.append(this_point)
    
    return target_points
    

# generalize this function at some point
def sphere_to_cart( coord ):
    print(coord)
    # x=rsinϕcosθ
    # y=rsinϕsinθ
    # z=rcosϕ
    return  [   [coord[0]*np.sin(0)*np.cos(0), coord[0]*np.sin(0)*np.sin(0),  coord[0]*np.cos(0)],
                [0,0,0],
                [coord[1]*np.sin(0)*np.cos(coord[2]), coord[1]*np.sin(0)*np.sin(coord[2]),  coord[1]*np.cos(2)],
                [coord[3]*np.sin(coord[5])*np.cos(coord[4]), coord[3]*np.sin(coord[5])*np.sin(coord[4]),  coord[3]*np.cos(coord[5])]
            ]


def job_watcher(target_points):
    print("Submitting Jobs.")

    for key, point in target_points.iterrows():
        print(key)
        mol = Molecule.Molecule({'Geometry' :sphere_to_cart(point),
                                'Atoms'     :['Mn','O','H','H']},name = 'SPE Calc')
        #if RUNNING < MAX_ RUNNING:
        #Set up the shotgut with point
        shotgun = DBShotgun.DBShotgun(mol)

        #For each geometry in queue call shotgun up until max
        shotgun.fire()



column_list = ['Mnr',
               'H1r','H1theta',
               'H2r','H2theta','H2phi']

target_points = pd.DataFrame(columns=column_list)
print(target_points)

for stationary_point in glob.glob("*.out"):
    print(stationary_point)
    parser = cclib.io.ccopen(stationary_point)
    stationary_point_data = parser.parse()

    print("There are {} atoms and {} MOs".format(stationary_point_data.natom,
                                        stationary_point_data.nmo))

    print("Stationary Point Energy: {}".format(stationary_point_data.scfenergies[-1]))

    target_points = add_points_to_queue(stationary_point,target_points)

print(target_points)
job_watcher(target_points)








