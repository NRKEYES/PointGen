# SPE energies      eV
# freqencies        1/cm
# displacements     delta angstrom

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
from Chem_Tools import *

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
    
    #coords now in:
    # r, theta, phi
    
    
    # convert all to minimized
    # remove Mn extras
    Mn_sc = np.subtract( Mn_sc, [0, Mn_sc[1], Mn_sc[2]] )
    H1_sc = np.subtract( H1_sc, [0, Mn_sc[1], Mn_sc[2]] )
    H2_sc = np.subtract( H2_sc, [0, Mn_sc[1], Mn_sc[2]] )
    #remove H1 Extra--- So H1 is on xz plane
    H1_sc = np.subtract( H1_sc,[0,0, H1_sc[2]] )
    H2_sc = np.subtract( H2_sc,[0,0, H1_sc[2]] )
    
    Mnr, H1r, H1theta, H2r, H2theta, H2phi = Mn_sc[0], H1_sc[0], H1_sc[1], H2_sc[0], H2_sc[1], H2_sc[2]
    
    if np.isnan(H1theta):
        H1theta=0
    if np.isnan(H2theta):
        H2theta=0
  
    return Mnr,H1r,H1theta,H2r,H2theta,H2phi



def distance_check(new_point, old_points):
    # sqrt(     𝑟2+𝑟′2−2𝑟𝑟′[sin(𝜃)sin(𝜃′)cos(𝜙−𝜙′)+cos(𝜃)cos(𝜃′)]    )
    total_distance = new_point.Mnr - old_points

    pass



#Put points in queue
def add_points_to_queue(stationary_point, target_points):
    print("\n Currently accesing: {} \n".format(stationary_point) )
    data = cclib.io.ccopen(stationary_point ).parse()

    print("There are {} atoms and {} MOs".format(data.natom,
                                        data.nmo))

    print("Stationary Point Energy: {}".format(data.scfenergies[-1]))

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
    print(len(orthoganal_modes))

    # Generate Jobs
    current_range = np.linspace(-0.5,1.5,4) # change this later to be more clever

    for p in itertools.product(current_range, repeat=len(orthoganal_modes)):
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
    target_points = add_points_to_queue(stationary_point,target_points)

print(target_points)
print(target_points.describe())
input("Press Enter to continue...")
job_watcher(target_points)








