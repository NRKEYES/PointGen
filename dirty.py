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






# For CLI processing
import argparse 
parser = argparse.ArgumentParser(description='Generate a Grid.')
# Set up all allowed arguments.
parser.add_argument('filename', 
                    action='append', 
                    default=[],
                    nargs='?') # allows this posistional argument to be empty

parser.add_argument('-s','--min_seperation', #s for short
                    action='store', 
                    default=1.5,
                    type=int, )

parser.add_argument('-l','--max_seperation', #l for long
                    action='store', 
                    default=20.0,
                    type=int)

parser.add_argument('-i','--iteration_max', 
                    action='store',
                    default=100,
                    type=int)

parser.add_argument('-c','--convergence_limit', 
                    action='store',
                    default=0.1,
                    type=int)                    

# Actually examine CLI input
args = parser.parse_args()
print(args)

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

from string import Template
from tinydb import TinyDB, Query
from tinydb.operations import delete



#Keep database of completed jobs. and Jobs to submit.
# TODO-List
queue = TinyDB('queue.json')
queue_query =  Query()
# DONE-List
completed = TinyDB('completed.json')
completed_query = Query()

# OPEN UP USER DESIGNATED FILE
parser = cclib.io.ccopen(args.filename)
stationary_point_data = parser.parse()

print("There are {} atoms and {} MOs".format(stationary_point_data.natom,
                                        stationary_point_data.nmo))
print("Stationary Point Energy: {}".format(stationary_point_data.scfenergies[-1]))




def cart_to_sphere ( coord ):
    x,y,z = coord[0], coord[1], coord[2]
    
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arctan(y/x)
    phi = np.arctan( np.sqrt(x**2 + y**2) / z )
    return [r,theta,phi]

def minimize_coords( coords ):
    # centering should be done first, hence the subtraction
    # convert all to spherical, centering on Oxygen
    Mn_sc = cart_to_sphere( np.subtract(parsed_data.atomcoords[-1][0] , parsed_data.atomcoords[-1][1] ))
    O_sc  = [0,0,0]
    H1_sc = cart_to_sphere( np.subtract(parsed_data.atomcoords[-1][2] , parsed_data.atomcoords[-1][1] ))
    H2_sc = cart_to_sphere( np.subtract(parsed_data.atomcoords[-1][3] , parsed_data.atomcoords[-1][1] ))
    
    # convert all to minimized
    # Mn Tail only
    Mn_sc = np.subtract( Mn_sc,[0, Mn_sc[1], Mn_sc[2]] )
    H1_sc = np.subtract( H1_sc,[0, Mn_sc[1], Mn_sc[2]] )
    H2_sc = np.subtract( H2_sc,[0, Mn_sc[1], Mn_sc[2]] )
    #flatten H1 onto plane
    H1_sc = np.subtract( H1_sc,[0,0, H1_sc[2]] )
    H2_sc = np.subtract( H2_sc,[0,0, H1_sc[2]] )
    
    this_point = pd.DataFrame({'Energy':parsed_data.scfenergies[-1],
                                                                        
                                'Mnr':  Mn_sc[0],
                               
                               'H1r' :  H1_sc[0],
                               'H1theta' : H1_sc[1],
                               
                               'H2r':   H1_sc[0],
                                'H2theta' : H2_sc[1],
                               'H2phi' : H2_sc[2],
                               
                              }, index=[0])



    for coord in coords:
        pass



#Put points in queue
def add_points_to_queue():
    for stationary_point in args.filename:
        print("\n Currently accesing: {} \n".format(stationary_point))
        
        parser = cclib.io.ccopen(stationary_point )
        stationary_point_data = parser.parse()
        print("There are %i atoms and %i MOs" %(stationary_point_data.natom,
                                                stationary_point_data.nmo))

        print( stationary_point_data.scfenergies[-1])
        # I am ignoring ZPE since most points wont have it available to them.

        for normal_freqs, normal_disp in zip(stationary_point_data.vibfreqs,
                                            stationary_point_data.vibdisps):
            print("-----------Mode-----------")
            print(normal_freqs)
            print("Displacements:\n{}".format(normal_disp))
            print("\nGeometry: \n{}".format(stationary_point_data.atomcoords[-1]))


            #Generate Job for queue
            current_range = np.linspace(-3,3,4) ## change this later to be more clever
            for p in itertools.product(current_range, repeat = 4):
                
                #Check that displacment isn't already in queue

                if  queue.count(queue_query.Displacements == list(p)) > 0     or\
                    completed.count(queue_query.Displacements == list(p)) > 0    :

                    print("Already in queue. Skipping.")
                    continue
                else:            # if actually new submit job  
                    displacement = np.multiply(p,np.transpose(normal_disp)).transpose()
                    
                    new_geom = np.add( displacement, stationary_point_data.atomcoords[-1])


                    print("Adding: \n Displacement Vector: \n{}\n Geometry: \n {}\n".format(p,new_geom)) 

                    queue.insert({  'Displacements': p,
                                    'Geometry': new_geom.tolist(),
                                    'Atoms':    get_element_symbol(stationary_point_data.atomnos )  })
# Submit jobs


def job_watcher():
    runningJobs = 0

    # ensure not too many shotguns are firing at once!
    if len(queue) > 0:
        for point in queue.search(queue_query.Displacements.exists()):
            print(point)


            mol = Molecule.Molecule(point,name = args.filename)

            #if RUNNING < MAX_ RUNNING:

            #Set up the shotgut with the queue database
            shotgun = DBShotgun.DBShotgun(mol)
            #For each geometry in queue call shotgun up until max
            shotgun.fire()

            queue.remove( doc_ids=[point.doc_id])


        
    else:
        #Wait to recheck
        time.sleep(10)
        return









add_points_to_queue()
job_watcher()








