from random import randint
import subprocess as sub
from string import Template
from timeit import default_timer as timer

import time
import glob
import sys
import os
import re
import cclib
import itertools

import numpy as np
import pandas as pd


import platform
platform.node()






class DBShotgun(object):
    def __init__(self,  _molecule, 
                        directoryName = 'Input_Files', 
                        functional = ["M06"], 
                        basisSet = ["6-311++G**"]):
        self.molecule = _molecule
        self.directoryName = directoryName
        
        #Create Directory for this molecule
        if not os.path.exists(self.directoryName):
            os.makedirs(self.directoryName)
            self.directoryName = directoryName
            print ('Directory Created: {}'.format(self.directoryName))
        else:
            print ('Directory Exists: {}'.format(self.directoryName))
            
        
        # Set up range of basis set and functional
        self.functional = functional
        self.basisSet = basisSet


    def filename(self):
        return time.strftime("%Y%m%d-%H%M%S")



    def write_input(self, _filename, _func, _basis):
        # Read in input template
    
        with open('input.template','r') as template:
            temp = template.read()
            s = Template(temp)

            # If template file opens write out input files is a molecular subdirectory
            with open( self.directoryName + '/' + _filename + '.inp' , 'w') as f:
                out = s.substitute( Name        = _filename,
                                    Level       = _func,
                                    Basis       = _basis,
                                    Charge      = self.molecule.charge,
                                    Multi       = self.molecule.multiplicity,
                                    XYZ         = self.molecule.get_xyz_string())
                f.write(out)

    def submit_local(self, _filename):
        print(sub.check_call(['pwd']))
        with open( self.directoryName + '/' + _filename + '.out' , 'w') as f:
            return sub.call(['orca', self.directoryName + '/' + _filename + '.inp'], stdout=f)



    def write_submit(self, _filename):
        try:
            os.remove('orca.pbs')
        except:
            print ("No orca file to delete. Continuing.")


        #Write Orca Submit File 
        with open('submit.template' , 'r') as template:
            contents = template.read()
            s = Template(contents)
            substituted = s.substitute( ID = _filename,
                                        DirectoryName =  self.directoryName ,
                                        InputName=  _filename + '.inp',
                                        OutputName =  _filename + '.out')

        with open('orca.pbs', 'w') as writingFile:
            writingFile.write(substituted)

        time.sleep(1)
        return


    def submit_job(self):
        return sub.check_output(['qsub','orca.pbs'])


    def fire(self):
        self.shotgun_graphic()
        current_filename = self.filename()

        for basis, func in itertools.product(self.basisSet, self.functional):
            self.write_input(current_filename, func, basis)
            #if local, no submit needed
            try:
                self.submit_local(current_filename)
            except:
                self.write_submit(current_filename)
                self.submit_job()
            
        return


    def shotgun_graphic(self):
        print(',______________________________________________,,_                              ')
        print('O_________________________,----------._ [____]    \~~~~_____        ______       ')
        print('                       (_(||||||||||||)___________/    __   \______/     /    ')
        print('                           `----------\'       \_\']\__  (__)              |   ')
        print('                                                     \__________         |    ')
        print('                                                                \________\     ')

