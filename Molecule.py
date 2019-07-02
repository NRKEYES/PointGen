class Molecule(object):
    def __init__(self, _data, _charge = 1, _spin = 5, name = ''):
        
        self.name = name
        self.charge = _charge
        self.multiplicity = _spin
        self.atoms = _data['Atoms']
        
        self.xyz = _data['Geometry']  # This is just lines of text.
        
    def get_xyz_string(self):
        print(self.xyz)
        print(self.atoms)

        xyzString = ""
        for i, coords in enumerate(self.xyz):
            xyzString = "{} {} {}\n".format(xyzString,  str(self.atoms[i]), (' '.join(str(i) for i in coords)))
        return xyzString
        
    def molecule_print(self):
        print ("    Name: " + self.name)
        print ("    Charge: " + str(self.charge))
        print ("    Multiplicity: " + str(self.multiplicity))
                
        for atom in self.xyz:
            print (u"\u269B  " +atom)
