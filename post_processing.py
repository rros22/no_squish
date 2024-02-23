import math
import decimal
"""
def truncate(number, digits) -> float:
    # Improve accuracy with floating point operations, to avoid truncate(16.4, 2) = 16.39 or truncate(-1.13, 2) = -1.12
    nbDecimals = len(str(number).split('.')[1])
    if nbDecimals <= digits:
        return number
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper
"""
def angstroms(meters):
    return meters*1E10
#takes molecule sites positions and writes them to pdb

def N2_to_pdb(mol, filename):
    #open or create file
    f = open(filename, "a")

    for i in [0,1]:

        #write to file
        f.write("ATOM")
        #calculate character length of atom cardinal number
        no_len = len(str(mol.sites[i].no))
        #write to file atom number
        f.write((11 - 4 - no_len)*" ")
        f.write(str(mol.sites[i].no))
        #calculate name length (S for site + site number i.e one digit for moelcules with less than 10 sites)
        name_len = no_len + 1
        #white spaces for justification
        f.write(2*" ")
        f.write("S" + str(mol.sites[i].no))
        #truncate position decimal places
        x = round(angstroms(mol.sites[i].g_coord[0]), 3)
        y = round(angstroms(mol.sites[i].g_coord[1]), 3)
        z = round(angstroms(mol.sites[i].g_coord[2]), 3)
        #get the length of the coordinate number string
        x_len = len(str(x))
        y_len = len(str(y))
        z_len = len(str(z))
        #white spaces, then coordinates
        f.write((38 - 13 - x_len - name_len)*" ")
        f.write(str(x))
        f.write((46 - 38 - y_len)*" ")
        f.write(str(y))
        f.write((54 - 46 - z_len)*" ")
        f.write(str(z))
        #element symbol
        sym_len = 1
        f.write((78 - 54 - sym_len)*" ")
        f.write(mol.sites[i].label)

        f.write("\n")

def molecules_to_pdb(molecules, filename):
    for mol in molecules:
        N2_to_pdb(mol, filename)

    f = open(filename, "a")
    f.write("ENDMDL")
    f.write("\n")
"""
def n2_to_csv(mol, filename):
    f = open(filename, "a")

    for i in [0,1]:
        f.write("Atom_" str(i))

        x = mol.sites[i].g_coord[0]
        y = mol.sites[i].g_coord[1])
        z = mol.sites[i].g_coord[2])

        f.write(",")
        f.write(str(x) + ',' + str(y) + ',' + str(z) + ',')
"""
