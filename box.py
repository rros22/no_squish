import math
import numpy as np
import post_processing as post
import linear_molecule as lin
import update_algorithm as updt

class box:

    def __init__(self, size, mol_no):
        self.size = size
        self.mol_no = mol_no
        self.molecules = []

        side = int(np.cbrt(mol_no))
        #create molecules, and initialise them with random positions and orientaions
        for i in range(0, mol_no):
            self.molecules.append(lin.linear_molecule(CoM = np.array([0, 0, 0]), velocity = np.array([0, 0, 0]), orientation = np.array([0, 0, 1, 0]), omega = np.array([0, 0, 0, 0]), no = i))

        #distribute molecules in box
        x = size[0]
        y = size[1]
        z = size[2]

        for i in range(0,side):
            for j in range(0,side):
                for k in range(0, side):
                    self.molecules[i + side*(j + side*k)].position = np.array([x/side + i*x/side, y/side + j*y/side, z/side + k*z/side])
                    #random orientation
                    angle = np.random.uniform(0,360)
                    axis = np.random.uniform(1,10000,3)
                    self.molecules[i + side*(j + side*k)].orientation = lin.define_quaternion(angle, axis)
                    #calculate rotation matrices
                    self.molecules[i + side*(j + side*k)].rot = lin.rotation_matrix_2(self.molecules[i + side*(j + side*k)].orientation)
                    self.molecules[i + side*(j + side*k)].site_global_coordinates()
                    #random velocities
                    self.molecules[i + side*(j + side*k)].velocity = np.random.uniform(-2000,2000,3)
                    #random angular velocities
                    omega = np.random.uniform(-358E12, 358E12, 4)
                    omega[0] = 0
                    self.molecules[i + side*(j + side*k)].p = lin.quaternion_momenta(self.molecules[i + side*(j + side*k)].m_inertia, omega, self.molecules[i + side*(j + side*k)].q)



#delete old data
f = open("data.pdb", "w")
f.close()

domain = box(np.array([1E-9, 1E-9, 1E-9]), 8)
post.molecules_to_pdb(domain.molecules, "data.pdb")

dt = 1E-15

for i in range(0,10000):
    updt.integrate_linear_dynamics(domain.molecules, dt)
    updt.update_angular_dynamics(domain.molecules, dt)

    post.molecules_to_pdb(domain.molecules, "data.pdb")
    print(i)
