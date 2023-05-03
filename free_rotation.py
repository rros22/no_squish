import math
import numpy as np
import post_processing as post
import update_algorithm as updt
import plot_dynamics as plt_dyn
#define starting orientation with quaternion, and angular velocity
#orientation = np.array()
omega = np.array([0,0,0,np.pi/4])

#VECTOR ARGUMENTS TO FUNCTIONS ARE NP ARRAYS!!

def normalise(vector):
    return vector/np.linalg.norm(vector)


def define_quaternion(angle, axis):
    angle = math.radians(angle)
    axis = normalise(axis)
    return np.hstack((math.cos(angle/2),axis*math.sin(angle/2)))

#converts initial angular velocities to initial quaternion momenta
def quaternion_momenta(m_inertia, omega, q):
    #s matrix
    row_1 = [q[0], -q[1], -q[2], -q[3]]
    row_2 = [q[1], q[0], -q[3], q[2]]
    row_3 = [q[2], q[3], q[0], -q[1]]
    row_4 = [q[3], -q[2], q[1], q[0]]

    S = np.array([row_1, row_2, row_3, row_4])

    #inertia matrix
    row_1 = [m_inertia[0], 0, 0, 0]
    row_2 = [0, m_inertia[1], 0, 0]
    row_3 = [0, 0, m_inertia[2], 0]
    row_4 = [0, 0, 0, m_inertia[3]]
    D = np.array([row_1, row_2, row_3, row_4])

    #convert angular velocities
    omega_r = [0, 0, 0, 0]

    for i in [1, 2, 3]:
        omega_r[i] = math.radians(omega[i])

    return 2*(S @ D @ omega_r)

def quaternion_torque(torque, q):
        #s matrix
        row_1 = [q[0], -q[1], -q[2], -q[3]]
        row_2 = [q[1], q[0], -q[3], q[2]]
        row_3 = [q[2], q[3], q[0], -q[1]]
        row_4 = [q[3], -q[2], q[1], q[0]]
        S = np.array([row_1, row_2, row_3, row_4])

        return 2*(S @ torque)

def rotation_matrix(angle, axis):
    q = define_quaternion(angle, axis)
    row_1 = [q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2, 2*(q[1]*q[2] + q[0]*q[3]), 2*(q[1]*q[3] - q[0]*q[2])]
    row_2 = [2*(q[1]*q[2] - q[0]*q[3]), q[0]**2 - q[1]**2 + q[2]**2 - q[3]**2, 2*(q[2]*q[3] + q[0]*q[1])]
    row_3 = [2*(q[1]*q[3] + q[0]*q[2]), 2*(q[2]*q[3] - q[0]*q[1]), q[0]**2 - q[1]**2 - q[2]**2 + q[3]**2]
    #have to transpose
    return np.transpose(np.array([row_1, row_2, row_3]))

def rotation_matrix_2(q):
    row_1 = [q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2, 2*(q[1]*q[2] + q[0]*q[3]), 2*(q[1]*q[3] - q[0]*q[2])]
    row_2 = [2*(q[1]*q[2] - q[0]*q[3]), q[0]**2 - q[1]**2 + q[2]**2 - q[3]**2, 2*(q[2]*q[3] + q[0]*q[1])]
    row_3 = [2*(q[1]*q[3] + q[0]*q[2]), 2*(q[2]*q[3] - q[0]*q[1]), q[0]**2 - q[1]**2 - q[2]**2 + q[3]**2]
    #have to transpose
    return np.transpose(np.array([row_1, row_2, row_3]))


#interaction site
class interaction_site:

    def __init__(self, l_coord = np.array([0, 0, 0]), g_coord = np.array([0, 0, 0]), force = np.array([0, 0, 0]), mass = 1, no = 0, label = ""):
        self.l_coord = l_coord
        self.g_coord = g_coord
        self.force = force
        self.mass = mass
        self.no = no
        self.label = label

    def print(self):
        print("SITE " + str(self.no) + ":" + "\n")
        print("l_coord : " + str(self.l_coord))
        print("g_coord : " + str(self.g_coord))
        print("force   : " + str(self.force))
        print("mass    : " + str(self.mass))
        print("label   : " + str(self.label))
        print("\n")

class linear_molecule:
    #label & local coordinates
    label = "L2"
    l_1 = np.array([-0.5, 0, 0])
    l_2 = np.array([0.5, 0, 0])
    #interaction sites
    def __init__(self, CoM = np.array([0, 0, 0]), orientation = np.array([0, 0, 1, 0]), omega = np.array([0, 0, 0, 0]), no = 0):
        #molecule parameters
        self.no = no
        self.position= CoM
        self.q = orientation
        self.rot = rotation_matrix_2(self.q)
        #define interaction sites & initialise
        self.sites = [interaction_site(l_coord = self.l_1, no = 2*self.no, label = "N"), interaction_site(l_coord = self.l_2, no = 2*self.no + 1, label = "N")]

        self.sites[0].g_coord = self.position + self.rot @ self.sites[0].l_coord
        self.sites[1].g_coord = self.position + self.rot @ self.sites[1].l_coord

        #inertial parameters
        self.mass = self.sites[0].mass + self.sites[1].mass
        self.m_inertia = np.array([1, 0.001, 1, 1]) #includes ficticious intertia
        self.force = self.sites[0].force + self.sites[1].force
        self.torque = np.array([0, 0, 0, 0])
        self.q_torque = np.array([0, 0, 0, 0])

        #quaternion momenta
        self.p = quaternion_momenta(self.m_inertia, omega, self.q)

    def set_orientation(self, q):
        self.q = q
        self.rot = rotation_matrix_2(self.q)
        #sites
        self.sites[0].g_coord = self.position + self.rot @ self.sites[0].l_coord
        self.sites[1].g_coord = self.position + self.rot @ self.sites[1].l_coord

    def set_torque(self, torque):
        self.torque = torque
        self.q_torque = quaternion_torque(torque, self.q)

    #rotate method // for debug purposes
    def rotate(self, angle, axis):

        new_rot = rotation_matrix(angle, axis)

        self.rot = new_rot @ self.rot

        self.sites[0].g_coord = self.position + self.rot @ self.sites[0].l_coord
        self.sites[1].g_coord = self.position + self.rot @ self.sites[1].l_coord


    #update algorithm (free rotation for now)
    def update_angular_dynamics(self, dt):
        x = [0, 0]# q, p
        #torque update 1/2 (half time-step implicit)
        x[1] = updt.propagator_p(self.p, self.q_torque, dt)
        self.p = x[1]
        #free rotation
        x = updt.propagator_k(3, self.m_inertia, self.p, self.q, dt/2)
        x = updt.propagator_k(2, self.m_inertia, x[1], x[0], dt/2)
        x = updt.propagator_k(1, self.m_inertia, x[1], x[0], dt)
        x = updt.propagator_k(2, self.m_inertia, x[1], x[0], dt/2)
        x = updt.propagator_k(3, self.m_inertia, x[1], x[0], dt/2)

        self.q = x[0]
        self.p = x[1]

        #re-orientate and calculate torques
        self.set_orientation(self.q)
        self.set_torque(self.torque)

        #torque update 1/2 (half time-step implicit)
        x[1] = updt.propagator_p(self.p, self.q_torque, dt) #use new torque
        self.p = x[1]

    def calculate_omega(self):
        #s matrix
        q = self.q
        m_inertia = self.m_inertia

        row_1 = [q[0], -q[1], -q[2], -q[3]]
        row_2 = [q[1], q[0], -q[3], q[2]]
        row_3 = [q[2], q[3], q[0], -q[1]]
        row_4 = [q[3], -q[2], q[1], q[0]]

        S = np.array([row_1, row_2, row_3, row_4])
        S_inv = np.linalg.inv(S)

        #inertia matrix
        row_1 = [m_inertia[0], 0, 0, 0]
        row_2 = [0, m_inertia[1], 0, 0]
        row_3 = [0, 0, m_inertia[2], 0]
        row_4 = [0, 0, 0, m_inertia[3]]
        D = np.array([row_1, row_2, row_3, row_4])
        D_inv = np.linalg.inv(D)

        return 1/2*(D_inv @ S_inv @ self.p)

    def print(self):
        print("MOL " + str(self.no) + ":" + "\n")
        self.sites[0].print()
        self.sites[1].print()

        print("Position   : " + str(self.position))
        print("Orientation: " + str(self.q))
        print("Mass       : " + str(self.mass))
        print("M_Inertia  : " + str(self.m_inertia))
        print("Force      : " + str(self.force))
        print("Quat Torque: " + str(self.q_torque))
        print("\n")


#test code
mol_1 = linear_molecule(CoM = np.array([0, 0, 0]), orientation = np.array([1, 0, 0, 0]), omega = np.array([0, 0, 0, 6]), no = 0)
mol_1.set_torque(np.array([0, 0, 0, 0.01]))


axis = [0, 0, 1]
post.N2_to_pdb(mol_1, "data.pdb")

omega_z = []
z = []

for i in range(1,10000):
    mol_1.update_angular_dynamics(0.01)
    post.N2_to_pdb(mol_1, "data.pdb", "a")
    omega_z.append(math.degrees(mol_1.calculate_omega()[3]))


plt_dyn.plot_omega(omega_z, 0.5)
