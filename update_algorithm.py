import math
import numpy as np


#linear dynamics
def lennard_jones_force(site_1, site_2, epsilon, sigma):
    #extract atom positions
    x_1 = site_1.g_coord[0]
    y_1 = site_1.g_coord[1]
    z_1 = site_1.g_coord[2]

    x_2 = site_2.g_coord[0]
    y_2 = site_2.g_coord[1]
    z_2 = site_2.g_coord[2]

    #force vector
    f_x = -4*epsilon*((12*sigma**12*(x_2-x_1)/((x_2-x_1)**2+(y_2-y_1)**2+(z_2-z_1)**2)**7)-(6*sigma**6*(x_2-x_1)/((x_2-x_1)**2+(y_2-y_1)**2+(z_2-z_1)**2)**4))
    f_y = -4*epsilon*((12*sigma**12*(y_2-y_1)/((x_2-x_1)**2+(y_2-y_1)**2+(z_2-z_1)**2)**7)-(6*sigma**6*(y_2-y_1)/((x_2-x_1)**2+(y_2-y_1)**2+(z_2-z_1)**2)**4))
    f_z = -4*epsilon*((12*sigma**12*(z_2-z_1)/((x_2-x_1)**2+(y_2-y_1)**2+(z_2-z_1)**2)**7)-(6*sigma**6*(z_2-z_1)/((x_2-x_1)**2+(y_2-y_1)**2+(z_2-z_1)**2)**4))

    f = np.array([f_x, f_y, f_z])
    f_neg =  np.array([-f_x, -f_y, -f_z])

    #set forces`
    site_1.force = site_1.force + f
    site_2.force = site_2.force + f_neg

def wall_force(site, side_length, epsilon, sigma):
    #define wall_limits
    x_min = -0.1*side_length*2
    x_max = 1.1*side_length*2

    y_min = -0.1*side_length
    y_max = 1.1*side_length

    z_min = -0.1*side_length
    z_max = 1.1*side_length
    #left/right wall force
    site.force[0] = site.force[0] + 48*epsilon*(sigma**12)/(abs(site.g_coord[0] - x_min)**13)
    site.force[0] = site.force[0] - 48*epsilon*(sigma**12)/(abs(x_max - site.g_coord[0])**13)

    #tomp/bottom wall force
    site.force[1] = site.force[1] + 48*epsilon*(sigma**12)/(abs(site.g_coord[1] - y_min)**13)
    site.force[1] = site.force[1] - 48*epsilon*(sigma**12)/(abs(y_max - site.g_coord[1])**13)

    #front/back wall force
    site.force[2] = site.force[2] + 48*epsilon*(sigma**12)/(abs(site.g_coord[2] - z_min)**13)
    site.force[2] = site.force[2] - 48*epsilon*(sigma**12)/(abs(z_max - site.g_coord[2])**13)


def set_forces_sites(molecules):
    #parameters
    epsilon = 6.57E-22
    sigma = 3.17E-10
    #reset fores
    for mol in molecules:
        for site in mol.sites:
            site.reset_force()

    #calculate pairwise forces
    for i in range(0, len(molecules)):
        for j in range(i + 1, len(molecules)):
            lennard_jones_force(molecules[i].sites[0], molecules[j].sites[0], epsilon, sigma)
            lennard_jones_force(molecules[i].sites[0], molecules[j].sites[1], epsilon, sigma)
            lennard_jones_force(molecules[i].sites[1], molecules[j].sites[0], epsilon, sigma)
            lennard_jones_force(molecules[i].sites[1], molecules[j].sites[1], epsilon, sigma)

        #test
        wall_force(molecules[i].sites[0], 4E-9, epsilon, sigma)
        wall_force(molecules[i].sites[1], 4E-9, epsilon, sigma)

def set_CoM_force(molecules):
    for mol in molecules:
        mol.set_CoM_force()

def set_CoM_force_n(molecules):
    for mol in molecules:
        mol.set_CoM_force_n()

def next_position(molecules, dt):

    for mol in molecules:
        mol.position = mol.position + mol.velocity*dt + mol.force*1/2/mol.mass*(dt**2)
def next_velocity(molecules, dt):
    for mol in molecules:
        mol.velocity =  mol.velocity + dt/2/mol.mass*(mol.force + mol.force_n)


def integrate_linear_dynamics(molecules, dt):

    set_forces_sites(molecules)
    set_CoM_force(molecules)
    next_position(molecules, dt)

    #set site site_global_coordinates again to recalculate forces
    for mol in molecules:
        mol.site_global_coordinates()

    set_forces_sites(molecules)
    set_CoM_force_n(molecules)
    next_velocity(molecules, dt)















def global_to_local(q):
    row_1 = [q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2, 2*(q[1]*q[2] + q[0]*q[3]), 2*(q[1]*q[3] - q[0]*q[2])]
    row_2 = [2*(q[1]*q[2] - q[0]*q[3]), q[0]**2 - q[1]**2 + q[2]**2 - q[3]**2, 2*(q[2]*q[3] + q[0]*q[1])]
    row_3 = [2*(q[1]*q[3] + q[0]*q[2]), 2*(q[2]*q[3] - q[0]*q[1]), q[0]**2 - q[1]**2 - q[2]**2 + q[3]**2]
    #have to transpose
    return np.array([row_1, row_2, row_3])












#angular dynamics
def permutations_k(k):
    #p_1
    row_1 = [0, -1, 0, 0]
    row_2 = [1, 0, 0, 0]
    row_3 = [0, 0, 0, 1]
    row_4 = [0, 0, -1, 0]
    p_1 = np.array([row_1, row_2, row_3, row_4])

    #p_2
    row_1 = [0, 0, -1, 0]
    row_2 = [0, 0, 0, -1]
    row_3 = [1, 0, 0, 0]
    row_4 = [0, 1, 0, 0]
    p_2 = np.array([row_1, row_2, row_3, row_4])

    #p_3
    row_1 = [0, 0, 0, -1]
    row_2 = [0, 0, 1, 0]
    row_3 = [0, -1, 0, 0]
    row_4 = [1, 0, 0, 0]
    p_3 = np.array([row_1, row_2, row_3, row_4])

    p = np.array([p_1, p_2, p_3])

    return p[k - 1]


def zeta_k(k, m_inertia, p, q):
    return 1/(4*m_inertia[k])*(np.transpose(p) @ permutations_k(k) @ q)

def propagator_k(k, m_inertia, p, q, dt):
    #external functions
    zeta = zeta_k(k, m_inertia, p, q)
    p_k = permutations_k(k)
    #propagator analytic expansion
    q_new = math.cos(zeta*dt)*q + math.sin(zeta*dt)*(p_k @ q)
    p_new = math.cos(zeta*dt)*p + math.sin(zeta*dt)*(p_k @ p)

    return [q_new, p_new]

def propagator_p(p, q_torque, dt):
    return p + dt/2*q_torque

def compute_torque(molecules):  #just recompute quaternion torque, no force recalculation
    for mol in molecules:
        #convert site forces to local frame of reference
        force_0 = global_to_local(mol.q) @ mol.sites[0].force
        force_1 = global_to_local(mol.q) @ mol.sites[1].force

        #calculate torques
        torque_3 = np.cross(mol.l_1, force_0) + np.cross(mol.l_2, force_1)

        #convert to 4-array
        torque = np.array([0, 0, 0, 0])
        for i in [1, 2, 3]:
            torque[i] = torque_3[i-1]
        #calculate quaternion torque and set new torque
        mol.set_torque(torque)

def update_angular_dynamics(molecules, dt):
    #compute torque
    compute_torque(molecules)

    for mol in molecules:
        x = [0, 0] # q, p
        #torque update 1/2 (half time-step implicit)
        x[1] = propagator_p(mol.p, mol.q_torque, dt)
        mol.p = x[1]
        #free rotation
        x = propagator_k(3, mol.m_inertia, mol.p, mol.q, dt/2)
        x = propagator_k(2, mol.m_inertia, x[1], x[0], dt/2)
        x = propagator_k(1, mol.m_inertia, x[1], x[0], dt)
        x = propagator_k(2, mol.m_inertia, x[1], x[0], dt/2)
        x = propagator_k(3, mol.m_inertia, x[1], x[0], dt/2)

        mol.q = x[0]
        mol.p = x[1]

        #re orientate
        mol.set_orientation(mol.q)

        #calculate new torques
        mol.set_torque(mol.torque)

        #torque update 1/2 (half time-step implicit)
        x[1] = propagator_p(mol.p, mol.q_torque, dt) #use new torque
        mol.p = x[1]
