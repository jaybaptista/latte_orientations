#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 13:05:58 2018

@author: dv

Edits for general use by R Sanderson 4 March 2020
"""


#%% START
print('====>START<====')


#%% Import packages

import numpy as np
#the following two packages are for reading and working with the snapshot data. 
#pull from https://bitbucket.org/awetzel/ if not using the binder interface (where they are preinstalled)
import gizmo_analysis as ga
import utilities as ut

################# Settings #######################

# Radius of axis ratio calculations (put into logspace) [kpc]
r_ini = 0.1
r_fin = 200

# Number of points calculated for radius in logspace
N_logspace = 50 #1000

# Put radius values into logspace
r_ini_log = np.log10(r_ini)
r_fin_log = np.log10(r_fin)
r = np.logspace(r_ini_log, r_fin_log, N_logspace, base=10)

# Tolerance in axis ratios and maximum iterations (using reduced inertia tensor method)
tolerance = 0.001
iters_max = 1000


# directory where simulations are stored (include a trailing slash)
file_dirc = '../../../data/latte_metaldiff/'

# name of simulation (equals name of directory inside the above where that sim is located)
sims_name = 'm12f_res7100'

#components to calculate the shape of
compts  = ['dark'] # possibilities: ['star','gas','dark','total']


################# Function to calculate axis ratio of a specified component #################
# !!! warning - this function uses global variables defined in the rest of the script.
def axis_ratio_calculation(c, r_max):

    # Set up arrays for x, y, z, m, and d_tot
    if (c == 'total'):
        x = np.hstack( [x_dict[c] for c in ['star','gas','dark']] )
        y = np.hstack( [y_dict[c] for c in ['star','gas','dark']] )
        z = np.hstack( [z_dict[c] for c in ['star','gas','dark']] )
        m = np.hstack( [m_dict[c] for c in ['star','gas','dark']] )
        d_tot = np.hstack( [d_tot_dict[c] for c in ['star','gas','dark']] )
    else:
        x = x_dict[c]
        y = y_dict[c]
        z = z_dict[c]
        m = m_dict[c]
        d_tot = d_tot_dict[c]

    # Initialize s = c/a, p = c/b, and q = b/a
    s = 1
    p = 1
    q = 1

    # Select all particles with total distances less than 2*r_max
    n_sel = (d_tot < 2*r_max)
    x = x[n_sel]
    y = y[n_sel]
    z = z[n_sel]
    m = m[n_sel]
    d = (x**2 + y**2/q**2 + z**2/s**2)**(1/2)

    # Initialize difference values
    ds = 2*tolerance
    dp = 2*tolerance
    dq = 2*tolerance

    # Loop-counter
    loop_counter = 0
    
    n_select = np.ones(4)

    # While-loop to check if
    while not( (ds < tolerance) and (dp < tolerance) and (dq < tolerance)):

        # If-statement for loop-counter
        loop_counter = loop_counter+1
        if loop_counter == iters_max+1:
            loop_counter = loop_counter-1
            break

        # Save current values of axis ratios and ellipsoidal distances
        s_old = s
        p_old = p
        q_old = q
        d_old = d

        # Select particles less than r_max using current ellipsoidal distance
        n_select = (d_old < r_max)
        
        ### JMB
        
        if n_select.sum() == 0:
            print("Selection region contains zero particles, ending iteration...")
            break
        
        ###
        
        
        x_select = x[n_select]
        y_select = y[n_select]
        z_select = z[n_select]
        m_select = m[n_select]

        # Calculate new ellipsoidal distance for this iteration
        d_select = (x_select**2 + y_select**2/q_old**2 + z_select**2/s_old**2)**(1/2)

        # Calculate moment of inertia tensor

        Ixx_terms = m_select*x_select*x_select
        Ixy_terms = m_select*x_select*y_select
        Ixz_terms = m_select*x_select*z_select
        Iyy_terms = m_select*y_select*y_select
        Iyz_terms = m_select*y_select*z_select
        Izz_terms = m_select*z_select*z_select

        Ixx = Ixx_terms.sum()
        Ixy = Ixy_terms.sum()
        Ixz = Ixz_terms.sum()
        Iyy = Iyy_terms.sum()
        Iyz = Iyz_terms.sum()
        Izz = Izz_terms.sum()
        I_m = np.array([[Ixx,Ixy,Ixz],[Ixy,Iyy,Iyz],[Ixz,Iyz,Izz]])/m_select.sum()

        # Eigen values and eigen vectors of reduced inertia tensor matrix
        eigen_values, eigen_vectors = np.linalg.eigh(I_m)

        # Sort eigen vectors to find new x, y, z (where x<->a, y<->b, z<->c)
        #teigen_vectors = np.transpose(eigen_vectors)
        #sorted_vectors = teigen_vectors[indices]
        XX = eigen_vectors[:,2]
        YY = eigen_vectors[:,1]
        ZZ = eigen_vectors[:,0]
        sorted_vectors = [XX, YY, ZZ]
        sorted_vectors = np.transpose(sorted_vectors)

        # check for non-orthogonal eigen vectors
        dot1 = np.around(np.dot(XX,YY), decimals=5)
        dot2 = np.around(np.dot(YY,ZZ), decimals=5)
        dot3 = np.around(np.dot(ZZ,XX), decimals=5)
        if (dot1 == 0) and (dot2 == 0) and (dot3 == 0):
            X = ut.coordinate.get_coordinates_rotated(np.transpose([x,y,z]), sorted_vectors)
            nonorthogonal_counter = 0
        else:
            nonorthogonal_counter = 1
            break

        # Sort eigen values to find new a, b, c (where a>=b>=c)
        #indices = np.unravel_index(np.argsort(eigen_values), eigen_values.shape)
        #sorted_values = eigen_values[indices]
        sorted_values = eigen_values**(1/2)
        a = sorted_values[2]
        b = sorted_values[1]
        c = sorted_values[0]

        #X = ut.coordinate.get_coordinates_rotated(np.transpose([x,y,z]), sorted_vectors)
        x_new = X[:,0]
        y_new = X[:,1]
        z_new = X[:,2]

        # New axis-ratios
        s = c/a
        p = c/b
        q = b/a
        d = (x_new**2 + y_new**2/q**2 + z_new**2/s**2)**(1/2)

        # Differences in axis-ratios
        ds = np.abs(s-s_old)
        dp = np.abs(p-p_old)
        dq = np.abs(q-q_old)

    # Calculate the triaxiality
    T = (a**2-b**2)/(a**2-c**2)

    #and the axis half-lengths (semimajor axis is always r_max)
    b = r_max * q
    c = r_max * s

    #determine the number of particles used
    n = n_select.sum()

    # Store radius, # of particles, axis ratios & lengths, triaxiality, number of iters, number of non-orthogonal systems, & eigenvectors in an array
    # note that eigenvectors are defined in the default principal-axis coordinate system set up during snapshot loading
    axis_ratios_comp = np.array([r_max, n, s, p, q, b, c,
                                 T, loop_counter, nonorthogonal_counter,
                                 sorted_vectors[0,0], sorted_vectors[0,1], sorted_vectors[0,2],
                                 sorted_vectors[1,0], sorted_vectors[1,1], sorted_vectors[1,2],
                                 sorted_vectors[2,0], sorted_vectors[2,1], sorted_vectors[2,2]])
    #print(axis_ratios_comp)

    return axis_ratios_comp



######### start of main script ########
#%% Read snapshot and calculate the positions relative to center, in a first estimate of the principal axis frame

# Open snapshot file
if 'total' in compts:
    compts_read =  ['star','gas','dark']
else:
    compts_read = compts

part = ga.io.Read.read_snapshots(compts_read, 'snapshot', 600, simulation_directory=file_dirc+sims_name, assign_hosts_rotation=True)


# x,y,z of particles from components - pre-calculate these to save time
x_dict={}
y_dict={}
z_dict={}
m_dict={}
d_tot_dict = {}

if 'total' in compts:
    for c in ['star','dark','gas']:
        dsel = (part[c].prop('host.distance.total')<2*r_fin)
        d_tot_dict[c] = part[c].prop('host.distance.total')[dsel]
        x_dict[c] = part[c].prop('host.distance.principal')[dsel,0]
        y_dict[c] = part[c].prop('host.distance.principal')[dsel,1]
        z_dict[c] = part[c].prop('host.distance.principal')[dsel,2]
        m_dict[c] = part[c]['mass'][dsel]    
else:
    for c in compts:
        dsel = (part[c].prop('host.distance.total')<2*r_fin)
        d_tot_dict[c] = part[c].prop('host.distance.total')[dsel]
        x_dict[c] = part[c].prop('host.distance.principal')[dsel,0]
        y_dict[c] = part[c].prop('host.distance.principal')[dsel,1]
        z_dict[c] = part[c].prop('host.distance.principal')[dsel,2]
        m_dict[c] = part[c]['mass'][dsel]


# release memory from the snapshot
del(part)

#%% Calculate axis ratios

# Initialize for-loop
axis_ratios = {}
for c in compts:
    print('calculating axis ratios for {0}'.format(c),flush=True)
    axis_ratios[c] = np.zeros([len(r), 19])  #look in the function defined above to see what is stored in this array
    # For-loop over all r where you want to calculate the shape
    for i,ri in enumerate(r):
        print('{0} kpc'.format(ri),flush=True) #this can take a while for large r, so print status updates
        axis_ratios[c][i,:] = axis_ratio_calculation(c, ri)
    print('Saving to file',flush=True)
    with open('file_axis_ratios_{0}_{1}_EVEN_MORE_WRONG.txt'.format(sims_name,c),'w') as f:
        f.write('#r n s p q b c T nloop nnorth e00 e01 e02 e10 e11 e12 e20 e21 e22 \n')
        np.savetxt(f, axis_ratios[c])


#%% END
print('====>END<====')