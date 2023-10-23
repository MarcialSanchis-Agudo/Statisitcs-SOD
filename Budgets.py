#%%

import numpy as np
import h5py as h5
import pandas as pd
import matplotlib as plt

from mpi4py import MPI
from pandas import HDFStore
import matplotlib.pyplot as plt



comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
filename ="/home/marcial/Simulacionprueba/rank_4/channel_equi/results_AVGRSS_channel-4_21.h5"
#filename1 = "Re180.078.h5"
dataset = {}
periodic = 1
tolerance = 1e-8

print("Pre-postproc by rank 0, reading hdf5 file {}".format(rank))
with h5.File(filename, "r") as f:
    # Print all root level object names (aka keys) 
    # these can be group or dataset names 
        
    print("Keys: %s" % f.keys())
    column_names = list(f.keys())


    # Create a datafnrame with column names as keys and empty lists as values
    stats = pd.DataFrame()
    
        
    # Read the data from the HDF5 file and populate the dataset
    for col_name in column_names:
        dataset[col_name] = f[col_name][:]
        stats[col_name] = f[col_name][:]
        
        #print(stats['avvel_x'])
    counts = stats['x'].value_counts()
    print(counts)

    #stats['y'] = stats['y'].astype(float)
    #stats['x'] = stats['x'].astype(float)
    #print(stats['avvel_x'])


        
#%% 
    # Reynolds stresses
stats['avR_u2'] = stats['avve2_x'] - (stats['avvel_x'] * stats['avvel_x'])
stats['avR_uv'] = stats['avvex_x'] - stats['avvel_x'] * stats['avvel_y']
stats['avR_uw'] = stats['avvex_y'] - stats['avvel_x'] * stats['avvel_z']
stats['avR_v2'] = stats['avve2_y'] - stats['avvel_y'] * stats['avvel_y']
stats['avR_vw'] = stats['avvex_z'] - stats['avvel_y'] * stats['avvel_z']
stats['avR_w2'] = stats['avve2_z'] - stats['avvel_z'] * stats['avvel_z']
stats['avvel2'] = stats['avvel_x'] * stats['avvel_x']


#%%
"""    
    # Production term
        def Production_budget_2D(stats):
            P_xx = -2*(stats['avR_u2'] * stats['gradU_xpress_mean_x'] + stats['avR_uv']*stats['gradU_xpress_mean_y'])
            P_yy = -2*(stats['avR_uv'] * stats['gradU_ypress_mean_x'] + stats['avR_v2']*stats['gradU_ypress_mean_y'])
            P_zz = -2*(stats['avR_uw']*stats['gradU_zpress_mean_x'] + stats['avR_vw']*stats['gradU_zpress_mean_y'])
            P_xy = -(stats['avR_u2']*stats['gradU_ypress_mean_x'] + stats['avR_uv']*stats['gradU_xpress_mean_x'] + stats['avR_uv']*stats['gradU_ypress_mean_y'] + stats['avR_v2']*stats['gradU_xpress_mean_y'])
            P_xz = -(stats['avR_u2']*stats['gradU_zpress_mean_x'] + stats['avR_uw']*stats['gradU_xpress_mean_x'] + stats['avR_uv']*stats['gradU_zpress_mean_y'] + stats['avR_vw']*stats['gradU_xpress_mean_y'])
            P_yz = -(stats['avR_uv']*stats['gradU_zpress_mean_x'] + stats['avR_uw']*stats['gradU_ypress_mean_x'] + stats['avR_v2']*stats['gradU_zpress_mean_y'] + stats['avR_vw']*stats['gradU_ypress_mean_y'])
            P_2D = P_xx + P_yy + P_zz + P_xy + P_xz + P_yz
            return P_2D

        print(Production_budget_2D(stats))
#%%
    # Dissipation term
        def Dissipation_budget_2D(stats):
            E_xx = -2*(stats['avegradU_budgets_x'] - stats['gradU_xpress_mean_x']**2 - stats['gradU_xpress_mean_y']**2 - stats['gradU_xpress_mean_z']**2)
            E_yy = -2*(stats['avegradU_budgets_y'] - stats['gradU_ypress_mean_x']**2 - stats['gradU_ypress_mean_y']**2 - stats['gradU_ypress_mean_z']**2)
            E_zz = -2*(stats['avegradU_budgets_y'] - stats['gradU_zpress_mean_x']**2 - stats['gradU_zpress_mean_y']**2 - stats['gradU_zpress_mean_z']**2)
            E_xy = -2*(stats['avegradU_budgets_xy'] - stats['gradU_xpress_mean_x']*stats['gradU_ypress_mean_x'] - stats['gradU_xpress_mean_y']*stats['gradU_ypress_mean_y'] - stats['gradU_ypress_mean_z']*stats['gradU_xpress_mean_z'])
            E_xz = -2*(stats['avegradU_budgets_xz'] - stats['gradU_xpress_mean_x']*stats['gradU_zpress_mean_x'] - stats['gradU_xpress_mean_y']*stats['gradU_zpress_mean_y'] - stats['gradU_xpress_mean_z']*stats['gradU_zpress_mean_z'])
            E_yz = -2*(stats['avegradU_budgets_yz'] - stats['gradU_ypress_mean_x']*stats['gradU_zpress_mean_x'] - stats['gradU_ypress_mean_y']*stats['gradU_zpress_mean_y'] - stats['gradU_ypress_mean_z']*stats['gradU_zpress_mean_z'])
            E_3D = E_xx + E_yy + E_zz + E_xy + E_xz + E_yz
            return E_3D
 #%%       
    #Viscous diffusion tensor
        def Viscous_budget_2D(stats):
            D_xx = viscosity *(stats[''])
            return V_2D

else:
    print("RSS calculations {}".format(rank))

comm.Barrier
"""
############# PLOTTING #################
#%%
fixed_x = stats['x'].iloc[69]
fixed_z = stats['z'].iloc[33]
print("x", fixed_x)
print("z", fixed_z)

matching_indices = np.isclose(stats['x'], fixed_x, atol=tolerance) & np.isclose(stats['z'] , fixed_z, atol=tolerance)

# Filter the data based on the matching indices
matching_x = [stats['x'][i] for i, match in enumerate(matching_indices) if match]
matching_z = [stats['z'][i] for i, match in enumerate(matching_indices) if match]
column_names = list(['avR_u2','avR_uv','avR_uw','avR_v2','avR_vw','avR_w2','avvel_x','avve2_x','avvel2'])
filtered_data = stats[matching_indices].sort_values(by=['y'])

for col_name in column_names:
    plt.scatter(filtered_data['y'], filtered_data[col_name], label='Scatter')

    # Create a line plot to join the dots smoothly
    plt.plot(filtered_data['y'], filtered_data[col_name], label='Line', color='red')

    plt.xlabel('y')
    plt.ylabel(col_name)
    plt.title('{} for x={}, z={}'.format(col_name, fixed_x, fixed_z))
    plt.show()

    #Save the plot as a .png file with the col_name as the file name
   # plt.savefig('/home/marcial/Simulacionprueba/Alvis/{}.png'.format(col_name))

    # Close the current plot to release resources
    plt.close()


# %%
