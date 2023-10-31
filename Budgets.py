#%%

import numpy as np
import h5py as h5
import pandas as pd
import matplotlib as plt

from mpi4py import MPI
from pandas import HDFStore
import matplotlib.pyplot as plt

"""

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
#filename ="/home/marcial/Simulacionprueba/rank_4/channel_val_elapsed/results_AVGRSS_channel-4_21.h5"
filename ="/home/marcial/Simulacionprueba/Alvis/results_AVGRSS_channel-4_200001.h5"
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

"""

# List of file paths to open iteratively
filenames = [
    "/home/mech/sanchis/sanchis/channel_gpu/channel_incom/180_6M/para/results_AVGRSS_channel-1_1.h5",
    "/home/mech/sanchis/sanchis/channel_gpu/channel_incom/180_6M/para/results_AVGRSS_channel-1_100001.h5"
]
N = 2

dataset = {}
periodic = 1
tolerance = 1e-8

# Initialize an empty dataframe to store the combined data
stats = pd.DataFrame()
stats_final = pd.DataFrame()

print("Pre-postproc by rank 0, reading HDF5 files iteratively")

for filename in filenames:
    print(f"Opening file: {filename}")

    with h5.File(filename, "r") as f:
        # Print all root level object names (aka keys)
        # these can be group or dataset names
        print("Keys: %s" % list(f.keys()))


        # Read the data from the HDF5 file and populate the dataset
        for col_name in list(f.keys()):
            #dataset[col_name] = f[col_name][:]
            stats[col_name] = f[col_name][:]
            
        stats_final = stats_final.add(stats, fill_value = 0)
        # Append the data from the current file to the combined_data dataframe
        #combined_data = combined_data.append(stats, ignore_index=True)



        
#%% 
    # Reynolds stresses

stats_final = stats_final.div(3)
print(stats_final['avvel_x'])
stats_final['avR_u2'] = stats_final['avve2_x'] - (stats_final['avvel_x'] * stats_final['avvel_x'])
stats_final['avR_uv'] = stats_final['avvex_x'] - stats_final['avvel_x'] * stats_final['avvel_y']
stats_final['avR_uw'] = stats_final['avvex_y'] - stats_final['avvel_x'] * stats_final['avvel_z']
stats_final['avR_v2'] = stats_final['avve2_y'] - stats_final['avvel_y'] * stats_final['avvel_y']
stats_final['avR_vw'] = stats_final['avvex_z'] - stats_final['avvel_y'] * stats_final['avvel_z']
stats_final['avR_w2'] = stats_final['avve2_z'] - stats_final['avvel_z'] * stats_final['avvel_z']
stats_final['avvel2'] = stats_final['avvel_x'] * stats_final['avvel_x']


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
#fixed_x = 0.3867764603626959 
#fixed_z = 1.506120399122544 
fixed_x = stats_final['x'].iloc[179994]
fixed_z = stats_final['z'].iloc[179994]
print("x", fixed_x)
print("z", fixed_z)

matching_indices = np.isclose(stats_final['x'], fixed_x, atol=tolerance) & np.isclose(stats_final['z'] , fixed_z, atol=tolerance)
#matching_indices = []
#for i in range(len(stats['x'])):
#    if fixed_x - tol < stats['x'].iloc[i] and stats['x'].iloc[i] < fixed_x + tol:
#        if fixed_z - tol < stats['z'].iloc[i] and stats['z'].iloc[i] < fixed_z + tol:
#            matching_indices.append(i)
#matching_indices_x = np.isclose(stats['x'], fixed_x, atol=tolerance)
#print('final', len(matching_indices), matching_indices)
filtered_data = stats_final[matching_indices].sort_values(by=['y'])

# Filter the data based on the matching indices
matching_x = [stats_final['x'][i] for i, match in enumerate(matching_indices) if match]
matching_z = [stats_final['z'][i] for i, match in enumerate(matching_indices) if match]
column_names = list(['avvel_x'])#,'avR_uv','avR_uw','avR_v2','avR_vw','avR_w2','avvel_x','avve2_x','avvel2'])

############# FILTERED DATA #########################
indices_with_zero_avvel_x = filtered_data.index[filtered_data['avvel_x'] == 0].tolist()
indices_with_zero_y = filtered_data.index[filtered_data['y'] == 0].tolist()
print(len(filtered_data['x']))


for col_name in column_names:
    plt.scatter(filtered_data['y'], filtered_data[col_name], label='Scatter')

    # Create a line plot to join the dots smoothly
    plt.plot(filtered_data['y'], filtered_data[col_name], label='Line', color='red')

    plt.xlabel('y')
    plt.ylabel(r"$\overline{ \tilde{w^2}} $")
    plt.title('{} for x={}, z={}'.format(col_name, fixed_x, fixed_z))
    plt.show()

    #Save the plot as a .png file with the col_name as the file name
   # plt.savefig('/home/marcial/Simulacionprueba/Alvis/{}.png'.format(col_name)),'avR_uv','avR_uw','avR_v2','avR_vw','avR_w2','avvel_x','avve2_x','avvel2'])

    # Close the current plot to release resources
    plt.close()


# %%
