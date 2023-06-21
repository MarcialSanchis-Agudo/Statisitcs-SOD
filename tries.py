import numpy as np
import h5py as h5
import pandas as pd
import matplotlib as plt

from mpi4py import MPI
from pandas import HDFStore



comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
filename ="results_AVGRSS_channel-3_1.h5"
filename1 = "Re180.078.h5"
dataset = {}
periodic = 1

if rank == 0:
    print("Pre-postproc by rank 0, reading hdf5 file {}".format(rank))
    with h5.File(filename, "r") as f:
        # Print all root level object names (aka keys) 
        # these can be group or dataset names 
        
        print("Keys: %s" % f.keys())
        column_names = list(f.keys())


    # Create a dictionary with column names as keys and empty lists as values
        stats = pd.DataFrame()
        dataset = {col_name: [] for col_name in column_names}
    
        
    # Read the data from the HDF5 file and populate the dataset
        for col_name in column_names:
            dataset[col_name] = f[col_name][:]
            stats[col_name] = f[col_name][:]
        
        print(stats['avvel_x'])
        counts = stats['x'].value_counts()
        print(counts)

        stats['y'] = stats['y'].astype(float)
        stats['x'] = stats['x'].astype(float)
        print(stats['avvel_x'])

    with h5.File(filename1, "r") as g:
        # Print all root level object names (aka keys) 
        # these can be group or dataset names 
        
        print("Keys: %s" % g.keys())
        column_names1 = list(g.keys())


    # Create a dictionary with column names as keys and empty lists as values
        stats1 = pd.DataFrame()
    
        
    # Read the data from the HDF5 file and populate the dataset
        for col_name1 in column_names1:
            stats1[col_name1] = g[col_name1][:]
        


### STRESSES and BUDGETS
# Print the result
# Calculate the mean for every column for rows with the same 'x' and 'y' values
        mean_by_xy = stats.groupby(['x', 'y']).mean()
        print(len(mean_by_xy))

# Print the result
        print(counts)
        
        if periodic == 1:
    # Homogeneous in z    
            for col_name in column_names:
                if col_name not in ['x', 'y', 'z']:
                    group_sums = stats.groupby(['x','y'])[col_name].transform('sum')
                    group_counts = stats.groupby(['x', 'y'])[col_name].transform('count')

    # Calculate the average for each group
                    group_averages = group_sums / group_counts

    # Replace the values in the DataFrame
                    stats[col_name] = group_averages
        
    # updated dataframe
            stats = stats.drop_duplicates(['x', 'y'])

        counts = stats['x'].value_counts()

# Print the result
        print(counts)
        print(stats['x'])


        
        


    # Reynolds stresses
        stats['avR_u2'] = stats['avve2_x'] - stats['avvel_x'] * stats['avvel_x']
        stats['avR_uv'] = stats['avvex_x'] - stats['avvel_x'] * stats['avvel_y']
        stats['avR_uw'] = stats['avvex_y'] - stats['avvel_x'] * stats['avvel_z']
        stats['avR_v2'] = stats['avve2_y'] - stats['avvel_y'] * stats['avvel_y']
        stats['avR_vw'] = stats['avvex_z'] - stats['avvel_y'] * stats['avvel_z']
        stats['avR_w2'] = stats['avve2_z'] - stats['avvel_z'] * stats['avvel_z']


    
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
        
    #Viscous diffusion tensor
        def Viscous_budget_2D(stats):
            return 1

else:
    print("RSS calculations {}".format(rank))

comm.Barrier