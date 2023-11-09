#%%

import numpy as np
import h5py as h5
import pandas as pd
import matplotlib as plt

from mpi4py import MPI
from pandas import HDFStore
import matplotlib.pyplot as plt

############# CASE SET UP ##################
Rtau = 180
Re =  np.exp((1.0/0.88)*np.log(Rtau/0.09))
mu = 2/Re
utau = (Rtau*mu)
Ny = 30
p = 4
##### VALIDATION LEE CHANNEL #####
data1 = np.loadtxt('/home/mech/sanchis/sanchis/channel_gpu/Validation/datachannel/LeeMoser_chan180/profiles/chan180.means', skiprows=1)
data = np.loadtxt('/home/mech/sanchis/sanchis/channel_gpu/Validation/datachannel/LeeMoser_chan180/profiles/chan180.reystress', skiprows=1)  # Skip the first row if it contains column names
y = data[:, 1]  # Extract the first column (y values)
r_uu = data[:, 2]
r_uv = data[:,5]
y1 = data1[:,1]
U = data1[:,2]
print(y)

# List of file paths to open iteratively

filenames = [
    "/home/mech/sanchis/sanchis/channel_gpu/channel_incom/180_6M/para_fine_2/results_equi_AVGRSS_channel-1_800001.h5",
    "/home/mech/sanchis/sanchis/channel_gpu/channel_incom/180_6M/para_fine_2/results_equi_AVGRSS_channel-1_900001.h5",
    "/home/mech/sanchis/sanchis/channel_gpu/channel_incom/180_6M/para_fine_2/results_equi_AVGRSS_channel-1_1000001.h5"
]

"""
filenames = [
    "/home/mech/sanchis/sanchis/channel_gpu/channel_incom/180_6M/para/results_equi_AVGRSS_channel-1_700001.h5",
    "/home/mech/sanchis/sanchis/channel_gpu/channel_incom/180_6M/para/results_equi_AVGRSS_channel-1_800001.h5",
    "/home/mech/sanchis/sanchis/channel_gpu/channel_incom/180_6M/para/results_equi_AVGRSS_channel-1_900001.h5"
]


filenames = [
    #"/home/mech/sanchis/sanchis/channel_gpu/Alvis/180_Refine/results_equi_AVGRSS_channel-1_100001.h5",
   #"/home/mech/sanchis/sanchis/channel_gpu/Alvis/180_Refine/results_AVGRSS_channel-1_200001.h5",
    "/home/mech/sanchis/sanchis/channel_gpu/Alvis/results_AVGRSS_channel-1_1.h5"
]
"""
N = 3

dataset = {}
periodic = 1
tol = 1e-6
tolerance = 1e-6

# Initialize an empty dataframe to store the combined data
stats = pd.DataFrame()
stats_final = pd.DataFrame()

print("Pre-postproc by rank 0, reading HDF5 files iteratively")

for filename in filenames:
    print(f"Opening file: {filename}")

    with h5.File(filename, "r") as f:
        # Print all root level object names (aka keys)
        # these can be group or dataset names
        print("Keys: %s" % list(f.keys()), len(f['y']))

        # Read the data from the HDF5 file and populate the dataset
        for col_name in list(f.keys()):
            #dataset[col_name] = f[col_name][:]
            stats[col_name] = f[col_name][:]
            
        stats_final = stats_final.add(stats, fill_value = 0)
        print('y 1', stats['y'].max())
        # Append the data from the current file to the combined_data dataframe
        #combined_data = combined_data.append(stats, ignore_index=True)



        
#%% 
    # Reynolds stresses

stats_final = stats_final.div(N)
stats_final.round(decimals = 8)
#stats_final *= Rtau
print(stats_final.y.nunique())
print('y', stats_final['y'].max(), len(stats_final['y']))

stats_final['avR_u2'] = stats_final['avve2_x'] - (stats_final['avvel_x'] * stats_final['avvel_x'])
stats_final['avR_uv'] = stats_final['avvex_x'] - stats_final['avvel_x'] * stats_final['avvel_y']
stats_final['avR_uw'] = stats_final['avvex_y'] - stats_final['avvel_x'] * stats_final['avvel_z']
stats_final['avR_v2'] = stats_final['avve2_y'] - stats_final['avvel_y'] * stats_final['avvel_y']
stats_final['avR_vw'] = stats_final['avvex_z'] - stats_final['avvel_y'] * stats_final['avvel_z']
stats_final['avR_w2'] = stats_final['avve2_z'] - stats_final['avvel_z'] * stats_final['avvel_z']
stats_final['avvel2'] = stats_final['avvel_x'] * stats_final['avvel_x']

############# AVERAGE SPACE ###################
averaged_data_final = pd.DataFrame()
averaged_data = stats_final.groupby(['y', 'x']).mean().reset_index()
averaged_data_1 = averaged_data.groupby(['y', 'z']).mean().reset_index()
averaged_data_final = averaged_data_1.groupby(['y']).mean().reset_index()

print(averaged_data_final['y'])

############ SYMETRIES ############
column_names = ['avegradU_budgets_x', 'avegradU_budgets_xy', 'avegradU_budgets_xz', 'avegradU_budgets_y', 'avegradU_budgets_yz', 'avegradU_budgets_z', 'avegradU_press_x1', 'avegradU_press_x2', 'avegradU_press_x3', 'avegradU_press_y1', 'avegradU_press_y2', 'avegradU_press_y3', 'avegradU_press_z1', 'avegradU_press_z2', 'avegradU_press_z3', 'avpre', 'avpre2', 'avpre4', 'avve2_x', 'avve2_y', 'avve2_z', 'avve3_x', 'avve3_y', 'avve3_z', 'avve4_x', 'avve4_y', 'avve4_z', 'avvel_x', 'avvel_y', 'avvel_z', 'avvelpr_x', 'avvelpr_y', 'avvelpr_z', 'avvex3_xy', 'avvex3_xyz', 'avvex3_xz', 'avvex3_yx', 'avvex3_yz', 'avvex3_zx', 'avvex3_zy', 'avvex_x', 'avvex_y', 'avvex_z', 'avvpre3']
averaged_data_final['y'] = averaged_data_final['y'].apply(lambda x: 2 - x if x >= 1 else x)
cols = ['avvel_x','avvel_y','avvel_z']
for col in cols:
    averaged_data_final[col] /= utau

cols = ['avR_u2', 'avR_v2', 'avR_w2', 'avve2_x', 'avvel2','avve2_y','avve2_z','avR_uw', 'avR_uv', 'avR_vw','avvex_x']
for col in cols:
    averaged_data_final[col] /= utau**2

cols = ['x','y','z']
for col in cols:
    averaged_data_final[col] *= Rtau

column_names_even = ['avR_u2', 'avR_v2', 'avR_w2', 'avvel_x', 'avve2_x', 'avvel2','avvel_y','avvel_z','avve2_y','avve2_z']
column_names_odd = ['avR_uw', 'avR_uv', 'avR_vw','avvex_x']
# Create a new DataFrame 'symmetric' with the same columns as 'averaged_data_final'
#symmetric = averaged_data_final.copy()
symmetric   = {}
symmetric_final = {}
target_keys = list(averaged_data_final.keys())
d = len(averaged_data_final['y'])
d2 = int(d/2+1)
d3 = int((Ny-1)*(p-1)+p)
for k in target_keys:
    symmetric[k] = np.empty((d2,))
    #symmetric_final[k] = np.empty((d2,))

symmetric = pd.DataFrame(symmetric)


print(d,d2)
# Iterate through each column and calculate symmetric averages
for col_name in column_names_even:
    for i in range(d2):
        symmetric.iloc[i, symmetric.columns.get_loc(col_name)] = (averaged_data_final[col_name][i] + averaged_data_final[col_name][d - i - 1]) / 2
    
for col_name in column_names_odd:
    for i in range(d2):
        symmetric.iloc[i, symmetric.columns.get_loc(col_name)] = (averaged_data_final[col_name][i] - averaged_data_final[col_name][d - i - 1]) / 2

for col_name in cols:
    for i in range(d2):
        symmetric.iloc[i, symmetric.columns.get_loc(col_name)] = (averaged_data_final[col_name][i])



    

print('s',symmetric['y'])
print(symmetric['y'].max())
print(len(symmetric['y']))

processed_indices = set()

# Iterate through the DataFrame
count = 0
cols = ['avegradU_budgets_x', 'avegradU_budgets_xy', 'avegradU_budgets_xz', 'avegradU_budgets_y', 'avegradU_budgets_yz', 'avegradU_budgets_z', 'avegradU_press_x1', 'avegradU_press_x2', 'avegradU_press_x3', 'avegradU_press_y1', 'avegradU_press_y2', 'avegradU_press_y3', 'avegradU_press_z1', 'avegradU_press_z2', 'avegradU_press_z3', 'avpre', 'avpre2', 'avpre4', 'avve2_x', 'avve2_y', 'avve2_z', 'avve3_x', 'avve3_y', 'avve3_z', 'avve4_x', 'avve4_y', 'avve4_z', 'avvel_x', 'avvel_y', 'avvel_z', 'avvelpr_x', 'avvelpr_y', 'avvelpr_z', 'avvex3_xy', 'avvex3_xyz', 'avvex3_xz', 'avvex3_yx', 'avvex3_yz', 'avvex3_zx', 'avvex3_zy', 'avvex_x', 'avvex_y', 'avvex_z', 'avvpre3']
#for col in cols:
processed_indices = set()
drop = set()
count = 1
k = 0
symmetric_final = symmetric.copy()
print('f',symmetric_final['y'])
for i in range(len(symmetric['y'])):
    if i not in processed_indices:
        for j in range(i + 1, len(symmetric['y'])):
            if (np.abs(symmetric['y'][i] - symmetric['y'][j]) <= tol and j not in drop):
                    count += 1
                    symmetric_final['avvel_x'][i] += symmetric['avvel_x'][j]
                    # Mark both indices as processed
                    processed_indices.add(i)
                    drop.add(j)

        symmetric_final['avvel_x'][i] /= count
        count = 1

symmetric_final = symmetric_final.drop(index=list(drop))
print('c',count,'ex',len(drop))
print(symmetric_final['avvel_x'])
print(symmetric['avvel_x'])
print('s2',len(symmetric_final['y']))
#print('0',)  




"""
#%% 
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


############# PLOTTING 3D #################
#%%
fixed_x = stats_final['x'].iloc[179994]
fixed_z = stats_final['z'].iloc[179994]
print("x", fixed_x)
print("z", fixed_z)

matching_indices = np.isclose(stats_final['x'], fixed_x, atol=tolerance) & np.isclose(stats_final['z'] , fixed_z, atol=tolerance)
filtered_data = stats_final[matching_indices].sort_values(by=['y'])

# Filter the data based on the matching indices
matching_x = [stats_final['x'][i] for i, match in enumerate(matching_indices) if match]
matching_z = [stats_final['z'][i] for i, match in enumerate(matching_indices) if match]
column_names = list(['avR_u2','avR_uv','avR_uw','avR_v2','avR_vw','avR_w2','avvel_x','avve2_x','avvel2'])

############# FILTERED DATA #########################
indices_with_zero_avvel_x = filtered_data.index[filtered_data['avvel_x'] == 0].tolist()
indices_with_zero_y = filtered_data.index[filtered_data['y'] == 0].tolist()
print('y', filtered_data['y'].max())

for col_name in column_names:
    plt.scatter(filtered_data['y'], filtered_data[col_name], label='Scatter')

    # Create a line plot to join the dots smoothly
    plt.plot(filtered_data['y'], filtered_data[col_name], label='Line', color='red')

    plt.xlabel('y')
    plt.ylabel(r"$\overline{w^2}$")
    #plt.ylabel(r"$V$")
    plt.title('{} for x={}, z={}'.format(col_name, fixed_x, fixed_z))
    plt.show()

    #Save the plot as a .png file with the col_name as the file name
    #plt.savefig('/home/mech/sanchis/sanchis/channel_gpu/Visualizations/{}.png'.format(col_name))

    # Close the current plot to release resources
    plt.close()
################## AVERAGED PLOTING ##########################
column_name = list(['avve2_x'])#,'avR_uv','avR_uw','avR_v2','avR_vw','avR_w2','avvel_x','avve2_x','avvex_y'])
for col_name in column_name:
    plt.scatter(averaged_data_final['y'], averaged_data_final[col_name], label='Scatter')

    # Create a line plot to join the dots smoothly
    plt.plot(averaged_data_final['y'], averaged_data_final[col_name], label='Line', color='red')

    plt.xlabel(r'$y^{+}$')
    plt.ylabel(r"$\overline{\tilde{vw}}$")
    #plt.ylabel(r"$W$")
    plt.show()

    #Save the plot as a .png file with the col_name as the file name
    #plt.savefig('/home/mech/sanchis/sanchis/channel_gpu/Visualizations/{}.png'.format(col_name))

    # Close the current plot to release resources
    plt.close()


################## AVERAGED PLOTING SYMETRIES ##########################
column_name = list(['avvel_x'])#,'avR_uv','avR_uw','avR_v2','avR_vw','avR_w2','avvel_x','avve2_x','avvel2'])

for col_name in column_name:
    plt.xscale("log")
    plt.scatter(symmetric.iloc[:d2,symmetric.columns.get_loc('y')], symmetric.iloc[:d2,symmetric.columns.get_loc(col_name)], label='Scatter')

    # Create a line plot to join the dots smoothly
    plt.plot(symmetric.iloc[:d2,symmetric.columns.get_loc('y')], symmetric.iloc[:d2,symmetric.columns.get_loc(col_name)], label='Line', color='red')

    plt.xlabel(r'$y^{+}$')
    plt.ylabel(r"$\overline{uw}$")
    #plt.ylabel(r"$W$")
    plt.show()

    #Save the plot as a .png file with the col_name as the file name
    #plt.savefig('/home/mech/sanchis/sanchis/channel_gpu/Visualizations/{}.png'.format(col_name))

    # Close the current plot to release resources
    plt.close()
"""
############ COMPARISON VS LEE ################
plt.figure(figsize=(8, 6))
plt.plot(y, r_uu, label='Lee&Moser', linestyle='-')  # You can adjust the linestyle as needed
# Scatter points on top of the line plot
plt.scatter(symmetric_final.iloc[:d2, symmetric.columns.get_loc('y')], symmetric_final['avR_u2'][:d2], label='Sod', color='red')
plt.xscale("log")
plt.xlabel('y+')
plt.ylabel(r"$\overline{u^2}^+$")  # Use the column name for the y-axis label
plt.title('Sod2d vs. Lee & Moser')
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(y, r_uv, label='Lee&Moser', linestyle='-')  # You can adjust the linestyle as needed
# Scatter points on top of the line plot
plt.scatter(symmetric_final.iloc[:d2, symmetric.columns.get_loc('y')], symmetric_final['avR_uv'][:d2], label='Sod', color='red')
plt.xscale("log")
plt.xlabel('y+')
plt.ylabel(r"$\overline{uv}^+$")  # Use the column name for the y-axis label
plt.title('Sod2d vs. Lee & Moser')
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(y1, U, label='Lee&Moser', linestyle='-')  # You can adjust the linestyle as needed
# Scatter points on top of the line plot
plt.scatter(symmetric_final.iloc[:d2, symmetric.columns.get_loc('y')], symmetric_final['avvel_x'][:d2], label='Sod', color='red')
plt.xscale("log")
plt.xlabel('y+')
plt.ylabel(r"$U^+$")  # Use the column name for the y-axis label
plt.title('Sod2d vs. Lee & Moser')
plt.legend()
plt.grid(True)
plt.show()

