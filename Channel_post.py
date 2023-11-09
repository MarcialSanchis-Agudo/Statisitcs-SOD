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
rp = 8
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
"""
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

"""
filenames = [
    #"/home/mech/sanchis/sanchis/channel_gpu/Alvis/180_Refine/results_equi_AVGRSS_channel-1_100001.h5",
   #"/home/mech/sanchis/sanchis/channel_gpu/Alvis/180_Refine/results_AVGRSS_channel-1_200001.h5",
    "/home/mech/sanchis/sanchis/channel_gpu/Alvis/results_AVGRSS_channel-1_1.h5"
]

filenames = [
    "/home/mech/sanchis/sanchis/channel_gpu/channel_incom/180_6M/Bump/results_AVGRSS_channel-1_1.h5",
    #"/home/mech/sanchis/sanchis/channel_gpu/channel_incom/180_6M_2/results_AVGRSS_channel-1_400001.h5"
]
"""
N = 3

dataset = {}
periodic = 1
tol = 1e-8
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
        print("Keys: %s" % list(f.keys()), len(f['y']))

        # Read the data from the HDF5 file and populate the dataset
        for col_name in list(f.keys()):
            #dataset[col_name] = f[col_name][:]
            stats[col_name] = f[col_name][:]
            
        stats_final = stats_final.add(stats, fill_value = 0)
        print('y 1', stats['y'].max())
        # Append the data from the current file to the combined_data dataframe
        #combined_data = combined_data.append(stats, ignore_index=True)



        

######################################################## Reynolds stresses ##########################################

stats_final = stats_final.div(N)
print(2*rp-rp/4)
stats_final = stats_final.round(decimals = (int(2*rp-rp/4)))
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

##################################################### AVERAGE SPACE ####################################################
averaged_data_final = pd.DataFrame()
averaged_data = stats_final.groupby(['y', 'x']).mean().reset_index()
averaged_data_final = averaged_data.groupby(['y', 'z']).mean().reset_index()
averaged_data_final = averaged_data_final.groupby(['y']).mean().reset_index()
averaged_data_final = averaged_data_final.sort_values(by=['y'])

print('y unique', stats_final.y.nunique())
print('x unique', stats_final.y.nunique())
print('z unique', stats_final.y.nunique())
print(' final y', averaged_data_final.y.nunique())

"""
count = 0
cols = ['avR_uv','avR_uw','avR_v2','avR_vw','avR_w2','avvel_x','avve2_x','avvel2','avR_u2']
#for col in cols:
processed_indices = set()
drop = set()
count = 1
k = 0
#symmetric_final = symmetric.copy()
print('f',averaged_data_final['y'])

for col in cols:
    for i in range(len(averaged_data_final['y'])):
        if i not in processed_indices:
            for j in range(i + 1, len(averaged_data_final['y'])):
                if (np.abs(averaged_data_final['y'][i] - averaged_data_final['y'][j]) <= tol and j not in drop):
                    count += 1
                    averaged_data_final[col][i] += averaged_data_final[col][j]
                    # Mark both indices as processed
                    processed_indices.add(i)
                    drop.add(j)

        averaged_data_final[col][i] /= count
        count = 1

averaged_data_final = averaged_data_final.drop(index=list(drop))
averaged_data_final = averaged_data_final.reset_index(drop=True)
print('c',count,'ex',len(drop))
print(averaged_data_final['avvel_x'])
print('s2',len(averaged_data_final['y']))
print(averaged_data_final['y'])
"""

################################################ SYMETRIES ######################################################

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

print(averaged_data_final['y'])
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


############ ################################  COMPARISON VS LEE ########################################################
plt.figure(figsize=(8, 6))
plt.plot(y, r_uu, label='Lee&Moser', linestyle='-')  # You can adjust the linestyle as needed
# Scatter points on top of the line plot
plt.scatter(symmetric.iloc[:d2, symmetric.columns.get_loc('y')], symmetric['avR_u2'][:d2], label='Sod', color='red')
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
plt.scatter(symmetric.iloc[:d2, symmetric.columns.get_loc('y')], symmetric['avR_uv'][:d2], label='Sod', color='red')
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
plt.scatter(symmetric.iloc[:d2, symmetric.columns.get_loc('y')], symmetric['avvel_x'][:d2], label='Sod', color='red')
plt.xscale("log")
plt.xlabel('y+')
plt.ylabel(r"$U^+$")  # Use the column name for the y-axis label
plt.title('Sod2d vs. Lee & Moser')
plt.legend()
plt.grid(True)
plt.show()