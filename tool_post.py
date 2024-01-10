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
Ny = 19
p = 4
rp = 8
##### VALIDATION LEE CHANNEL #####
data1 = np.loadtxt('/home/mech/sanchis/sanchis/channel_gpu/Validation/datachannel/LeeMoser_chan180/profiles/chan180.means', skiprows=1)
data = np.loadtxt('/home/mech/sanchis/sanchis/channel_gpu/Validation/datachannel/LeeMoser_chan180/profiles/chan180.reystress', skiprows=1)  # Skip the first row if it contains column names
y = data[:, 1]  # Extract the first column (y values)
r_uu = data[:, 2]
r_uv = data[:,5]
r_vv = data[:, 3]
r_ww = data[:, 4]
r_uw = data[:,6]
r_vw = data[:,7]
y1 = data1[:,1]
U = data1[:,2]
W = data1[:,4]
tol = 1e-15
print(y)
# List of file paths to open iteratively

"""
filenames = [
    "/home/mech/sanchis/sanchis/channel_gpu/channel_incom/180_6M/para_fine/results_AVGRSS_final_channel-1_0.h5",
]
"""
filenames = [
    "/home/mech/sanchis/sanchis/channel_gpu/Berzelius/data/inter/results_AVGRSS_final_channel-1_0.h5",
]
print("Pre-postproc by rank 0, reading HDF5 files iteratively")

stats = pd.DataFrame()
stats_final = pd.DataFrame()
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
        print('x 1', stats['x'].max())
        print('z 1', stats['z'].max())
        # Append the data from the current file to the combined_data dataframe
        #combined_data = combined_data.append(stats, ignore_index=True)

print(stats_final.y.nunique())
print('y', stats_final['y'].max(), len(stats_final['y']),'avvel_x',stats_final['avvel_x'].max())

stats_final['avR_u2'] = stats_final['avve2_x'] - (stats_final['avvel_x'] * stats_final['avvel_x'])
stats_final['avR_uv'] = stats_final['avvex_x'] - stats_final['avvel_x'] * stats_final['avvel_y']
stats_final['avR_uw'] = stats_final['avvex_y'] - stats_final['avvel_x'] * stats_final['avvel_z']
stats_final['avR_v2'] = stats_final['avve2_y'] - stats_final['avvel_y'] * stats_final['avvel_y']
stats_final['avR_vw'] = stats_final['avvex_z'] - stats_final['avvel_y'] * stats_final['avvel_z']
stats_final['avR_w2'] = stats_final['avve2_z'] - stats_final['avvel_z'] * stats_final['avvel_z']
stats_final['avvel2'] = stats_final['avvel_x'] * stats_final['avvel_x']

print('check', stats_final['avvel_x'][1],stats_final['avvel_x'][145])
pd.set_option("display.max_rows", None)
stats_final = stats_final.sort_values("y", ascending=True)
stats_final = stats_final.reset_index(drop=True)
print("before",stats_final['y'])

column_name = list(['avR_u2','avR_uv','avR_uw','avR_v2','avR_vw','avR_w2','avvel_x','avve2_x','avvex_y'])
for col_name in column_name:
    # Create a line plot to join the dots smoothly
    plt.plot(stats_final['y'][:], stats_final[col_name][:], label='Line', color='red')

    plt.xlabel(r'$y^{+}$')
    plt.ylabel(r"$\overline{\tilde{vw}}$")
    #plt.ylabel(r"$W$")
    plt.show()

    #Save the plot as a .png file with the col_name as the file name
    #plt.savefig('/home/mech/sanchis/sanchis/channel_gpu/Visualizations/{}.png'.format(col_name))

    # Close the current plot to release resources
    plt.close()

############### SYMMETRIES ###################
stats_final['y'] = stats_final['y'].apply(lambda x: 2 - x if x >= 1 else x)
cols = ['avvel_x','avvel_y','avvel_z']
for col in cols:
    stats_final[col] /= utau

cols = ['avR_u2', 'avR_v2', 'avR_w2', 'avve2_x', 'avvel2','avve2_y','avve2_z','avR_uw', 'avR_uv', 'avR_vw','avvex_x']
for col in cols:
    stats_final[col] /= utau**2

cols = ['x','y','z']
for col in cols:
    stats_final[col] *= Rtau

column_names_even = ['avR_u2', 'avR_v2', 'avR_w2', 'avvel_x', 'avve2_x', 'avvel2','avvel_y','avvel_z','avve2_y','avve2_z']
column_names_odd = ['avR_uw', 'avR_uv', 'avR_vw','avvex_x']
# Create a new DataFrame 'symmetric' with the same columns as 'stats_final'
#symmetric = stats_final.copy()
symmetric   = {}
symmetric_final = {}
target_keys = list(stats_final.keys())
d = len(stats_final['y'])
d2 = int(d/2+1)
d3 = int((Ny-1)*(p-1)+p)
for k in target_keys:
    symmetric[k] = np.empty((d2,))
    #symmetric_final[k] = np.empty((d2,))

symmetric = pd.DataFrame(symmetric)

print(stats_final['y'])
print(d,d2)
# Iterate through each column and calculate symmetric averages
for col_name in column_names_even:
    for i in range(d2):
        symmetric.iloc[i, symmetric.columns.get_loc(col_name)] = (stats_final[col_name][i] + stats_final[col_name][d - i - 1]) / 2
    
for col_name in column_names_odd:
    for i in range(d2):
        symmetric.iloc[i, symmetric.columns.get_loc(col_name)] = (stats_final[col_name][i] - stats_final[col_name][d - i - 1]) / 2

print(cols)
for col_name in cols:
    for i in range(d2):
        symmetric.iloc[i, symmetric.columns.get_loc(col_name)] = (stats_final[col_name][i])



    

print('s',symmetric['y'])
maxValueIndex = symmetric['y'].idxmax()
print(symmetric['y'].max(),maxValueIndex)
print(len(symmetric['y']))

############## NO SYMM #########################
column_name = list(['avvel_x'])#,'avR_uv','avR_uw','avR_v2','avR_vw','avR_w2','avvel_x','avve2_x','avvex_y'])
for col_name in column_name:
    plt.scatter(symmetric.iloc[:d2, symmetric.columns.get_loc('y')], symmetric['avvel_x'][:d2], label='Sod', color='blue')

    # Create a line plot to join the dots smoothly
    plt.plot(stats_final['y'][:d2], stats_final[col_name][:d2], label='Line', color='red')

    plt.xlabel(r'$y^{+}$')
    plt.ylabel(r"$\overline{\tilde{vw}}$")
    #plt.ylabel(r"$W$")
    plt.show()

    #Save the plot as a .png file with the col_name as the file name
    #plt.savefig('/home/mech/sanchis/sanchis/channel_gpu/Visualizations/{}.png'.format(col_name))

    # Close the current plot to release resources
    plt.close()




############ ################################  COMPARISON VS LEE ########################################################
plt.figure(figsize=(8, 6))
plt.plot(y1, U, label='Lee&Moser', linestyle='-')  # You can adjust the linestyle as needed
# Scatter points on top of the line plot
plt.scatter(symmetric.iloc[:, symmetric.columns.get_loc('y')], symmetric['avvel_x'][:], label='Sod', color='red')
plt.xscale("log")
plt.xlabel('y+')
plt.ylabel(r"$U^+$")  # Use the column name for the y-axis label
plt.title('Sod2d vs. Lee & Moser')
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(y1, W, label='Lee&Moser', linestyle='-')  # You can adjust the linestyle as needed
# Scatter points on top of the line plot
plt.scatter(symmetric.iloc[:, symmetric.columns.get_loc('y')], symmetric['avvel_z'][:], label='Sod', color='red')
plt.xscale("log")
plt.xlabel('y+')
plt.ylabel(r"$W^+$")  # Use the column name for the y-axis label
plt.title('Sod2d vs. Lee & Moser')
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(y, r_uu, label='Lee&Moser', linestyle='-')  # You can adjust the linestyle as needed
# Scatter points on top of the line plot
plt.scatter(symmetric.iloc[:, symmetric.columns.get_loc('y')], symmetric['avR_u2'][:], label='Sod', color='red')
plt.xscale("log")
plt.xlabel('y+')
plt.ylabel(r"$\overline{u^2}^+$")  # Use the column name for the y-axis label
plt.title('Sod2d vs. Lee & Moser')
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(y, r_vv, label='Lee&Moser', linestyle='-')  # You can adjust the linestyle as needed
# Scatter points on top of the line plot
plt.scatter(symmetric.iloc[:d2, symmetric.columns.get_loc('y')], symmetric['avR_v2'][:d2], label='Sod', color='red')
plt.xscale("log")
plt.xlabel('y+')
plt.ylabel(r"$\overline{v^2}^+$")  # Use the column name for the y-axis label
plt.title('Sod2d vs. Lee & Moser')
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(y, r_ww, label='Lee&Moser', linestyle='-')  # You can adjust the linestyle as needed
# Scatter points on top of the line plot
plt.scatter(symmetric.iloc[:d2, symmetric.columns.get_loc('y')], symmetric['avR_w2'][:d2], label='Sod', color='red')
plt.xscale("log")
plt.xlabel('y+')
plt.ylabel(r"$\overline{w^2}^+$")  # Use the column name for the y-axis label
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
plt.plot(y, r_uw, label='Lee&Moser', linestyle='-')  # You can adjust the linestyle as needed
# Scatter points on top of the line plot
plt.scatter(symmetric.iloc[:d2, symmetric.columns.get_loc('y')], symmetric['avR_uw'][:d2], label='Sod', color='red')
plt.xscale("log")
plt.xlabel('y+')
plt.ylabel(r"$\overline{uv}^+$")  # Use the column name for the y-axis label
plt.title('Sod2d vs. Lee & Moser')
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(y, r_vw, label='Lee&Moser', linestyle='-')  # You can adjust the linestyle as needed
# Scatter points on top of the line plot
plt.scatter(symmetric.iloc[:d2, symmetric.columns.get_loc('y')], symmetric['avR_vw'][:d2], label='Sod', color='red')
plt.xscale("log")
plt.xlabel('y+')
plt.ylabel(r"$\overline{vw}^+$")  # Use the column name for the y-axis label
plt.title('Sod2d vs. Lee & Moser')
plt.legend()
plt.grid(True)
plt.show()