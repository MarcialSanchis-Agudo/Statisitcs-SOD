import numpy as np
import h5py as h5
import pandas as pd
import matplotlib as plt
import math
from pathlib import Path

from mpi4py import MPI
from pandas import HDFStore
import matplotlib.pyplot as plt

############# CASE SET UP ##################
Rtau = 180
Re =  np.exp((1.0/0.88)*np.log(Rtau/0.09))
mu = 2/Re
utau = (Rtau*mu)
print('mu',mu,'expected utau',utau)
Ny = 16
p = 4
rp = 8
##### VALIDATION LEE CHANNEL #####
data1 = np.loadtxt('/home/mech/sanchis/sanchis/channel_gpu/Validation/datachannel/LeeMoser_chan180/profiles/chan180.means', skiprows=1)
data = np.loadtxt('/home/mech/sanchis/sanchis/channel_gpu/Validation/datachannel/LeeMoser_chan180/profiles/chan180.reystress', skiprows=1)
fig_path  = "Figs/"
Path(fig_path).mkdir(exist_ok=True)
try:
    # Load data using np.genfromtxt, skipping the first two rows
    data_jimenez = np.loadtxt('/home/mech/sanchis/sanchis/channel_gpu/Validation/datachannel/Jimenez_hoyas/Re950/profiles/Re950.prof')
    # Print the loaded data
except ValueError as e:
    print("Error occurred while loading the file:", e)
    
"""
y = data_jimenez[:, 0]  # Extract the first column (y values)
r_uu = data_jimenez[:, 3]
r_uv = data_jimenez[:,10]
r_vv = data_jimenez[:, 4]
r_ww = data_jimenez[:, 5]
r_uw = data_jimenez[:,11]
r_vw = data_jimenez[:,12]
y1 = data_jimenez[:,1]
U = data_jimenez[:,2]
W = data_jimenez[:,4]
"""
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
"""
tol = 1e-15
print('y',y[1],'y1',y1[1],'U',U[1])
# List of file paths to open iteratively


# Good dissipative
filenames = [
    "/home/mech/sanchis/sanchis/channel_gpu/Berzelius/data/tol/6_n/data_new_10/results_AVGRSS_final_channel-1_32678.h5",
]

filenames = [
    "/home/mech/sanchis/sanchis/channel_gpu/Berzelius/data/tol/6_n/data/results_AVGRSS_final_channel-1_32607.h5",
]

filenames = [
    "/home/mech/sanchis/sanchis/channel_gpu/Berzelius/data/tol/6_n/les_16/results_AVGRSS_final_channel-1_32541.h5",
]

filenames = [
    "/home/mech/sanchis/sanchis/channel_gpu/Berzelius/data/tol/6_n_2/dat_les/results_AVGRSS_final_channel-1_32762.h5",
]

########### VISHAL ###########
filenames = [
    "/home/mech/sanchis/sanchis/channel_gpu/Berzelius/data/tol/6_n_3/data_large/results_AVGRSS_final_channel-1_32564.h5"
]

"""
############ MARCIAL ###################
filenames = [
    "/home/mech/sanchis/sanchis/channel_gpu/Berzelius/DNS_180_30/data/results_AVGRSS_final_channel-1_5243.h5"
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

"""
stats_final['avR_u2'] = stats_final['avve2_x'] - (stats_final['avvel_x'] * stats_final['avvel_x'])
stats_final['avR_uv'] = stats_final['avvex_x'] - stats_final['avvel_x'] * stats_final['avvel_y']
stats_final['avR_uw'] = stats_final['avvex_y'] - stats_final['avvel_x'] * stats_final['avvel_z']
stats_final['avR_v2'] = stats_final['avve2_y'] - stats_final['avvel_y'] * stats_final['avvel_y']
stats_final['avR_vw'] = stats_final['avvex_z'] - stats_final['avvel_y'] * stats_final['avvel_z']
stats_final['avR_w2'] = stats_final['avve2_z'] - stats_final['avvel_z'] * stats_final['avvel_z']
stats_final['avvel2'] = stats_final['avvel_x'] * stats_final['avvel_x']
"""

print('check', stats_final['avvel_x'][1],stats_final['avvel_x'][50])
pd.set_option("display.max_rows", None)
############### SORTING AND AVERAGIN NODES #########################
stats_final = stats_final.sort_values("y", ascending=True)
stats_final = stats_final.reset_index(drop=True)
stats_final = stats_final.groupby(['y']).mean().reset_index()
##################################################################
print("before",stats_final['y'])
print(stats_final.y.nunique())
print('y', stats_final['y'].max(), len(stats_final['y']),'avvel_x',stats_final['avvel_x'].max())
stats_final['utau'] = np.sqrt(mu*stats_final['gradU_xpress_mean_y'])
utau = stats_final['utau'][0]
Rtau = utau/mu
print('utau', utau, 'Rtau', Rtau)

column_name = list(['avvel_x','avve2_x','avR_u2','avR_uv','Production','Dissipation'])#,'avR_v2','avR_vw','avR_w2','avvel_x','avve2_x','avvex_y'])
for col_name in column_name:
    # Create a line plot to join the dots smoothly
    plt.plot(stats_final['y'][:], stats_final[col_name][:], label='Line', color='red')

    plt.xlabel(r'$y^{+}$')
    plt.ylabel(r"$\overline{u^2}$")
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

cols = ['avR_u2', 'avR_v2', 'avR_w2', 'avve2_x','avve2_y','avve2_z','avR_uw', 'avR_uv', 'avR_vw','avvex_x']
for col in cols:
    stats_final[col] /= utau**2

cols = ['gradU_xpress_mean_y']
for col in cols:
    stats_final[col] *= (utau**2/mu)
cols = ['x','y','z']
for col in cols:
    stats_final[col] *= Rtau

column_names_even = ['avR_u2', 'avR_v2', 'avR_w2', 'avvel_x', 'avve2_x','avvel_y','avvel_z','avve2_y','avve2_z']
column_names_odd = ['avR_uw', 'avR_uv', 'avR_vw','avvex_x','gradU_xpress_mean_y']

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


################ convergence ########################
symmetric['convergence'] = symmetric['gradU_xpress_mean_y'] + symmetric['avR_uv'] - 1 + symmetric['y']
plt.plot(symmetric['y'][:d2], symmetric['gradU_xpress_mean_y'][:d2], label='molecular')
plt.plot(symmetric['y'][:d2],- symmetric['avR_uv'][:d2], label='turbulent')
plt.plot(symmetric['y'][:d2], 1 - stats_final['y'][:d2]/Rtau, label='total')
plt.xscale('log')
# Adding labels and a legend
plt.xlabel(r'$y^{+}$')
plt.ylabel('Convergence terms')
plt.title('Multivariable plot')
plt.legend()

# Show the plot
plt.show()

plt.plot(symmetric['y'][:d2],symmetric['gradU_xpress_mean_y'][:d2]/Rtau - symmetric['avR_uv'][:d2] - 1 + symmetric['y'][:d2]/Rtau, label='convergence')
plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r'$y^{+}$')
plt.ylabel(r"$\frac{dU^{+}}{dy^{+}} - \overline{uv}^+ - 1 + y^+$")
plt.title('Convergence')
plt.legend()
plt.show()

######################################################




############## NO SYMM #########################
column_name = list(['gradU_xpress_mean_y'])#,'avR_uv','avR_uw','avR_v2','avR_vw','avR_w2','avvel_x','avve2_x','avvex_y'])
for col_name in column_name:
    plt.scatter(symmetric.iloc[:d2, symmetric.columns.get_loc('y')], symmetric['convergence'][:d2], label='Sod', color='blue')

    # Create a line plot to join the dots smoothly
    #plt.plot(stats_final['y'][:d2], stats_final[col_name][:d2], label='Line', color='red')

    plt.xlabel(r'$y^{+}$')
    #plt.ylabel(r"$\overline{\tilde{vw}}$")
    plt.ylabel(r"Convergence")
    plt.show()

    #Save the plot as a .png file with the col_name as the file name
    #plt.savefig('/home/mech/sanchis/sanchis/channel_gpu/Visualizations/{}.png'.format(col_name))

    # Close the current plot to release resources
    plt.close()

############################################################################


############ ################################  COMPARISON VS LEE VS Jimenez ########################################################
plt.figure(figsize=(8, 6))
plt.plot(y1, U, label='Lee&Moser', linestyle='-')  # You can adjust the linestyle as needed
# Scatter points on top of the line plot
plt.scatter(symmetric.iloc[:, symmetric.columns.get_loc('y')], symmetric['avvel_x'][:], label='Sod', color='red')
plt.xscale("log")
#plt.yscale('log')
plt.xlabel('y+')
plt.ylabel(r"$U^+$")  # Use the column name for the y-axis label
plt.title('Sod2d vs. Lee & Moser')
plt.legend()
plt.grid(True)
plt.savefig(fig_path + f'U.jpg',bbox_inches='tight',dpi=1000)
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
plt.savefig(fig_path + f'W.jpg',bbox_inches='tight',dpi=1000)
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
plt.savefig(fig_path + f'u2.jpg',bbox_inches='tight',dpi=1000)
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
plt.savefig(fig_path + f'v2.jpg',bbox_inches='tight',dpi=1000)
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
plt.savefig(fig_path + f'w2.jpg',bbox_inches='tight',dpi=1000)
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
plt.savefig(fig_path + f'uv.jpg',bbox_inches='tight',dpi=1000)
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(y, r_uw, label='Lee&Moser', linestyle='-')  # You can adjust the linestyle as needed
# Scatter points on top of the line plot
plt.scatter(symmetric.iloc[:d2, symmetric.columns.get_loc('y')], symmetric['avR_uw'][:d2], label='Sod', color='red')
plt.xscale("log")
plt.xlabel('y+')
plt.ylabel(r"$\overline{uw}^+$")  # Use the column name for the y-axis label
plt.title('Sod2d vs. Lee & Moser')
plt.legend()
plt.grid(True)
plt.savefig(fig_path + f'uw.jpg',bbox_inches='tight',dpi=1000)
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
plt.savefig(fig_path + f'vw.jpg',bbox_inches='tight',dpi=1000)
plt.show()
