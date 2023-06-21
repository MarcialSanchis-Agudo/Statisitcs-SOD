#%%
import numpy as np
import h5py as h5
import pandas as pd
import matplotlib.pyplot as plt
import pymech.exadata
import pymech.neksuite as ns


from mpi4py import MPI
from pandas import HDFStore
from scipy import interpolate

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
existing_filename = "results_AVGRSS_channel-1_141.h5"
new_filename = "b08channelnek0.f00001"
existing_dataset = {}
new_dataset = {}
periodic = 0

#%%
if rank == 0:
    print("Pre-postproc by rank 0, reading existing HDF5 file: {}".format(existing_filename))
    with h5.File(existing_filename, "r") as existing_file:
        # Print all root level object names (aka keys)
        # these can be group or dataset names
        print("Keys: %s" % existing_file.keys())
        column_names = list(existing_file.keys())

        # Create a dictionary with column names as keys and empty lists as values
        existing_stats = pd.DataFrame()
        existing_dataset = {col_name: [] for col_name in column_names}

        # Read the data from the existing HDF5 file and populate the dataset
        for col_name in column_names:
            existing_dataset[col_name] = existing_file[col_name][:]
            existing_stats[col_name] = existing_file[col_name][:]

        counts = existing_stats['avegradU_press_x1'].value_counts()

        # Print the result
        print(counts)
        print(existing_stats['avve3_x'])

        plt.plot(existing_stats['avve3_x'])
        plt.xlabel('Index')
        plt.ylabel('Value')
        plt.title('avve3_x')
        plt.show()
#%%
print("\nReading new file: {}".format(new_filename))
#%%
with open(new_filename, "r") as new_file:
        # Process the new file data and compare with the existing data
        # Add your code to process the new file here 
        new_data = ns.readnek(new_filename)  # Load the data from the new file
        print("Read nek done")

    # Process the new data and compare with the existing data
        new_stats = pd.DataFrame(data=new_data, columns=column_names) 
        print("Keys nek: %s" % new_stats.keys())
# %%
