
#%%
import numpy as np
import h5py as h5
import numpy as np
import pandas as pd
import math
import pyAlya




#%%
#%%

#### SET UP DATA
CASESTR = 'channel2dpy.h5'
stats = "/home/marcial/Simulacionprueba/Alvis/results_AVGRSS_final_channel-1_14001.h5"
equi_mesh = "/home/marcial/TGV/mesh/big_603060/channel_equi-1.hdf"
mesh  = pyAlya.Mesh.load(CASESTR,compute_massMatrix=False)
stats_file = h5.File(stats, "r")
print(stats_file.keys())


#### ORDER STATS
def find_row(matrix, value1, value2, tolerance=1e-6):
    for i, row in enumerate(matrix):
        if math.isclose(row[0], value1, rel_tol=tolerance) and math.isclose(row[1], value2, rel_tol=tolerance):
            return i
    return -1 # Return -1 if the values are not found in any row


with h5.File(equi_mesh, "r") as f:
        # Print all root level object names (aka keys) 
        # these can be group or dataset names 

    column_names =  [name for name in stats_file.keys()]
    stats2D = pd.DataFrame()
    for col_name in column_names:
        stats2D[col_name] = stats_file[col_name][:]


    # Create a dictionary with column names as keys and empty lists as values

    #### EQUIDISTANT XYZ
    xyz = np.zeros((len(stats_file["x"]),2),np.double)
    nodes = np.zeros(len(stats_file["x"]))
    print(range(len(stats_file["x"])),stats_file["z"][0] )

    for i in range(len(stats_file["x"])):
        xyz[i,0] = np.array(f["VTKHDF"]["Points"][int(stats2D["nodes"][i])-1,0])
        xyz[i,1] = np.array(f["VTKHDF"]["Points"][int(stats2D["nodes"][i])-1,1])
        nodes[i] = find_row(mesh.xyz,xyz[i,0],xyz[i,1])
    
    #print(nodes)
    #stats_file.close()

    #with h5.File(stats,"a") as stats_file:
        
    #    stats_file.create_dataset(name= "nodes_order", data= np.array(nodes), dtype=int )
    
    #stats_file.close()
    #with h5.File(stats,"r") as stats_file:
        
    #    stats = stats_file.sort_values(by=['nodes_order'])
    #    print(stats_file.keys())
    
    #stats_file.close()
    stats2D["nodes_order"] = nodes
    stats = stats2D.sort_values(by=['nodes_order'])
    #print(stats["x"][:])



column_names = list(stats.keys())
print("Keys: ", column_names)

Var_dict = {}
for var in column_names:
    #if var in ["x","y","z"]: continue  

    Var_dict[var] = np.array(stats[var])

print("Var :", Var_dict["x"])


with h5.File(CASESTR, "r+") as f:
        # Print all root level object names (aka keys) 
        # these can be group or dataset names 
      
    print(f.keys())


    # Create a dictionary with column names as keys and empty lists as values
    f["MESH"]["xyz"][:,0] = Var_dict["x"][:]
    f["MESH"]["xyz"][:,1] = Var_dict["y"][:]

#print('xy mesh:', mesh.xyz)
#print('xy stats:', Var_dict["x"],Var_dict["y"] )



stats2D_field = pyAlya.Field(xyz=mesh.xyz, serial = True, **Var_dict)
stats2D_field.save(CASESTR)
pyAlya.pprint(0,stats2D_field)

#%%
CASESTR = 'channel3dpy.h5'
stats = "/home/marcial/Simulacionprueba/Alvis/results_AVGRSS_final_channel-1_100001.h5"
equi_mesh = "/home/marcial/TGV/mesh/big_603060/channel_equi-1.hdf"
mesh  = pyAlya.Mesh.load(CASESTR,compute_massMatrix=False)
stats_file = h5.File(stats, "r")

#%%
with h5.File(CASESTR, "r+") as f:
        # Print all root level object names (aka keys) 
        # these can be group or dataset names 
      
    print(f.keys())


    # Create a dictionary with column names as keys and empty lists as values
    f["MESH"]["xyz"][:,0] = Var_dict["x"][:]
    f["MESH"]["xyz"][:,1] = Var_dict["y"][:]
    
#%% #####################STATS 3D ############################
x0 = stats2D['x'].iloc[0]
z0 = stats2D['z'].iloc[1]
x1 = stats2D['x'].iloc[0]
z1 = stats2D['z'].iloc[1]
y0 = stats2D['y'].min(), 
y1 = stats2D['y'].max()

p1, p2 = pyAlya.Geom.Point(x0,y0,z0), pyAlya.Geom.Point(x1,y1,z1)
line   = pyAlya.Geom.Line(p1,p2,npoints=100)
fieldp = mesh.interpolate(line.xyz,**Var_dict,global_max_iter=3)
fieldp.clearNaN()
if pyAlya.utils.is_rank_or_serial(0): fieldp.save('line.h5',mpio=False)




