import numpy as np
#import matplotlib.pyplot as plt
import math
from scipy.special import legendre


# Specify the filename for the Gmsh .geo file
geo_filename = "channel_p5.geo"
Re=180

Lx = 2*math.pi; Lz = math.pi; Ly = 2.0; d0=1.0
Nx = 40; Nz = 40; Ny = 34; porder = 5; refine=3  #0: no-stretching, 1:tanh, 2:cos, 3:power-law  
# Important: Ny has to be even !!

# Wall-normal coordinates
if refine==0:
    y=np.linspace(0,Ly,ny)
    p = 17; per = Ny - p + 1; ny = p + per
elif refine==1:
    beta = 3.8 # Stretching factor
    y = np.tanh(beta * (np.linspace(0, 1, ny) - 0.5)) / np.tanh(beta * 0.5) + 1.0
    p = 17; per = Ny - p + 1; ny = p + per
elif refine==2:
    theta = np.linspace(0, np.pi, p + per)
    y = -np.cos(theta)+1
    p = 17; per = Ny - p + 1; ny = p + per
elif refine==3:
    def distribute_points_power_law(n_points, exponent):
        points = np.linspace(0, 1, n_points)
        scaled_points = points ** exponent
        scaled_points /= np.max(scaled_points)  # Normalize to ensure maximum value is 1
        return scaled_points
    # Parameters
    ny = Ny//2+1
    exponent = 1.2
    # Generate grid points
    y_1 = distribute_points_power_law(ny, exponent)
    y_2 = 2-np.flip(y_1)
    y = np.unique(np.concatenate((y_1, y_2)))
    p = 17; per = len(y) - p

print(len(y))
print(f"Number of elements in y Ny={len(np.diff(y))}:")

delta_z= (Lz/Nz)*Re
delta_x= (Lx/Nx)*Re
delta_wall = np.diff(y)[0]*Re
delta_mid = max(np.diff(y))*Re

print('-----------')
print('delta_wall = ',delta_wall)
print('delta_mid  = ',delta_mid)

print('-----------')
print('aspect_ratio_xy = ',delta_mid/delta_x)
print('aspect_ratio_xz = ',delta_x/delta_z)
print('aspect_ratio_yz = ',delta_mid/delta_z)
print('aspect_ratio_xy_wall = ',delta_x/delta_wall)
print('aspect_ratio_yz_wall = ',delta_z/delta_wall)

def gauss_lobatto_legendre_points(porder):
    degree=porder-1
    # Calculate the roots of the derivative of the Legendre polynomial of order n
    roots, _ = np.polynomial.legendre.leggauss(degree) # The function expects the degree i.e. order-1
    # Add the endpoints -1 and 1
    points = np.concatenate(([-1.0], roots, [1.0]))
    return points

# Calculate Gauss-Lobatto-Legendre points
points_x = (1 - gauss_lobatto_legendre_points(porder))*delta_x/2
points_z = (1 - gauss_lobatto_legendre_points(porder))*delta_z/2
points_wall = (1 - gauss_lobatto_legendre_points(porder))*delta_wall/2
points_mid  = (1 - gauss_lobatto_legendre_points(porder))*delta_mid/2

print('-----------')
print(f'Points coming out of GLL: {gauss_lobatto_legendre_points(porder)}')
print(f'Points transformed to [0 to 1] and extended to the desired length: {(1+gauss_lobatto_legendre_points(porder))*delta_x/2}')

print('-----------')
print(f"Gauss-Lobatto-Legendre points for n={porder}:")
print('x',points_x)
print('z',points_z)
print('mid',points_mid)
print('wall',points_wall)

print('-----------')
print('spacings_x (min): ',min(abs(np.diff(points_x))),    ' spacings_x (max): ',max(abs(np.diff(points_x))))
print('spacings_z (min): ',min(abs(np.diff(points_x))),    ' spacings_z (max): ',max(abs(np.diff(points_z))))
print('spacings_y (min): ',min(abs(np.diff(points_wall))), ' spacings_y (max): ',max(abs(np.diff(points_wall))))

# Open the .geo file in write mode
with open(geo_filename, "w") as geo_file:
    # Write the additional code
    geo_file.write("// Additional parameters\n")
    geo_file.write("d0 = "+str(d0)+";"+" pi = "+str(math.pi)+";"+" Lx = "+str(Lx)+";"+" Ly = "+str(Ly)+";"+" Lz = "+str(Lz)+";\n")
    geo_file.write("Nx = "+str(Nx)+";"+" Nz = "+str(Nz)+";"+" s = d0; //Main box\n")
    
    # Write the Gmsh script for the 1D line
    geo_file.write("// 1D Line\n")

    # Write points on the 1D line
    for i, yi in enumerate(y):
        geo_file.write(f"Point({i + 1}) = {{0, {yi}, 0, 1.0}};\n")

    for i in range(p+per-1):
        geo_file.write(f"Line({i + 1}) = {{{i+1}, {i+2}}};\n")

    # Connect the points to form a line
    geo_file.write("Transfinite Curve {1")
    for i in range(2, len(y)):
        geo_file.write(f", {i}")
    geo_file.write("} = 1;\n")

    # Physical line for boundary condition or grouping
    geo_file.write("// Change layer to increase x subdivision\n")
    geo_file.write("Extrude {Lx, 0, 0} { Line{1")
    for i in range(2, len(y)):
        geo_file.write(f", {i}")
    geo_file.write("}; Layers{(Nx)*d0}; Recombine;}\n")
    geo_file.write("// Change layer to increase z subdivision\n")
    geo_file.write("Extrude {0, 0, Lz} { Surface{")
    for i in range(0,len(y) - 1):
        if i == 0:
            geo_file.write(f"{80+(per*5) - 4*i}")
        else:
            geo_file.write(f",{80+(per*5) - 4*i}")
    geo_file.write("}; Layers{Nz*d0}; Recombine;}\n")

    geo_file.write("Physical Surface(\"wall\") = {")
    for i in range(1,3):
        if i == 1:
            geo_file.write(f"{93+(per*5)}")
        else:
            geo_file.write(f",{431+(per*27)}")
    geo_file.write("};\n")

    geo_file.write("Physical Surface(\"Periodic\") = {")
    for i in range(0,len(y) - 1):
        if i == 0:
            geo_file.write(f"{80+(per*5) - 4*i},{102+(per*5) + 22*i},{89+(per*5)+22*i},{97+(per*5)+22*i}")
        else:
            geo_file.write(f",{80+(per*5) - 4*i},{102+(per*5) + 22*i},{89+(per*5)+22*i},{97+(per*5)+22*i}")
    geo_file.write("};\n")
    geo_file.write("Physical Volume(\"fluid\") = {1")
    for i in range(2, len(y)):
        geo_file.write(f", {i}")
    geo_file.write("};\n")
    geo_file.write("Mesh.MshFileVersion = 2.2;\n")
    geo_file.write(f"Mesh.ElementOrder = {porder};\n")
    geo_file.write("Mesh 3;\n")
    for i in range(len(y)-1):
        geo_file.write(f"Periodic Surface  {{{102+(per*5) + 22*i}}} = {{{80+(per*5) - 4*i}}} Translate {{{0},{0},{Lz}}};\n")
    for i in range(len(y)-1):
        geo_file.write(f"Periodic Surface {{{97+(per*5)+22*i}}} = {{{89+(per*5)+22*i}}} Translate {{{Lx},{0},{0}}};\n")

print('')
print(f"Gmsh .geo file '{geo_filename}' has been created.")











