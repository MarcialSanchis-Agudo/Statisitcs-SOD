import numpy as np
import matplotlib.pyplot as plt
import math

# Define the points on the 1D line

p = 17
per = 13
theta = np.linspace(0, np.pi, p + per)
y = -np.cos(theta)+1
Lx = 2*(math.pi)
Lz = (math.pi)
print(range(p))
# Plot the points
plt.plot(np.zeros_like(y), y, marker='o', linestyle='-', color='b', label='1D Line Points')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('1D Line Points')
plt.legend()
plt.grid(True)
plt.show()

# Specify the filename for the Gmsh .geo file
geo_filename = "line.geo"

# Open the .geo file in write mode
with open(geo_filename, "w") as geo_file:
    # Write the additional code
    geo_file.write("// Additional parameters\n")
    geo_file.write("d0 = 1.0; pi = 3.14159265359; Lx = 2*pi; Ly = 2; Lz = (pi)*d0;\n")
    geo_file.write("Nx = 15; Nz = 25; s = d0; //Main box\n")
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
    geo_file.write("Mesh.ElementOrder = 8;\n")
    geo_file.write("Mesh 3;\n")
    for i in range(len(y)-1):
        geo_file.write(f"Periodic Surface  {{{102+(per*5) + 22*i}}} = {{{80+(per*5) - 4*i}}} Translate {{{0},{0},{Lz}}};\n")
    for i in range(len(y)-1):
        geo_file.write(f"Periodic Surface {{{97+(per*5)+22*i}}} = {{{89+(per*5)+22*i}}} Translate {{{Lx},{0},{0}}};\n")
print(f"Gmsh .geo file '{geo_filename}' has been created.")
