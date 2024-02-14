import numpy as np
from scipy.special import legendre
import math

Re=180
p = 17
per = 2
porder = 5
#theta = np.linspace(0, np.pi, p + per)
#y = -np.cos(theta)+1
# Wall-normal coordinates
beta = 3.8 #3.8  # Stretching factor
ny=p+per
y = np.tanh(beta * (np.linspace(0, 1, ny) - 0.5)) / np.tanh(beta * 0.5) + 1.0

Nx = 23
Nz = 25
Lx = math.pi
Lz = (math.pi)/2
delta_z= (Lz/Nz)*Re
delta_x= (Lx/Nx)*Re
delta_wall = (0.04)*180
delta_mid = (2/265)*180


def gauss_lobatto_legendre_points(porder):
    # Calculate the roots of the derivative of the Legendre polynomial of order n
    roots, _ = np.polynomial.legendre.leggauss(porder - 1)

    # Add the endpoints -1 and 1
    points = np.concatenate(([-1.0], roots, [1.0]))

    return points

# Example: Calculate Gauss-Lobatto-Legendre points for n=4
points_x = (1 - gauss_lobatto_legendre_points(porder))*delta_x/2
points_z = (1 - gauss_lobatto_legendre_points(porder))*delta_z/2
points_wall = (1 - gauss_lobatto_legendre_points(porder))*delta_wall/2
points_mid = (1 - gauss_lobatto_legendre_points(porder))*delta_mid/2

print(f"Gauss-Lobatto-Legendre points for n={porder}:")
print('x',points_x)
print('z',points_z)
print('mid',points_mid)
print('wall',points_wall)
print('-----------')
print('spacings_x (min): ',min(abs(np.diff(points_x))),    ' spacings_x (max): ',max(abs(np.diff(points_x))))
print('spacings_z (min): ',min(abs(np.diff(points_x))),    ' spacings_z (max): ',max(abs(np.diff(points_z))))
print('spacings_mid (min): ',min(abs(np.diff(points_mid))), ' spacings_y (max): ',max(abs(np.diff(points_mid))))
print('spacings_y (min): ',min(abs(np.diff(points_wall))), ' spacings_y (max): ',max(abs(np.diff(points_wall))))