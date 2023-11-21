import numpy as np
from scipy.special import legendre
import math

Re=180
p = 17
per = 13
theta = np.linspace(0, np.pi, p + per)
y = -np.cos(theta)+1
Nx = 15
Nz = 25
Lx = math.pi
Lz = (math.pi)/2
delta_z= (Lz/Nz)*Re
delta_x= (Lx/Nx)*Re
delta_wall = y[1]*180
delta_mid = (y[15]-y[14])*180


def gauss_lobatto_legendre_points(n):
    # Calculate the roots of the derivative of the Legendre polynomial of order n
    roots, _ = np.polynomial.legendre.leggauss(n - 1)

    # Add the endpoints -1 and 1
    points = np.concatenate(([-1.0], roots, [1.0]))

    return points

# Example: Calculate Gauss-Lobatto-Legendre points for n=4
n = 8
points_x = (1 - gauss_lobatto_legendre_points(n))*delta_x/2
points_z = (1 - gauss_lobatto_legendre_points(n))*delta_z/2
points_wall = (1 - gauss_lobatto_legendre_points(n))*delta_wall/2
points_mid = (1 - gauss_lobatto_legendre_points(n))*delta_mid/2

print(f"Gauss-Lobatto-Legendre points for n={n}:")
print('x',points_x)
print('z',points_z)
print('mid',points_mid)
print('wall',points_wall)
