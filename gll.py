import numpy as np
from scipy.special import legendre

delta_z= 37
def gauss_lobatto_legendre_points(n):
    # Calculate the roots of the derivative of the Legendre polynomial of order n
    roots, _ = np.polynomial.legendre.leggauss(n - 1)

    # Add the endpoints -1 and 1
    points = np.concatenate(([-1.0], roots, [1.0]))

    return points

# Example: Calculate Gauss-Lobatto-Legendre points for n=4
n = 8
points = (1 - gauss_lobatto_legendre_points(n))*delta_z/2

print(f"Gauss-Lobatto-Legendre points for n={n}:")
print(points)
