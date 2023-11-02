import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

import scipy.special as special



# Define numerical constants
Rtau = 180
Re =  np.exp((1.0/0.88)*np.log(Rtau/0.09))
mu = 2/Re
utau = (Rtau*mu)  # Replace with the actual value
#y_max = (utau*2)/mu # Replace with your desired 'y' value
y_max = 2
# Define the function in a Python-friendly form
def f(yp):
    term1 = utau * ((1.0 / 0.41) * np.log(1.0 + 0.41 * yp))
    term2 = 7.8 * (1.0 - np.exp(-yp / 11.0) - (yp / 11.0) * np.exp(-yp / 3.0))
    return term1 + term2

def l(yp):
    term = 1.5*(-(yp**2) + 2*yp)
    #term = ((1.5)/Rtau)*(-(yp*yp)/Rtau + 2*yp)
    return term

result = integrate.quad(lambda x: l(x), 0, 2*Rtau)
umax = Rtau/(result[0])
print(result)

# Define the range of y-values from 0 to y_max
y_values = np.linspace(0, y_max, 1000)  # Adjust the number of points as needed

# Calculate the corresponding function values
function_values = l(y_values)

# Create the plot
plt.plot(y_values, function_values)  # Add a description for your function
plt.xlabel('y')
plt.ylabel('V_l(y)')
plt.title('Parabolic laminar profile')
plt.grid(True)
plt.show()
