import math
import sys, os
import json
import shutil
import numpy as np
from copy import deepcopy
import csv
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})


U_inf = 1
L_ref = 1
Re = 10000

point_x = []
point_y = []
vel_x = []



with open('test.csv', 'r')  as csv_file:
    csv_reader = csv.DictReader(csv_file)
    for col in csv_reader:
        point_x.append(col['Points:0'])
        point_y.append(col['Points:1'])
        vel_x.append(col['velocity_average:0'])


point_x = np.asarray(point_x)   
point_x = point_x.astype(np.float64)
point_x[abs(point_x)<=10e-7] = 0

point_y = np.asarray(point_y) 
point_y = point_y.astype(np.float64)  

vel_x = np.asarray(vel_x)
vel_x = vel_x.astype(np.float64)   

x_position = np.unique(point_x)
x_position = np.sort(x_position)

delta_star_vec = np.zeros(x_position.size)
theta_vec = np.zeros(x_position.size)
delta_99_vec = np.zeros(x_position.size)

i = 0
for x in x_position:
    if x==x_position[223]:
        print('Test')

    index_position = point_x==x
    index_position_sort = np.argsort(point_y[index_position])
    vel_x_local = vel_x[index_position]
    index_vel_x_above_99 =  vel_x_local >= 0.99*U_inf
    index_vel_x_bellow_99 =  vel_x_local < 0.99*U_inf
    point_y_local = point_y[index_position]
    y_above_99 = point_y_local[index_vel_x_above_99]
    y_bellow_99 = point_y_local[index_vel_x_bellow_99]
    delta_99_above = min(y_above_99)
    delta_99_bellow = max(y_bellow_99)
    delta_99_vec[i] = (delta_99_above + delta_99_bellow)/2
    delta_star = np.trapz(1-(vel_x_local[index_position_sort]/U_inf),\
         point_y_local[index_position_sort]) 

    theta = np.trapz((vel_x_local[index_position_sort]/U_inf)*(1-(vel_x_local[index_position_sort]/U_inf)),\
         point_y_local[index_position_sort]) 

    delta_star_vec[i] = delta_star
    theta_vec[i] = theta
    i = i + 1


f1 = plt.figure(1)
plt.plot(x_position, delta_star_vec)
plt.xlim(-10,6)
plt.grid()
plt.xlabel('Position x')
plt.ylabel(r'$\delta^*$')
plt.title(r'$\delta^*(x)$')

f2 = plt.figure(2)
plt.plot(x_position, theta_vec)
plt.xlim(-10,6)
plt.grid()
plt.xlabel('Position x')
plt.ylabel(r'$\theta$')
plt.title(r'$\theta(x)$')


nu = U_inf*L_ref/Re
Re_delta_star = U_inf*delta_star_vec/nu
Re_theta = U_inf*theta_vec/nu

f3 = plt.figure(3)
plt.grid()
plt.plot(x_position, Re_delta_star)
plt.xlim(-10,6)
plt.xlabel('Position x')
plt.ylabel(r'$Re_{\delta^*}$')
plt.title(r'$Re_{\delta^*}(x)$')

f4 = plt.figure(4)
plt.grid()
plt.plot(x_position, Re_theta)
plt.xlim(-10,6)
plt.xlabel('Position x')
plt.ylabel(r'$Re_{\theta}$')
plt.title(r'$Re_{\theta}(x)$')

f5 = plt.figure(5)
plt.grid()
plt.plot(x_position, delta_99_vec)
plt.xlim(-10,6)
plt.xlabel('Position x')
plt.ylabel(r'$\delta_{99}$')
plt.title(r'$\delta_{99}(x)$')

plt.show()


