import numpy as np
import matplotlib.pyplot as plt
import sys


fileName = sys.argv[1]
uTau   = float(sys.argv[2])
nu     = float(sys.argv[3])
porder = int(sys.argv[4])
Re = uTau/nu


print ("--|| ALYA :: CHECKING THE GRID")
print ("----|| ALYA :: CALCULATIONS FOR PORDER =", porder)
print ("----|| ALYA :: CALCULATIONS FOR UTAU =", uTau)
print ("----|| ALYA :: CALCULATIONS FOR RE_TAU=", Re)
print ("----|| ALYA :: CALCULATIONS FOR VISCOSITY =", nu)

data = np.loadtxt(fileName+'Mesh.coord')
f = open('yCoord.dat','w')
g = open('xCoord.dat','w')
h = open('zCoord.dat','w')

ff = open('chan.box','w')

x0 = data[:,1]
y0 = data[:,2]
z = data[:,3]

x = x0[np.where(y0==np.amax(y0))]
y = y0[np.where(x0==np.amax(x0))]
xuni = sorted(np.unique(np.around(x,decimals=6)))
yuni = sorted(np.unique(np.around(y,decimals=6)))
zuni = sorted(np.unique(np.around(z,decimals=6)))

print ("----|| ALYA :: Nx=",len(xuni))
print ("----|| ALYA :: Ny=",len(yuni))
print ("----|| ALYA :: Nz=",len(zuni))

for i in range(len(xuni)):
  g.write("%.8f\n" % xuni[i])
g.close()
print ("----|| ALYA :: X-POINTS WRITTEN TO FILE")

for i in range(len(yuni)):
  f.write("%.8f\n" % yuni[i])
f.close()

print ("----|| ALYA :: Y-POINTS WRITTEN TO FILE")

for i in range(len(zuni)):
  h.write("%.8f\n" % zuni[i])
h.close()

print ("----|| ALYA :: Z-POINTS WRITTEN TO FILE")

ff.write("base.rea \n")
ff.write("-3 \t\t\t\t\t spatial dimension (if negative dump binary re2 and .rea file)\n")
ff.write("1  \t\t\t\t\t number of fields\n")
ff.write("#===========================================================\n")
ff.write("Box  \t\t\t\t\t any string with 1st character NOT EQUAL to C,M,Y,c,m,y\n")
ff.write("%d  %d  %d \t\t\t\t\t nelx,nely,nelz (if <0, length to be equally divided)\n" \
         %(-(len(xuni)-1),len(yuni)-1,-(len(zuni)-1)))
ff.write("0.0  1.0  1.0 \t \t x_0  x_Nelx ratio")
for i in range(len(yuni)):
  ff.write("\n %.7f " % (yuni[i]-1.0))
ff.write("\t\t\t\t\t y_0  y_1 ....  y_Nelx\n")
ff.write("0.0  1.0  1.0 \t\t\t\t\t z_0  z_Nelx ratio\n")
ff.write("P  ,P  ,W  ,W  ,P  ,P   \t\t\t\t\t  ! W,E,S,N,B,T (fixed 3CHAR format)\n")
ff.close()

print ("----|| ALYA :: GENBOX WRITTEN TO FILE")

dy = np.diff(yuni)
for i,dyi in enumerate(dy):
  if dyi < 1e-8:
    print('DELETE small dy: {}'.format(dyi))
    dy[i] = np.nan



mindy = min(dy)/porder
maxdy = max(dy)/porder

dx = (xuni[1]-xuni[0])/porder
dz = (zuni[1]-zuni[0])/porder


print ('----|| ALYA :: min(dy): {}'.format(mindy))
print ('----|| ALYA :: max(dy): {}'.format(maxdy))
print ('----|| ALYA :: dx: {}'.format(dx))
print ('----|| ALYA :: dz: {}'.format(dz))

# y+ = (y*u_tau)/nu
# Re_tau = (h*u_tau)/nu
# y+ = Re_tau * y / h
print ('----|| ALYA :: min(y+): {}'.format(mindy*uTau/nu)) 
print ('----|| ALYA :: max(y+): {}'.format(maxdy*uTau/nu)) 
print ('----|| ALYA :: x+: {}'.format(dx*uTau/nu)) 
print ('----|| ALYA :: z+: {}'.format(dz*uTau/nu)) 
print ('--|| ALYA :: END GRID CHECK') 
