import numpy as np
import os
import sys

fileName        =       sys.argv[1]
#geofile		= 	fileName+".geo.dat"
geofile		= 	fileName+".msh"
coord		=	fileName+"Mesh.coord"

f=open(geofile,'r')
g=open(coord,'w')


line = True
inCoord = False

while line:
	line = f.readline()
	
	if inCoord:
		#if 'END' in line:
		if '$EndNodes' in line:
			print('--|| ALYA : DONE.')
		#	print(line)
			inCoord = False
			line = False
		else:
			g.write(line)

	if line:
		#if 'COORD' in line:
		if '$Nodes' in line:
			line = f.readline()
			inCoord = True
			print('--|| ALYA : START WRITING')
		#	print(line)


f.close()
g.close()

