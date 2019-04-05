#coding: utf-8
import scipy.constants as sci
import numpy as np
import matplotlib.pyplot as pl
import math
from mpl_toolkits.mplot3d import Axes3D 

from matplotlib import cm
from numpy.linalg import inv


#Def Constants
v1 = 100
v2 = 10
v3 = 0

#matriz completa, onde (-1) são os nós livres
def fullMtz(n):
	mtz = np.empty([(n+1), (n+1)])*0 - 1
	nx = np.linspace(0,12,n+1)
	ny = np.linspace(12,0,n+1)
	for i in range (0,(n+1)):
			for j in range (0,(n+1)):
				#Quadrado Externo
				if nx[i] == 0 or nx[i] == 12:
					mtz[i,j] = v3
				if nx[j] == 0 or nx[j] == 12:
					mtz[i,j] = v3

				#Paredes Horizontais
				if ny[i] == 3 and nx[j] <= 6:
					mtz[i,j] = v3
				if ny[i] == 6 and nx[j] >=3 and nx[j] <= 9:
					mtz[i,j] = v3
				if ny[i] == 9 and nx[j] >=3:
					mtz[i,j] = v3
				#Paredes Verticais
				if nx[j] == 3 and ny[i] <= 9 and ny[i] >= 6:
					mtz[i,j] = v3
				if nx[j] == 9 and ny[i] <= 6 and ny[i] >= 3:
					mtz[i,j] = v3
				#Entrada V = 10
				if i==n and nx[j] >= 6 and nx[j] <= 9:
					mtz[i,j] = v2
				#Saída V = 100
				if j==n and ny[i] >= 9:
					mtz[i,j] = v1
	return mtz
	
#Definindo nós livres, agora fazer a matriz
def freeNodes(n,mtz):
	count = 0
	for i in range (0,(n+1)):
		for j in range (0,(n+1)):
				if mtz[i,j] not in [v1,v2,v3]:
					count = count + 1
	return count

#Aqui só conta a matriz dos nós livres
def mtzCount(n, mtz):
	count = 0
	for i in range (0,(n+1)):
		for j in range (0,(n+1)):
			if mtz[i,j] in [v1,v2,v3]:
				mtz[i,j] = 0
			if mtz[i][j] == (-1):
				count = count + 1
				mtz[i][j] = count
	return mtz
			
def move(i, j, direction):
	"""Move from index (i,j) in direction ’up’, ’down’, ’left’ or ’right’.
	"""
	if direction == 'up':
		return i - 1, j
	if direction == 'down':
		return i + 1, j
	if direction == 'left':
		return i, j - 1
	if direction == 'right':
		return i, j + 1
	# Unknown direction
	raise ValueError("Unknown direction %s" % direction)
	
#Pegar valores dos vizinhos, hora Indice, hora Potencial					
def getBrod(i, j, n, mtz):
	"""Return the k-indices of the (i,j) neighbours
	(k_up, k_right, k_down, k_left)
	going clockwise from the neighbour above node (i,j).
	"""
	klst = []
	for direction in ['up', 'left', 'down', 'right']:
		idir, jdir = move(i, j, direction)
		kdir = mtz[idir,jdir]
		klst.append(kdir)
	return klst
	
#Retorna a matriz de nós_livre x nós_livres, diagonal -4.
def tryMtz(n, mtz, fnodes):
	sparse_mtz = np.empty([fnodes, fnodes])*0
	count = 0
	#print mtz
	for i in range (0,(n+1)):
			for j in range (0,(n+1)):
				lista = []
				if mtz[i,j] != 0:
					
					lista = getBrod(i,j,n,mtz)
					
					for k in range (0,(fnodes)):
						
						if k in lista and k != 0:
							sparse_mtz[count][k-1] = 1
					sparse_mtz[count][count] = -4
					count = count + 1
					#print "Free node ", count , "Brods ", lista
	
	return sparse_mtz
#Vetor B, será o somatório negativo dos nós vizinhos. desconhecidos = 0
def vetB(n, mtz, fnodes):
	B = np.zeros((fnodes,1))
	count = 0
	for i in range (0,(n+1)):
		for j in range (0,(n+1)):
			somat = 0
			lista = []
			if mtz[i,j] == -1:
				lista = getBrod(i,j,n,mtz)
				for k in range (0, (4)):
					if lista[k] != -1:
						somat = somat - lista[k]
			
				B[count] = somat
				count = count + 1
	return B

#Finalmente o Vetor aclamado				
def vetPot (mtz, bvet):
	invMtz = inv(mtz)
	vetPot = np.dot(invMtz,bvet)
	return vetPot
	
#Mtz original
def posArray (n, mtz, fnodes):
	nx = np.linspace(0,12,n+1)
	ny = np.linspace(12,0,n+1)
	posX = np.empty([fnodes])
	posY = np.empty([fnodes])
	count = 0
	for i in range (0,(n+1)):
		for j in range (0,(n+1)):
			if mtz[i,j] == -1:
				posX[count] = nx[j]
				posY[count] = ny[i]
				count = count + 1
	
	#print posX
	#print posY
	return (posX, posY)

#plotar Surf recebe original (-1), todos os pontos
def potMtz(n, mtz, bvet):
	count = 0
	for i in range (0,(n+1)):
		for j in range (0,(n+1)):
			if mtz[i,j] == -1:
				mtz[i,j] = bvet[count]
				count = count + 1
	
	return mtz
#plotar scatter, todos os pontos
def newScat(n, MtzX, MtzY, MtzV):
	scatX = np.empty([(n+1)**2])
	scatY = np.empty([(n+1)**2])	
	scatV = np.empty([(n+1)**2])
	k = 0
	for i in range (0,(n+1)):
		for j in range (0,(n+1)):
			scatX[k] = MtzX[i,j]
			scatY[k] = MtzY[i,j]
			scatV[k] = MtzV[i,j]
			k = k + 1
	
	return scatX, scatY, scatV
	
def printV(V, nodes):
	for i in range (0,(nodes)):
		print "Nó ", i+1, "-> ", V[i], "\n"
		
def plotLab(n,teste):
	
	nx = np.linspace(0,12,n+1)
	ny = np.linspace(12,0,n+1)
	[MtzX,MtzY] = np.meshgrid(nx,ny)

	mtzOrigi = fullMtz(n)
	mtzCopia = fullMtz (n)
	nodes = freeNodes(n,mtzCopia)

	
	mtz_count = mtzCount(n,mtzCopia)
	sparse = tryMtz(n,mtzCopia,nodes)
	Bvet = vetB(n,mtzOrigi,nodes)
	V = vetPot(sparse,Bvet)
	#printV(V,nodes)
	
	posX, posY = posArray(n,mtzOrigi, nodes)
	MtzPot = potMtz(n, mtzOrigi, V)
	scatX,scatY,scatV = newScat(n,MtzX, MtzY,MtzPot)	
	
	print nodes
	fig = pl.figure()
	ax = fig.gca(projection='3d') 
	
	if teste == 1:
		surf = ax.scatter(posX,posY, V)
	if teste == 2:
		surf = ax.scatter(scatX,scatY, scatV)	
	if teste == 3:
		surf = ax.plot_surface(MtzX, MtzY,MtzPot, cmap=cm.jet,
		                   rstride=1, cstride=1)
	ax.w_xaxis.set_pane_color((0.5, 0.5, 1.0, 0.5))
	pl.show()
				
def main():
	print "Digite o valor de N, >= 8 e multiplo de 4"
	n = int(raw_input())
	print "Digite o valor de Graf (1,2,3) "
	teste = int(raw_input())

	plotLab(n,teste)

if __name__ == '__main__':
    main()
