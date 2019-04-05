import scipy.constants as sci
import numpy as np
import matplotlib.pyplot as pl
import math
from mpl_toolkits.mplot3d import Axes3D 

from matplotlib import cm
from numpy.linalg import inv

#define constants
V0 = 1.0
E0 = sci.epsilon_0
Pi = math.pi
#N = 10
L = 10e-2
#Delta = L/N


def posArray(L,N):
	Delta = L/N
	vet = Delta*np.linspace(0.5 , N-0.5 , N)
	return (vet,vet)

def impedanceMtz(MtzX,MtzY, L, N):
	Delta = L/N
	Zmn = np.empty([N**2, N**2])
	s = -1
	for a in range (0,N):
		for b in range (0,N):
			t = -1
			s = s+1
			for i in range (0,N):
				for j in range (0,N):
					t = t+1
					if a == i and b ==j:
						Zmn[s][t] = (Delta*math.log( math.sqrt(2) + 1) )/(Pi*E0)
					else:
						R2 = math.sqrt( (MtzX[a][b] - MtzX[i][j])**2 + 
													(MtzY[a][b]-MtzY[i][j])**2) 
						Zmn[s][t] = (Delta**2)/(4*Pi*E0*R2)
						
	return Zmn

def rhoMtz(Rhoro, N):
	MtzRho = np.empty([N,N])
	h = 0
	for i in range (0,N):
		for j in range (0,N):
			MtzRho[i][j] = Rhoro[h]
			h = h + 1
	
	return MtzRho
def totalCharge(Rhoro, L,N):
	Delta = L/N
	h = 0
	QT = 0
	for i in range (0,N**2):
		QT = QT + Rhoro[h]
		h = h + 1
	
	return QT*(Delta**2)

def plateCapacitance(QT, V0):
	return QT/V0
	
def plateMoM(V0,L,N):
	PosX,PosY = posArray(L,N)
	[MtzX,MtzY] = np.meshgrid(PosX,PosY)
	Zmn = impedanceMtz(MtzX,MtzY,L,N)
	
	V = V0*np.ones((N**2,1))

	Zmninv = inv(Zmn)
	Rhoro = np.dot(Zmninv, V)


	MtzRho = rhoMtz(Rhoro, N)
	QT = totalCharge(Rhoro,L,N)
	C = plateCapacitance(QT,V0)

	print "Capacitancia", C, "Com N = ", N, "\n"
	print "Minimo Rho", min(Rhoro), "Com N = ", N
	print "Maximo Rho", max(Rhoro), "Com N = ", N
	
	fig = pl.figure()
	ax = fig.gca(projection='3d') 

	surf = ax.plot_surface(MtzX, MtzY,MtzRho, cmap=cm.jet,
		                   rstride=1, cstride=1) #argumentos CHAVES PORRA

	ax.w_xaxis.set_pane_color((0.5, 0.5, 1.0, 0.5))
	#ax.set_xlabel('X (m)', fontsize=16, fontweight='bold')
	
	pl.show()
	
def main():
	'''print "Digite o valor da tensao na placa"
	V0 = float(raw_input())
	print "Digite o valor do comprimento da placa"
	L = float(raw_input())'''
	print "Digite o valor das divisoes N"
	N = int(raw_input())
	plateMoM(V0,L,N)

if __name__ == '__main__':
    main()


