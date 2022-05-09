# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 10:09:22 2021

@author: Kedar
"""

"""
 G1 and G2 are the spectral splitings in GHz , G1 betw / peak a and d ,
 G2 betw / peak b and c .
 tg , pg and te , pe are freeparameters.
 tg = 0 means , all the splitting corresponds to spin−orbit splitting
 te = pi/2 means , all the splitting is due to the Jahn teller effect .
 tB is theta_B , the relative angle betw / SiV high symmetry axis and the
     external b−field . f is the factor that diminishes the orbital g−factor
 g_L ( c . f . function " Hamiltonian " at the end of this file ) .
 eps is the value of uniaxial stress in GPa , we employ a scaling factor
     of the NV center as a reference . For arbitrary stress expressed in
     parameters alpha , beta ( c f . Sec . 2 . 2 . 5 ) , set eps = alpha + i ∗ beta .

"""
"""
B−field vector , same dimension as the expermental data vector to plot
the two in the same graph
"""
from numpy import ndarray

class myarray(ndarray):    
    @property
    def H(self):
        return self.conj().T
from numpy.linalg import eig, inv
from numpy import conj, ones
import numpy as np


def SiVModel (G1, G2, tg , pg , te , pe , tB , f , eps ):
    # constant and conversions
    t=tB*np.pi/180
    p = 45*np.pi/180
    GHz = 1e9
    c = 3e8
    G1 = G1*GHz
    G2 = G2*GHz

    k=2
    eps = eps*242*k*GHz
    # The value 242 GHz/GPa is the corresponding stress
    # response of the NV center. The factor
    # k is the free parameter for the SiV center 

    # Level splitting in excited (EEe) and ground state (EEg) 
    EEe = (G1+G2)/4
    EEg = (G1+G2)/4
    
    B=np.arange(0,7.1,0.1)
    
    #ground states
    Dg = 2*EEg*np.cos(tg)
    Qgx = EEg*np.sin(tg)*np.cos(pg)
    Qgy = EEg*np.sin(tg)*np.sin(pg)
    Qg = np.sqrt(Qgx**2+Qgy**2)
    DeltaEg = 2*np.sqrt(Qg**2+(Dg**2)/4)/GHz
    
    #excited states
    De = 2*EEe*np.cos(te)
    Qex = EEe*np.sin(te)*np.cos(pe)
    Qey = EEe*np.sin(te)*np.sin(pe)
    Qe = np.sqrt(Qex**2+Qey**2)
    DeltaEe = 2*np.sqrt(Qe**2+(De**2)/4)/GHz
    print(De, Qgx, Qgy ,Dg/GHz, Qgx/GHz, Qgy/GHz,end="/n")
        
    for b in range(0,np.size(B)-1):
        # B−Field as a vector
        Bv = np.array([[np.sin(t)*np.cos(p)] , [np.sin(t)*np.sin(p)] , [np.cos (t)]])*B[b]
        # Calculate the Hamiltonian matrices , for the uniaxial stress
        # measurements of Sternschlte 1995 , the results fit best with the
        # excited state stress parameter scaled by 1.3
        Hg = Hamiltonian (Dg, Bv , Qgx , Qgy , f , eps ) 
        He = Hamiltonian (De , Bv , Qex , Qey , f , 1.3*eps ) 
        # Calculate the eigenstates and −values for ground state . . .
        [Vg, Ega] = eig(Hg) 
        Eg = []
        Eg.append(np.diag(Ega))
        # and excited state
        [Ve , Eea] = eig(He)
        Ee = []
        Ee.append(np.diag(Eea))
        # Calculate the transition frequencies ( the PEAKS in the spectrum )
        # and the peak intensities , ILHWP is for the polarization analysis
        Ta , Ia,  Axa, Aya, Aza = Transitions(Eg , Ee , Vg , Ve ) 
        T = Ta
        # T = T[b]
        I = Ia
        # I[b] = Ia 
        Ax = Axa
        Ay = Aya
        Az = Aza
        # Ax[:,b] = Axa 
        # Ay[:,b] = Aya 
        # Az[:,b] = Aza 
        # ILHWP = np.transpose(ILHWP) 
        # Trans formation matrix for eigenstates
        # transform to e+/e− states ( eigenbasis for L_z operator )
        Tv = 1/np.sqrt(2)*np.matrix([[-1,0,-1j,0],[0,-1,0,-1j],[1,0,-1j,0],[0,1,0,-1j]])
        Vg = np.matmul(inv(Tv.H),Vg)
        Ve = np.matmul(inv(Tv.H),Ve)
        
    #constants
    dw = 10e9
    NE = 1500
    minE = -1500*GHz
    maxE = 1500*GHz   
    Ef =  (np.arange(minE,maxE,(maxE-minE)/(NE-1)))
    Bf = ((np.ones(np.size(Ef))).H)*B
    E = (Ef.T)*np.ones(np.size(B))
    dimension = np.ones(np.size(Ef)).H
    Exp = dimension*np.ones(np.size(B))
    I = I/np.max(np.max(I))
    for t in range(0,np.size(T,1)):
        ExpT = dimension*3.7e03*I[t,:]/(1+(E-dimension*T[t,:])**2)/(dw**2)
        Exp = Exp + ExpT
    return Eg, Ee, T, Ef, Exp

def StressUniaxial(eps):
    I2 = np.diag([1,1])
    
    HStress2dim = eps*[[-1*np.sqrt(3)],[np.srqt(3)]]
    alpha = np.real(eps)
    beta = np.imag(eps)
    HStress2dim = [[alpha,beta],[beta,-alpha]]
    HStresU = np.kron(HStress2dim,I2)
    
def Hamiltonian(D,B,a,b,f,eps):
    ge = 28e9
    gL = ge/2
    
    Lx = [[0,1],[1,0]]
    Ly = [[0,1j],[-1j,0]]
    Lz = np.diag([1 , -1])
    
    Sx = (1/2)*np.array([[0,1],[1,0]]) 
    Sy = (1j/2)*np.array([[0,-1],[1,0]])
    Sz = (1/2)*np.array([[1,0],[0,-1]])
    
    I2 = np.diag([1,1])
    Bx = B[0]
    By = B[1]
    Bz = B[2]
    
    H = -D*np.kron(Ly,Sz) #spin orbit coupling
    H += f*gL*np.kron(Bz*Ly,I2) #Orbital Zeeman effect
    H += ge*np.kron( I2 , Bz*Sz ) # spin Zeeman effect , x−component
    H += ge*np.kron( I2 , Bx*Sx ) # spin Zeeman effect , y−component
    H += ge*np.kron( I2 , By*Sy ) # spin Zeeman effect , z−component
    H += a*np.kron( Lz , I2 ) # Jahn−Teller in x−direction
    H += b*np.kron( Lx , I2 ) # # Jahn−Teller in y−direction
    
    # if we add the uniaxial stress , uncomment the following line
    # Hstress = StressUniaxial(eps)
    # H += Hstress
    return H
def Transitions(Eg ,Ee, Vg, Ve):
    kB = 1.38065e-23/6.6261e-34
    etaX = 1 
    etaY = 0.642
    etaZ = 0.823
    
    # dipole matrix
    px = [[1,0],[0,-1]]
    py = [[0,-1],[-1,0]]
    pz = 4*np.array([[1,0],[0,1]])
    
    #spin space representation
    I2 = np.eye(2)
    Px = np.kron(px,I2)
    Py = np.kron(py,I2)
    Pz = np.kron(pz,I2)
    
    T=[]
    Tlbl=[]
    Ax=[]
    Ay=[]
    Az=[]
    Ix=[]
    Iy=[]
    Iz=[]
    I=[]
    
    
    
    for g in range(0,np.size(Eg)-1):
        for e in range(0,np.size(Ee)-1):
            t = (np.size(Ee)-1)*(g)+e
            print(np.shape(Ee),np.shape(Eg),t)
            T.append(Ee[0][e]-Eg[0][g])
            vg = Vg[g]/np.linalg.norm(Vg[g])
            ve = Ve[e]/np.linalg.norm(Ve[e])
            Ax.append(conj(vg)*Px*ve) 
            Ay.append(conj(vg)*Py*ve) 
            Az.append(conj(vg)*Pz*ve)
            print(np.shape(Ax))
            Ix.append(np.abs(Ax[t])**2)
            Iy.append(np.abs(Ay[t])**2)
            Iz.append(np.abs(Az[t])**2)
            I.append((etaX*Ix[t]+etaY*Iy[t]+etaZ*Iz[t])*np.exp(-(Ee[0][e]-Ee[0][0]))/(kB*12))
    return T, I, Ax, Ay, Az       
 
print(SiVModel(0,0,0,0,0,0,0,0,0 + 0j))