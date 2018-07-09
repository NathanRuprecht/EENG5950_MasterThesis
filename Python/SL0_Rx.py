# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import time
from scipy.fftpack import dct, idct
from scipy.io import wavfile
from decimal import Decimal
import matplotlib.pyplot as plt

#Variables needed for SL0
L=3
sigma_min=1e-5
sigma_decrease_factor=0.5
mu_0=2
Fs=8e3
phi="Phi1024.txt"

def main(): 
    #Import Phi from saved .txt file
    f=open(phi,"r")
    data=f.readlines()
    f.close
    data=[x.strip() for x in data]
    N=int(np.sqrt(len(data)))
    Phi=np.zeros((N**2,1,))
    for i in range(len(data)):
        Phi[i]=data[i]
    Phi=np.matrix(np.reshape(Phi,(N,N)))
    
    #Create inverse of Psi
    invPsi=np.matrix(idct(np.identity(N)))
    
    #Read y for observations and M used
    f=open("y.txt", "r")
    data=f.readlines()
    f.close  
    data=[x.strip() for x in data]
    data=[x.split() for x in data]
    length=len(data)
    y=np.zeros((length,1,))
    M=np.zeros((length,1,))
    for i in range(length):
        y[i]=data[i][1]
        M[i]=data[i][0]
    y=np.matrix(y)
    M=np.matrix(M)
    
    #Read x for RMS calculations
    f=open("x.txt", "r")
    data=f.readlines()
    f.close
    data=[x.strip() for x in data]
    length=len(data)
    x=np.zeros((length,1,))
    xp=np.zeros((length,1,))
    for i in range(length):
        x[i]=data[i]
    x=np.matrix(x)
    length=len(x)
    dct_test=dct(x)
    
    #Initialize variables for while loop and reconstruct xp
    indx=0
    indy=0
    j=0
    do=1
    p=np.zeros((length,1,))
    while do:
        if (indx+N) > length:
            N = length-indx
            invPsi = np.matrix(idct(np.identity(N)))
        
        start=time.time()
        ytemp = y[indy:indy+int(M[j])]
        indy = indy+int(M[j])
        
        Phitemp = Phi[:int(M[j]),:N]
        xp[indx:indx+N] = SL0(Phitemp, ytemp, sigma_min, sigma_decrease_factor, mu_0, L)
        xp[indx:indx+N] = invPsi*xp[indx:indx+N]
        p[j]=time.time()-start
        indx = indx+N
        
        j=j+1
        
        if (np.sum(M[j+1:len(M)])==0.0): do=0
    x=np.divide(x, max(np.fabs(x)))
    xp= np.divide( xp, max(np.fabs(xp)) )
    add=np.sum(p[:j])
    RMS = np.sqrt( np.mean( np.square(x-xp)))
    
    wavfile.write('xp.wav',int(Fs), np.asarray(xp, dtype=np.int16))
    #xp=pd.DataFrame(data={'X': xp})
    #np.savetxt(r'xp.txt', xp.values, fmt='%1.4e')

def SL0(A, x, sigma_min, sigma_decrease_factor, mu_0, L):
    A_pinv=np.linalg.pinv(A)

    s=np.matrix(A_pinv*x)
    sigma = 2*max(np.fabs(s))

    while sigma>sigma_min:
        for i in range(L):
            delta = np.matrix(OurDelta(s,sigma))
            s = s-mu_0*delta
            s = s- A_pinv*(A*s-x)

        sigma = sigma * sigma_decrease_factor

    return s

def OurDelta(s, sigma):
    delta=np.multiply(s,np.exp(-np.square(np.fabs(s))
    / sigma**2))

    return delta

#Window main parent and title
if __name__ == "__main__":
    main()