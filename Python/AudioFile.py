import time
import math
import numpy as np
from scipy.io import wavfile
from scipy.fftpack import dct
import pandas as pd

import threading

reducedFs=8e3
M=1
N=1024

Phi=[]
Psi=[]

xo = [] #Original signal
n = [] #Calibration signal for noise
input_buffer = []
output_buffer = []
Ms=[]
data_ready=0
data_done=0
stop=0
timeC=[]


class data_buffer:
    def __init__(self, buff, M, start, end):
        self.buff = buff
        self.M = M
        self.start = start
        self.end = end

def processData(data):
    global data_done
    
    if data_done<NUM_DB:
        global start
        start=time.time()
        
        f=Psi*x[data.start:data.end, :]
        K = detect_K(f)
        M = int(round(m*K))
        if M>N: M=N
        data.M=M
        
        Phitemp=np.matrix(Phi[:data.M,:N])
        data.buff=np.asarray(np.matrix(Phitemp*f))
        
        with threading.Lock():
            data_done = data_done + 1
    
    return

def outputData():
    global output_buffer, Ms, data_done, data_ready, timeC, stop
    if data_done>=NUM_DB:
        with threading.Lock():
            data_ready=0
        
        for i in range(NUM_DB):
            output_buffer=np.append(output_buffer, db[i].buff)
            Ms=np.append(Ms, db[i].M)
        
        timeC=np.append(timeC, time.time()-start)
        with threading.Lock():
            data_done=0
            stop = 1
            
    return
    
def detect_Noise(f):
    length=len(f)
    noise=np.zeros((length,1,))

    for i in range(int(length)):
        noise[i]=20*np.log(np.fabs(f[i,0]))
        if noise[i]==float('-inf'):
            noise[i]=0.0
            
    return (np.mean(noise))

def detect_K(f):
    K=0
    length=len(f)

    for i in range(int(length)):
        if( 20*np.log(np.fabs(f[i,0])) > noise ):
            K=K+1
          
    return K
    
def main(filename, n):
    global Phi, Psi, stream, pa, noise, Ms, output_buffer, t, db, x, NUM_DB, m, N
    
    N=n
    m=np.log(N)*(0.35**2)

    if N==128:
        phi="Phi128.txt"
    elif N==256:
        phi="Phi256.txt"
    elif N==512:
        phi="Phi512.txt"
    elif N==1024:
        phi="Phi1024.txt"

    [Fo, xo]=wavfile.read(filename)
    xo=np.transpose(np.matrix(xo[:,0]))
    print("File\t",filename)

    #Trim the signal so we can down sample
    downratio=math.floor(Fo/reducedFs)
    trimmer= len(xo) % downratio
    length=len(xo)-trimmer
    x=np.zeros((length,1,))
    for i in range(len(x)):
        x[i] = xo[i]

    #Actual down sample
    Fs=Fo/downratio
    x1=np.zeros((int(len(x)/downratio),1))
    for i in range(len(x1)):
        x1[i] = x[i*downratio]
    x=x1
    x= np.divide( x, max(np.fabs(x)) )
    NUM_DB = int(len(x)/N)
    
    t=[]
    db=[]   
    for i in range(NUM_DB):
        db=np.append(db, data_buffer([], M, i*N, (i+1)*N))
        t=np.append(t, threading.Thread(target=processData, args=(db[i],), daemon=True))
    tout = threading.Thread(target=outputData, daemon=True)
    
    #CS variable Phi and Psi
    I = np.matrix(np.identity(N))
    Psi=np.matrix(dct(I))
    Phi=np.matrix(np.reshape( np.loadtxt(phi),(N,N)))
    
    noise = detect_Noise(dct(np.matrix(x)))
    print("Noise: ",noise,"dB")
    
    for i in range(NUM_DB):
        t[i].start()
    tout.start()    
    
    while not stop:
        for i in range(NUM_DB):
            if not t[i].isAlive():
                t[i] = threading.Thread(target=processData, args=(db[i],), daemon=True)
                t[i].start()
        if not tout.isAlive():
            tout = threading.Thread(target=outputData, daemon=True)
            tout.start()
    
    for i in range(NUM_DB):
        t[i].join()
    tout.join()    
    
    start=len(Ms)
    end=len(output_buffer)
    if end>start:
        for i in range(end-start):
            Ms=np.append(Ms,0)
    else:
        for i in range(start-end):
            output_buffer=np.append(output_buffer,0)  
    y=pd.DataFrame(data={'Y': output_buffer, 'M': Ms})
    np.savetxt(r'y.txt', y.values, fmt='%d')
    np.savetxt(r'x.txt', xo, fmt='%1.4e')


#Window main parent and title
if __name__ == "__main__":
    main()
    