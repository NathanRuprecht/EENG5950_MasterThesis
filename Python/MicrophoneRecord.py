import time
import math
import numpy as np
from scipy.fftpack import dct
import pyaudio
#import msvcrt #cannot use on Pi
import pandas as pd

import threading

#from pyaudio import PyAudio, paFloat32
pa = pyaudio.PyAudio()

FORMAT = pyaudio.paFloat32
WIDTH = 2
CHANNELS = 1
DEVICE=3
Fo = 8000
TIME_TO_CALIBRATE=3

reducedFs=8e3

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
       
def inputData(in_data, frame_count, time_info, flag):
    global xo, input_buffer, data_ready
    if data_ready: print("Error: Buffer Overrun")
    input_buffer = np.matrix(np.frombuffer(in_data, dtype=np.float32))
    xo = np.append(xo,input_buffer) #saves filtered data in an array
    with threading.Lock():
        data_ready=1

    return (input_buffer, pyaudio.paContinue)

def processData(data):
    global data_done
    
    if ((data_ready) & (data_done<NUM_DB)):
        global start
        start=time.time()
        if Fo != reducedFs: downsample(data)
        
        f=Psi*np.transpose(input_buffer[:,data.start:data.end])
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
    global output_buffer, Ms, data_done, data_ready, timeC
    if data_done>=NUM_DB:
        with threading.Lock():
            data_ready=0
        
        for i in range(NUM_DB):
            output_buffer=np.append(output_buffer, db[i].buff)
            Ms=np.append(Ms, db[i].M)
        
        timeC=np.append(timeC, time.time()-start)
        with threading.Lock():
            data_done=0            

    return

def stopData():
    global stop
    
    time.sleep(timeInterval)
    #input("Press ENTER in console to STOP")
    #msvcrt.getch()
    stop=1

def calibrate(in_data, frame_count, time_info, flag):
    global n
    if flag:
        print("Calibration Error")
    ntemp = np.frombuffer(in_data, dtype=np.float32)
    n = np.append(n,ntemp) #saves filtered data in an array
    return (ntemp, pyaudio.paContinue)

def downsample(data):
    #Trim the signal so we can down sample
    data.buff=np.transpose(input_buffer)
    
    trimmer= len(data.buff) % downratio
    length=len(data.buff)-trimmer
    x=[]
    for i in range(length):
        x = np.append(x, data.buff[i])
    #Actual down sample
    Fs=Fo/downratio
    data.buff=[]
    for i in range(int(length/downratio)):
        data.buff = np.append(data.buff, x[i*downratio])
    
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
    
def main(sec, n, f):
    global Phi, Psi, stream, pa, noise, Ms, output_buffer, t, db, timeInterval, N, Fo, m, downratio, NUM_DB
    
    timeInterval=sec
    N=n
    Fo=int(f)
    
    downratio=math.floor(Fo/reducedFs)
    N = 1024
    M = 128
    SAMPLES = N*2*downratio
    NUM_DB = int(SAMPLES/N/downratio)
    m=np.log(N)*(0.35**2)

    if N==128:
        phi="Phi128.txt"
    elif N==256:
        phi="Phi256.txt"
    elif N==512:
        phi="Phi512.txt"
    elif N==1024:
        phi="Phi1024.txt"
    
    t=[]
    db=[]   
    for i in range(NUM_DB):
        db=np.append(db, data_buffer([], M, i*N, (i+1)*N))
        t=np.append(t, threading.Thread(target=processData, args=(db[i],), daemon=True))
    tout = threading.Thread(target=outputData, daemon=True)
    tstop = threading.Thread(target=stopData, daemon=True)
    
    #CS variable Phi and Psi
    I = np.matrix(np.identity(N))
    Psi=np.matrix(dct(I))
    Phi=np.matrix(np.reshape( np.loadtxt(phi),(N,N)))
    
    print("Calibrating - Stay quiet to measure noise")
    calib = pa.open(format=FORMAT,
                    channels=CHANNELS,
                    rate=Fo,
                    input_device_index=DEVICE,
                    output=False,
                    input=True,
                    frames_per_buffer = SAMPLES,
                    stream_callback=calibrate)
    calib.start_stream()
    
    while calib.is_active():
        time.sleep(TIME_TO_CALIBRATE)
        calib.stop_stream()
    calib.close()
    pa.terminate()
    pa = pyaudio.PyAudio()
    
    noise = detect_Noise(dct(np.matrix(n)))
    print("Noise: ",noise,"dB")
    
    #input("Press Enter key in console to START")
    #msvcrt.getch()
    
    stream = pa.open(format=FORMAT,
                 channels=CHANNELS,
                 rate=Fo,
                 input_device_index=DEVICE,
                 output=False,
                 input=True,
                 frames_per_buffer = SAMPLES,
                 stream_callback=inputData) #if error on "invalid" output, toggle using pa.terminate() above
     
    stream.start_stream()
    for i in range(NUM_DB):
        t[i].start()
    tout.start()
    tstop.start()
    
    
    while not stop:
        for i in range(NUM_DB):
            if not t[i].isAlive():
                t[i] = threading.Thread(target=processData, args=(db[i],), daemon=True)
                t[i].start()
        if not tout.isAlive():
            tout = threading.Thread(target=outputData, daemon=True)
            tout.start()
              
    stream.stop_stream()
    stream.close() 
    pa.terminate()
    
    for i in range(NUM_DB):
        t[i].join()
    tout.join()
    tstop.join()
    
    #Creat y matrix with output_buffer and Ms, write to txt file for Tx     
    start=len(Ms)
    end=len(output_buffer)
    if end>start:
        for i in range(end-start):
            Ms=np.append(Ms,0)
    else:
        for i in range(start-end):
            output_buffer=np.append(output_buffer,0)  
    y=pd.DataFrame(data={'Y': output_buffer, 'M': Ms})
    x=pd.DataFrame(data={'X': xo})
    np.savetxt(r'y.txt', y.values, fmt='%d')
    np.savetxt(r'x.txt', x.values, fmt='%1.4e')

#Window main parent and title
if __name__ == "__main__":
    main()
    