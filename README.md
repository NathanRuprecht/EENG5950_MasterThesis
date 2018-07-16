# EENG5950 - Master Thesis

Taken at UNT during Fall 2017 and Spring 2018.

Faculty Advisor: Dr. Xinrong Li  

### Python Files
  
Python 3.6 running on Spyder via Anaconda
  
* SL0_Tx - transmitting (leaf) node using CS on user specified source  
  * AudioFile - do CS on audio file  
  * MicrophoneRecord - do CS on microphone recording for user specified length of time  
* SL0_Rx - receiver (sink) node using SL0 to reconstruct and calculate RMS error  
* Phi text files - the different N-sized Phi imported and used  
  
Audio files not provided

### Matlab Files

Version R2016a  

* create_Phi - create an N-sized Phi and saves to text file  

* Matlab_SL0_TxRx - main file running CS Tx and Rx system using SL0 reconstruction 
  * detect_Fmax - detect maximum frequency component in source to compare our system to traditional system sampling at Nyquist Rate  
  * detect_K - detect K-sparsity of source  
  * detect_Mu - detect coherence coefficient between used Phi and Psi  
  * detect_Noise - detects noise threshold of source  
  * estimate_SNR - additional error metric signal-to-noise ratio  

* BP - Basic Pursuit reconstruction algorithm  
* SL0 - Smoothed L0 norm reconstructions algorithm  
  
Audio files not provided

