# EENG5950 - Master's Thesis

Taken at UNT during Fall 2017 and Spring 2018.

Faculty Advisor: Dr. Xinrong Li

#### Python Files
The "Phi..." text files are the different N-sized Phi imported and used. SL0_Tx (Node) imports MicrophoneRecord and AudioFile python files depending on input source, while SL0_Rx (Sink) reconstructs and calculated RMS error.

#### Matlab Files

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

