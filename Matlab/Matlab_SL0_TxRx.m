clearvars
close all

%%%%%%%%%%GLOBAL VARIABLES%%%%%%%%%%%%%%%%%%%
% needed for SL0
L = 3;
sigma_min = 1e-5;
sigma_decrease_factor = 0.5;
mu_0 = 2;

% needed for recording
T='DCT';
source='file'; %'file' or 'record'
filename='ClarityClip.wav';
n=1024;
phi='Phi1024.txt';
reducedFs = 8e3; %Target sample frequency to downsample sig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%CREATE ORIGINAL SIGNAL%%%%%%%%%%%%%%%
switch source
    case 'record'
        Fo=48e3;
        timeInterval=5;
        BperSamp=16;
        numChannels=1;
        
        %Record your voice for timeInterval seconds.
        recObj = audiorecorder(Fo, BperSamp, numChannels);
        disp('Start Recording.')
        recordblocking(recObj, timeInterval);
        disp('End of Recording.');
        xo = getaudiodata(recObj);
    case 'file'
        [xo, Fo] = audioread(filename);
        xo=xo(:,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N     = 250;    % Order
Fpass = 3250;   % Passband Frequency
Fstop = 4250;  % Stopband Frequency
Wpass = 1;      % Passband Weight
Wstop = 1;      % Stopband Weight
dens  = 20;     % Density Factor

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, [0 Fpass Fstop Fo/2]/(Fo/2), [1 1 0 0], [Wpass Wstop], ...
           {dens});
xo=filter(b, 1, xo);

downratio=floor(Fo/reducedFs);
trimmer=mod(length(xo),downratio);
len=length(xo)-trimmer;
x=zeros(len,1);
for i=1:length(x)
    x(i,1)=xo(i,1);
end

Fs=Fo/downratio;
x1=zeros(length(x)/downratio,1);
for i=1:length(x1)
    x1(i,1)=x(i*downratio,1);
end
x=x1;
%%%%%%%%%%PROCESSING RESULTS%%%%%%%%%%%%%%%%%
x=x./max(abs(x));

t=(1:length(x))./Fs; % define a time vector
ssf=(ceil(-length(x)/2):ceil(length(x)/2)-1)/(length(x)./Fs); % frequency vector
switch T
    case 'DCT'
        fx=dct(x(1:length(x))); % do DFT/FFT
    case 'FFT'
        fx=fft(x(1:length(x))); % do DFT/FFT        
end

% needed for recording
freqz(x);
Ys=detect_Fmax(fft(x(1:length(x))), Fs);
NyquistRate=round(2*Ys*n/Fs);
noise = detect_Noise(fx, T);

I=eye(n);
Phi='Normal';
switch T
    case 'DCT'
        Psi1=dct(128);
        Psi=dct(I);
    case 'FFT'
        Psi1=fft(128);
        Psi=fft(I);
        
end
mu=detect_Mu(Phi, Psi1, 128);
invPsi=inv(Psi);
A=importdata(phi);
A=reshape(A,n,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j=0;
mTotal=0;
KTotal=0;
timeR=0;
timeC=0;
while j < length(x)-n %go from 0 to length of song in seconds
    clc
    fprintf('%f%%\n',j/length(x)*100);
    xtemp = x(j+1:j+n);

    tic
    switch T
        case 'DCT'
            f=dct(xtemp(1:length(xtemp))); % do DFT/FFT
        case 'FFT'
            f=fft(xtemp(1:length(xtemp))); % do DFT/FFT        
    end
    K = detect_K(f, noise, T);
    KTotal=KTotal+K;
    m_form=round(mu^2*K*log(length(xtemp)));
    if strcmp(T,'FFT')
        m=round(m_form/2);
    elseif strcmp(T,'DCT')
        m=m_form;
    end
    if m>n
        m=n;
    end
    mTotal=mTotal+m;

    Phi=A(1:m,1:n);
    ytemp=Phi*Psi*xtemp;
    timeC=timeC+toc;

    tic;
    A_pinv = pinv(Phi);
    xp(j+1:j+n, 1)= SL0(Phi, ytemp, sigma_min, sigma_decrease_factor, mu_0, L, A_pinv);
    xp(j+1:j+n, 1) = invPsi*xp(j+1:j+n, 1);
    timeR=timeR+toc;

    j=j+length(xtemp);
end %end while loop

len=length(x);
xtemp = x(j+1:len);

tic
switch T
    case 'DCT'
        f=dct(xtemp(1:length(xtemp))); % do DFT/FFT
    case 'FFT'
        f=fft(xtemp(1:length(xtemp))); % do DFT/FFT        
end
K = detect_K(f, noise, T);
KTotal=KTotal+K;
m_form=round(0.35^2*K*log(length(xtemp)));
if strcmp(T,'FFT')
    m=round(m_form/2);
elseif strcmp(T,'DCT')
    m=m_form;
end
mTotal=mTotal+m;
Phi=A(1:m,1:(len-j));
I=eye(len-j);
switch T
    case 'DCT'
        Psi=dct(I);
    case 'FFT'
        Psi=fft(I);
end

invPsi=inv(Psi);
ytemp=Phi*Psi*xtemp;
timeC=timeC+toc;

tic
A_pinv = pinv(Phi);
xp(j+1:len, 1)= SL0(Phi, ytemp, sigma_min, sigma_decrease_factor, mu_0, L, A_pinv);
xp(j+1:len, 1) = invPsi*xp(j+1:len, 1);
timeR=timeR+toc;

xp=xp./max(abs(xp));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,1), plot(t,x)                  % plot the waveform
title('x in time');
xlabel('seconds'); ylabel('amplitude')     % label the axes
switch T
    case 'DCT'
        fxs=fx;
        ssf=0:1:length(x)-1;
    case 'FFT'
        fxs=fftshift(fx); % shift it for plotting
end
subplot(2,2,2), plot(ssf,abs(fxs))         % plot magnitude spectrum
title('x in freq');
xlabel('frequency'); ylabel('magnitude')   % label the axes

switch T
    case 'DCT'
        fxp=dct(xp(1:length(xp)));
        fxps=fxp;
    case 'FFT'
        fxp=fft(xp(1:length(xp)));
        fxps=fftshift(fxp); % shift it for plotting
end
subplot(2,2,3), plot(t,real(xp)) % plot the waveform
title('xp in time');
xlabel('seconds'); ylabel('amplitude') % label the axes
subplot(2,2,4), plot(ssf,abs(fxps)) % plot magnitude spectrum
title('xp in freq');
xlabel('frequency'); ylabel('magnitude') % label the axes

if strcmp(source, 'record')
    filename='Voice';
end
RMStot = sqrt( sum( (x-xp).^2 ) /length(x) );
fprintf('%s\t%s\t%s\t%f\t%f\t%f\t%i\t%f\n',...
    'Matlab', 'PC', filename, RMStot, timeC, timeR, n, mTotal/length(x)*100)
%fprintf('Audio Source\tLength(x) in sec\tCompress Time(s)\tReconstruct Time(s)\tTransform\tRMS\tFmax\tNyquist M\tAvg M Used per Frame\tMu\tOverall CR(%%)\tK\tN used\n%s\t%i\t%i\t%i\t%s\t%f\t%i\t%i\t%i\t%f\t%i\t%i\t%i\n',...
    %filename, length(x)/Fs, timeC, timeR, T, RMStot, Ys, NyquistRate, mTotal/(length(x)/n), mu, mTotal/length(x)*100,KTotal/(length(x)/n), n);

%
% function SNR=estimate_SNR(estim_s,true_s)
% 
% err = true_s - estim_s;
% SNR = 10*log10(sum(abs(true_s).^2)/sum(abs(err).^2));
% end
%close all