function noise = detect_Noise(x, T)

noisefloor=0;

switch T
    case 'DCT'
        ourlength=length(x);
        ratio=0;
        dB=2;
        db=real(20.*log(x));
        f=length(x);
    case 'FFT'
        ourlength=length(x)/2;
        ratio=0.8;
        dB=0.01;
        
        db=real(20.*log(x(1:length(x)/2)));
        f=length(x)/2;
end

for i=1:ourlength
    if i > (ratio*ourlength)
       noisefloor = noisefloor + 20*log(abs(x(i,1)));
    end
end
noise = noisefloor / ((1-ratio)*ourlength);
noise = noise+20*log(dB);
end%end function