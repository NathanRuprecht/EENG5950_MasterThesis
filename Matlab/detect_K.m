function K = detect_K(x, noise, T)

switch T
    case 'DCT'
        ourlength=length(x);
    case 'FFT'
        ourlength=length(x)/2;
end

K=0;
for i=1:ourlength
    if (20*log(abs(x(i,1)))) > noise
        K = K+1;
    end
end

end%end function