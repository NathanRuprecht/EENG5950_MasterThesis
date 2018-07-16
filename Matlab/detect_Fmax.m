function true_fmax = detect_Fmax(x, Fs)

totalEnergy=0;
ourlength=floor(length(x)/2);
for i=1:ourlength
    totalEnergy = totalEnergy+(abs(x(i,1))).^2;
end
cutoff=0.05*totalEnergy;

Energy=0;
index=0;
for i=ourlength:-1:1
    if Energy <= cutoff
        Energy=Energy+(abs(x(i,1))).^2;
        index=i;
    else
        true_fmax = (index-1) / ( (length(x)/2) / (Fs/2) );
    end
end

end %end function