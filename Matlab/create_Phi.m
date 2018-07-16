clc
clearvars

T='DCT';
a='Binomial';
n=1024;
I=eye(n);

Mu=[];
for loop=1:1
    clc
    disp(loop);
    if strcmp(a, 'Binomial')
        dist=makedist('Binomial', 'N', 1, 'p', 0.5);
    else
        dist=makedist(a);
    end
    
    A=zeros(n);
    for i=1:n
        for j=1:n
            A(i,j)=random(dist);
        end
    end
    
    switch T
        case 'DCT'
            Psi=dct(I);
        case 'FFT'
            Psi=fft(I);
    end
    
    A=reshape(A,n^2,1);
    save('Phi1024.txt','A','-ascii');
end

%Mu=mean(Mu)