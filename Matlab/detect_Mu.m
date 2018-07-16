function mu = detect_Mu(a, B, n)
mu=0;
%Uncomment for Bournolli Dist
% if strcmp(a,'Binomial')
%     dist = makedist('Binomial', 'N', 1, 'p', 0.5);
% else
%     dist=makedist(a);
% end
dist=makedist(a);

A=zeros(n);
for i=1:n
    for j=1:n
        A(i,j)=random(dist);
    end
end

A=A*B;
for i=1:n
    for j=1:n
        temp = abs( (dot(A(i,:),A(:,j))) );
        temp = temp./[norm(A(i,:)).*norm(A(:,j))];
        mu = max( mu, temp);
    end
end

end%end function