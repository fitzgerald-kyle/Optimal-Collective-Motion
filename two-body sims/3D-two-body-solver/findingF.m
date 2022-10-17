N = 3;

syms s
syms mu dmu [N 6]
syms u [N 3]
syms F [6*N 6*N]

for i=1:N
    for j=1:3
        u(i,j) = 1/(N*s+1) * (mu(i,j) + s*sum(mu(:,j)));
    end
    
    dmu(i,1:3) = cross(mu(i,1:3), u(i,:)) + cross(mu(i,4:6), [1 0 0]);
    dmu(i,4:6) = cross(mu(i,4:6), u(i,:));
end

for i=1:6*N
    for j=1:6*N
        F(i,j) = diff(dmu(ceil(i/6),mod(i-1,6)+1), mu(ceil(j/6),mod(j-1,6)+1));
    end
end

simplify((N*s+1)*F)