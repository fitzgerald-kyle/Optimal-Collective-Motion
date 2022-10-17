syms F [6 6] 
syms G [6 6]
syms H [6 6]
syms mu(t) [6 1]
syms s

mu = mu(t);

C = zeros(6,6,6);
C(3,1,2) = 1;
C(1,3,2) = -1;
C(3,2,1) = -1;
C(2,3,1) = 1;
C(6,4,5) = 1;
C(4,6,5) = -1;
C(6,5,4) = -1;
C(5,6,4) = 1;

u1 = ((1+s)*mu(3)+s*mu(6))/(2*s+1);
u2 = (s*mu(3)+(1+s)*mu(6))/(2*s+1);

h = mu(1) + mu(3)*u1 + mu(4) + mu(6)*u2 - 1/2*(u1^2+u2^2+s*(u1-u2)^2);

for i= 1:6
    for j= 1:6
        sumF = 0;
        sumH = 0;
        for r = 1:6
            sumH = sumH + C(r,j,i)*functionalDerivative(h, mu(r));
            for s = 1:6
                sumF = sumF + C(i,r,s)*functionalDerivative(h, mu(r))*mu(s);
            end
        end
        F(i,j) = simplify(-functionalDerivative(sumF, mu(j)));
        G(i,j) = simplify(functionalDerivative(functionalDerivative(h, mu(i)), mu(j)));
        H(i,j) = simplify(-sumH);
    end
end