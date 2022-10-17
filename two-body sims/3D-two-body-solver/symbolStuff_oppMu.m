global t
syms u s t positive

%%%%%%%%%%%%%%%%%% u2 = -u5 %%%%%%%%%%%%%%%%%%%%

%
mu = u*(2*s+1);
MU = [0 mu 0];
v = [1 0 0];

R1 = simplify(R([0 mu/(2*s+1) 0])); R2 = simplify(R([0 -mu/(2*s+1) 0]));

M11 = simplify( s/(2*s+1)*t*Hat(MU) + eye(3) );
M12 = simplify( -s/(2*s+1)*Hat(MU)*int(int(Hat(v)*R1',t,0,t),t,0,t) ...
     - int(Hat(v)*R1',t,0,t) );
M13 = simplify( s/(2*s+1)*t*Hat(MU) );
M14 = simplify( -s/(2*s+1)*Hat(MU)*Hat(v)*int(int(R2',t,0,t),t,0,t) );
M31 = simplify( -s/(2*s+1)*t*Hat(MU) );
M32 = simplify( s/(2*s+1)*Hat(MU)*Hat(v)*int(int(R1',t,0,t),t,0,t) );
M33 = simplify( -s/(2*s+1)*t*Hat(MU) + eye(3) );
M34 = simplify( s/(2*s+1)*Hat(MU)*Hat(v)*int(int(R2',t,0,t),t,0,t) ...
    -Hat(v)*int(R2',t,0,t) );
%
%
M1 = simplify( (1+s)*[M11 M12 M13 M14]+s*[M31 M32 M33 M34] );
M3 = simplify( s*[M11 M12 M13 M14]+(1+s)*[M31 M32 M33 M34] );

syms J [12 12] 
syms I [6 12]

for i= 1:4
    I(1:3,3*i-2:3*i) = simplify( R1'*int(R1*M1(:,3*i-2:3*i),t,0,t) );
    disp(['1' num2str(i)])
    I(4:6,3*i-2:3*i) = simplify( R2'*int(R2*M3(:,3*i-2:3*i),t,0,t) );
    disp(['2' num2str(i)])
end
%
%
for i= 1:4
    J(1:3,3*i-2:3*i) = simplify( 1/(2*s+1) * I(1:3,3*i-2:3*i) );
    disp(['1' num2str(i)])
    J(4:6,3*i-2:3*i) = simplify( -1/(2*s+1) * R1'*int(R1*Hat(v)*I(1:3,3*i-2:3*i),t,0,t) );
    disp(['2' num2str(i)])
    J(7:9,3*i-2:3*i) = simplify( 1/(2*s+1) * I(4:6,3*i-2:3*i) );
    disp(['3' num2str(i)])
    J(10:12,3*i-2:3*i) = simplify( -1/(2*s+1) * R2'*int(R2*Hat(v)*I(4:6,3*i-2:3*i),t,0,t) );
    disp(['4' num2str(i)])
end
%}

Jat1 = simplify(subs(J, t, 1));
%Jat1_sInf = simplify(taylor(Jat1_with_u, s, inf, 'Order', 2));
det_Jat1 = simplify(det(Jat1));

function Rans = R(u)
global t
    normu = norm(u);
    Rans = eye(3) + Hat(u)/normu*sin(normu*t) + Hat(u)^2/normu^2*(1-cos(normu*t));
end

function uhat = Hat(u)
% This function computes the angular rotation matrix corresponding to
% a 3-dimensional vector u.
% Function Inputs:
%   u : 1x3 or 3x1 vector
% Function Outputs:
%   uhat : 3x3 angular rotation matrix corresponding to u

uhat = [0     -u(3) u(2);...
     u(3)  0     -u(1);...
     -u(2) u(1)  0];

end