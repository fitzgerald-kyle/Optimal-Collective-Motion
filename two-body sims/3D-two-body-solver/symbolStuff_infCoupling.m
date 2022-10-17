global t
syms mu2 real
syms s t positive

%%%%%%%%%%%%%%%%%% u5 = u2 + 2*pi,  u2 positive %%%%%%%%%%%%%%%%%%%%

%
mu8 = mu2+2*pi*(2*s+1);
mu_2 = [0 mu2 0]; mu_8 = [0 mu8 0];
syms u2 real
u5 = u2+2*pi;
v = [1 0 0];

syms U real
%U = s*(mu2+mu8)/(2*s+1);

R1 = simplify(R([0 u2 0])); R2 = simplify(R([0 u5 0])); S = simplify(R([0 U 0]));
%
M11 = simplify( S' + s/(2*s+1)*S'*int(S*Hat(mu_2),t,0,t) );
M12 = simplify( -S'*int(S*Hat(v)*R1',t,0,t) - ...
    s/(2*s+1)*S'*int(S*Hat(mu_2)*int(Hat(v)*R1',t,0,t),t,0,t) );
M13 = simplify( -S' - s/(2*s+1)*S'*int(S*Hat(mu_8),t,0,t) + eye(3) );
M14 = simplify( S'*int(S*Hat(v)*R2',t,0,t) + ...
    s/(2*s+1)*S'*int(S*Hat(mu_8)*int(Hat(v)*R2',t,0,t),t,0,t) - int(Hat(v)*R2',t,0,t) );
M31 = simplify( -M11 + eye(3) );
M32 = simplify( -M12 - int(Hat(v)*R1',t,0,t) );
M33 = simplify( -M13 + eye(3) );
M34 = simplify( -M14 - int(Hat(v)*R2',t,0,t) );
%}
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
%}
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