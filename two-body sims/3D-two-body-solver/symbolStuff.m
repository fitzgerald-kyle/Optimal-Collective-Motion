%%%%% NOTE %%%%%
% YOU ARE CURRENTLY ASSUMING THAT mu2 + mu8 > 0!!!

global t
syms mu2 mu8 real
syms s t positive
%syms mu1 mu2 mu3 mu4 mu5 mu6 mu7 mu8 mu9 mu10 mu11 mu12 s ...
%    dmu1 dmu2 dmu3 dmu4 dmu5 dmu6 dmu7 dmu8 dmu9 dmu10 dmu11 dmu12 ...
%    K T t
%{
dmu1 = s/(2*s+1)*(mu2*mu9-mu3*mu8);
dmu2 = s/(2*s+1)*(mu3*mu7-mu1*mu9)+mu6;
dmu3 = s/(2*s+1)*(mu1*mu8-mu2*mu7)-mu5;
dmu4 = 1/(2*s+1)*((1+s)*mu3*mu5+s*mu9*mu5-(1+s)*mu2*mu6-s*mu8*mu6);
dmu5 = 1/(2*s+1)*((1+s)*mu1*mu6+s*mu7*mu6-(1+s)*mu3*mu4-s*mu9*mu4);
dmu6 = 1/(2*s+1)*((1+s)*mu2*mu4+s*mu8*mu4-(1+s)*mu1*mu5-s*mu7*mu5);
dmu7 = -dmu1;
dmu8 = -dmu2+mu6+mu12;
dmu9 = -dmu3-mu5-mu11;
dmu10 = 1/(2*s+1)*((1+s)*mu9*mu11+s*mu3*mu11-(1+s)*mu8*mu12-s*mu2*mu12);
dmu11 = 1/(2*s+1)*((1+s)*mu7*mu12+s*mu1*mu12-(1+s)*mu9*mu10-s*mu3*mu10);
dmu12 = 1/(2*s+1)*((1+s)*mu8*mu10+s*mu2*mu10-(1+s)*mu7*mu11-s*mu1*mu11);

K = 1/(2*s+1)*sqrt(((1+s)*mu2+s*mu8)^2+((1+s)*mu3+s*mu9)^2);
T = 1/(2*s+1)*((1+s)*mu1+s*mu7) - (((1+s)*mu6+s*mu12)*((1+s)*mu3+s*mu9) ...
    + ((1+s)*mu5+s*mu11)*((1+s)*mu2+s*mu8) + s*(mu3*mu7-mu1*mu9+mu2*mu7-mu1*mu8)) / (2*s+1)^2*K^2;

diffK = ((1+s)*mu2+s*mu8)*((1+s)*dmu2+s*dmu8) + ((1+s)*mu3+s*mu9)*((1+s)*dmu3+s*dmu9);

eqn1 = ((1+s)*mu6+s*mu12)*((1+s)*mu2+s*mu8) - ((1+s)*mu5+s*mu11)*((1+s)*mu3+s*mu9);
eqn2 = ((1+s)*mu6+s*mu12)*((1+s)*mu3+s*mu9) + ((1+s)*mu5+s*mu11)*((1+s)*mu2+s*mu8);

diffeqn1 = ((1+s)*mu6+s*mu12)*((1+s)*dmu2+s*dmu8)+((1+s)*dmu6+s*dmu12)*((1+s)*mu2+s*mu8) - ...
    ((1+s)*mu5+s*mu11)*((1+s)*dmu3+s*dmu9)-((1+s)*dmu5+s*dmu11)*((1+s)*mu3+s*mu9);
diffeqn2 = ((1+s)*mu6+s*mu12)*((1+s)*dmu3+s*dmu9)+((1+s)*dmu6+s*dmu12)*((1+s)*mu3+s*mu9) + ...
    ((1+s)*mu5+s*mu11)*((1+s)*dmu2+s*dmu8)+((1+s)*dmu5+s*dmu11)*((1+s)*mu2+s*mu8);
%}
%
mu_2 = [0 mu2 0]; mu_8 = [0 mu8 0];
syms u2 u5 real
%u2 = ((1+s)*mu2+s*mu8)/(2*s+1); u5 = ((1+s)*mu8+s*mu2)/(2*s+1);
v = [1 0 0];

syms U positive
%U = s*(mu2+mu8)/(2*s+1);

R1 = simplify(R([0 u2 0])); R2 = simplify(R([0 u5 0])); S = simplify(R([0 U 0]));
%{
M11 = simplify( S' + s/(2*s+1)*S'*int(S*Hat(mu_2),t,0,t) );
M12 = simplify( -S'*int(S*Hat(v)*R1',t,0,t) - ...
    s/(2*s+1)*S'*int(S*Hat(mu_2)*int(Hat(v)*R1',t,0,t),t,0,t) );
M13 = simplify( -S' - s/(2*s+1)*S'*int(S*Hat(mu_8),t,0,t) + eye(3) );
M14 = simplify( S'*int(S*Hat(v)*R2',t,0,t) + ...
    s/(2*s+1)*S'*int(S*Hat(mu_8)*int(Hat(v)*R2',t,0,t),t,0,t) - int(Hat(v)*R2',t,0,t) );
M31 = simplify( -S' - s/(2*s+1)*S'*int(S*Hat(mu_2),t,0,t) + eye(3) );
M32 = simplify( S'*int(S*Hat(v)*R1',t,0,t) + ...
    s/(2*s+1)*S'*int(S*Hat(mu_2)*int(Hat(v)*R1',t,0,t),t,0,t) - int(Hat(v)*R1',t,0,t) );
M33 = simplify( S' + s/(2*s+1)*S'*int(S*Hat(mu_8),t,0,t) );
M34 = simplify( -S'*int(S*Hat(v)*R2',t,0,t) - ...
    s/(2*s+1)*S'*int(S*Hat(mu_8)*int(Hat(v)*R2',t,0,t),t,0,t) );
%}
%{
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
%{
for i= 1:4 % problem elements w/ (mu2,mu8)=(2pi,2pi): J(10,4), J(10,6), others
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