H= sym(zeros(6,6));
syms H13 H23 H46 H56
H(1,3)= H13;
H(2,3)= H23;
H(4,6)= H46;
H(5,6)= H56;

J= sym('J',[6 6]);


disp(simplify(H*J))