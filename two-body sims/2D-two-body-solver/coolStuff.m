figure;
n=1;
%v= linspace(-2,2,n);
%u= v.^3;
%p3 = 12.25*ones(1,n)+u;
p3= 13
for i = 1:n
    output_IVP = solve_IVP([0 0 0 0 0 0], 20*[-900 191.3 p3(i) -900 -191.3 -p3(i)], ...
        1, parameters(200), 0, 1);
    plot_function(output_IVP, output_IVP.x(end,:));
end