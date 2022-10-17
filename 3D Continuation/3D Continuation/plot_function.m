function plot_function(output,qf)

% This function creates a plot of the solution contained in the variable
% output
% Function Inputs:
%   output : structure containing the current solution
%   qf : desired boundary condition for q(tf)
% Function Outputs:
%   none

% Plot solution
x = reshape(output.q(1,4,:),1,length(output.t));
y = reshape(output.q(2,4,:),1,length(output.t));
z = reshape(output.q(3,4,:),1,length(output.t));
plot3(x,y,z,'b-')

hold on

% Plot initial condition for x at t=0
x0 = output.q(1,4,1);
y0 = output.q(2,4,1);
zf = output.q(3,4,1);
plot3(x0,y0,zf,'go')

% Plot desired boundary condition for x at t=tf
xf = qf(1,4);
yf = qf(2,4);
zf = qf(3,4);
plot3(xf,yf,zf,'ro')

hold off

% Set axis limits
axis([-0.5 1 -1 1 -1 1])

% Axis labels
xlabel('x')
ylabel('y')
zlabel('z')

% Set axis aspect ratio
daspect([1 1 1])

% Set font size
set(gca,'color','w','FontSize',20)

% Set plot frame color to white
set(gcf,'color','w')

% Draw plot now
drawnow

end