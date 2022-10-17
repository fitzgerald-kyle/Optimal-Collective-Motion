function plot_function(output,qf)

% This function creates a plot of the solution contained in the variable
% output
% Function Inputs:
%   output : structure containing the current solution
%   qf : desired boundary condition for q(tf)
% Function Outputs:
%   none

% Plot solution
x1 = reshape(output.q(1,4,:),1,length(output.t));
y1 = reshape(output.q(2,4,:),1,length(output.t));
z1 = reshape(output.q(3,4,:),1,length(output.t));
plot3(x1,y1,z1,'b-')

hold on

x2 = reshape(output.q(5,8,:),1,length(output.t));
y2 = reshape(output.q(6,8,:),1,length(output.t));
z2 = reshape(output.q(7,8,:),1,length(output.t));
plot3(x2,y2,z2,'c-')

% Plot initial condition for x at t=0
x01 = output.q(1,4,1);
y01 = output.q(2,4,1);
z01 = output.q(3,4,1);
plot3(x01,y01,z01,'go')

x02 = output.q(5,8,1);
y02 = output.q(6,8,1);
z02 = output.q(7,8,1);
plot3(x02,y02,z02,'go')

% Plot desired boundary condition for x at t=tf
xf1 = qf(1,4);
yf1 = qf(2,4);
zf1 = qf(3,4);
plot3(xf1,yf1,zf1,'ro')

xf2 = qf(5,8);
yf2 = qf(6,8);
zf2 = qf(7,8);
plot3(xf2,yf2,zf2,'ro')

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