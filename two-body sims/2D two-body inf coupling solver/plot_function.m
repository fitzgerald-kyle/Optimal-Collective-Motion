function plot_function(output,xf)

% This function creates a plot of the solution contained in the variable
% output
% Function Inputs:
%   output : structure containing the current solution
%   xf : desired boundary condition for x(tf)
% Function Outputs:
%   none

% Plot solution
plot(output.x(:,1),output.x(:,2),'b-')
hold on
plot(output.x(:,4),output.x(:,5),'c-')

% Plot initial condition for x at t=0
plot(output.x(1,1),output.x(1,2),'go')
plot(output.x(1,4),output.x(1,5),'go')

% Plot desired boundary condition for x at t=tf
plot(xf(1),xf(2),'ro')
plot(xf(4),xf(5),'ro')

hold off

% Set axis limits
axis([-0.5 1 -1 1])

% Axis labels
xlabel('x_1')
ylabel('x_2')

% Set axis aspect ratio
daspect([1 1 1])

% Set font size
set(gca,'color','w','FontSize',20)

% Set plot frame color to white
set(gcf,'color','w')

% Draw plot now
drawnow

end