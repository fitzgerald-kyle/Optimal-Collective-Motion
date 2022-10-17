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

% Plot initial condition for x at t=0
plot(output.x(1,1),output.x(1,2),'go')

% Plot desired boundary condition for x at t=tf
plot(xf(1),xf(2),'ro')

hold off

% Set axis limits
axis([min(output.x(:,1)) max(output.x(:,1)) ...
    min(output.x(:,2)) max(output.x(:,2))])
%axis([-0.5 1 -1 1])

% Axis labels
xlabel('x_1')
ylabel('x_2')

% Set axis aspect ratio
%daspect([1 1 1])

% Set font size
set(gca,'color','w','FontSize',20)

% Set plot frame color to white
set(gcf,'color','w')

% Draw plot now
drawnow

end