function plot_function(output,xf,mu)

% This function creates a plot of the solution contained in the variable
% output
% Function Inputs:
%   output : structure containing the current solution
%   xf : desired boundary condition for x(tf)
% Function Outputs:
%   none

figure('visible', 'off');

% Plot solution
plot(output.x(:,1),output.x(:,2),'b-','LineWidth',5)
hold on
plot(output.x(:,4),output.x(:,5),'c-','LineWidth',5)

% Plot initial condition for x at t=0
plot(output.x(1,1),output.x(1,2),'go')
plot(output.x(1,4),output.x(1,5),'go')

% Plot desired boundary condition for x at t=tf
plot(xf(1),xf(2),'ro')
plot(xf(4),xf(5),'ro')

hold off

% Set axis limits
axis([min( min(output.x(:,1)), min(output.x(:,4)) ), ...
    max( max(output.x(:,1)), max(output.x(:,4)) ), ...
    min( min(output.x(:,2)), min(output.x(:,5)) ), ...
    max( max(output.x(:,2)), max(output.x(:,5)) )])
%axis([-0.5 1 -1 1])


% Axis labels
xlabel('x_1')
ylabel('x_2')

% Set axis aspect ratio
daspect([1 1 1])

% Set font size
%set(gca,'color','w','FontSize',20)
set(gca, 'Color', 'none');
set(gca,'Visible','off');

% Set plot frame color to white
%set(gcf,'color','w')

[folderName, fileName] = createFolderAndFileNames(output.x(1,:), xf, output.p(1,:), ...
    mu);
fileName = [fileName(1:end-4), sprintf('(seg=%.2f).png',mu)];
cdir= cd('altmany-export_fig-b1a7288');
export_fig(['../',folderName,'_IMAGES/',fileName], '-transparent');
cd(cdir);

% Draw plot now
drawnow

end