function [u2Vec, u5Vec] = circleOptimalityCoplanar(numPts, s)
% This function creates a graph of ordered pairs (u2(0), u5(0)) at which
% the first conjugate point occurs for coplanar circular trajectories with
% u2(0) and u5(0) as initial conditions. This is achieved by exploiting
% the scalability of the system.

% Inputs:
%   s : scalar coupling strength

    %close all
    
    theta= linspace(-pi/4 ,pi-pi/4, numPts); % for sampling points on a circle
                                             % in the mu2(0)-mu8(0) plane
    
    % vectors that keep track of the first conj pt locations
    mu2Vec= [];    
    mu8Vec= [];
    r= 8*pi;    
    for i= 1:length(theta)
        if numPts > 4 && mod(theta(i),pi/4) == 0 % display at every 25% of calculations (for sanity)
            disp(['done with ', num2str((theta(i)+pi/4)*4/pi),'/4'])
        end
        
        u2= r*cos(theta(i)); u5= r*sin(theta(i));
        
        if u2==0 || u5==0
            continue
        end
        
        temp= [1+s -s; -s 1+s]*[u2; u5];
        mu2= temp(1); mu8= temp(2);
        
        output_IVP= main_IVP_solver([0 mu2 0 0 0 0 0 mu8 0 0 0 0],2,s,false);
        
        if isempty(output_IVP.tconj)
            % display a warning if the solver didn't find any conj pts; you
            % probably didn't integrate far enough along in time (increase
            % tf in main_IVP_solver)
            warning(['Didn''t capture conj pt for (mu2,mu8)=(' num2str(mu2) ...
                ',' num2str(mu8) ')']);
        else
            % multiply the initial condition by the conj pt time in order
            % to determine the initial condition at which the first conj
            % pt occurs
            mu2Vec(end+1)= mu2*output_IVP.tconj(1);
            mu8Vec(end+1)= mu8*output_IVP.tconj(1);
        end
    end
    
    u2Vec = ((1+s)*mu2Vec+s*mu8Vec)/(2*s+1);
    u5Vec = (s*mu2Vec+(1+s)*mu8Vec)/(2*s+1);
    
    if numPts == 1
        return;
    end
    
    figure; hold on
    plot(u2Vec, u5Vec, '.r', 'MarkerSize', 15);
    hold on; grid on 
    % exploit the system's origin symmetry
    plot(-u2Vec, -u5Vec, '.r', 'MarkerSize', 15);
    title(['Location of First Conjugate Point wrt u_2(0) and u_5(0) for \lambda = ', ...
        num2str(s)]);
    xlabel('u_2(0)'); ylabel('u_5(0)');
    daspect([1 1 1]);
end