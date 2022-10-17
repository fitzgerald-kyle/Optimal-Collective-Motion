function circleOptimalityNonCoplanar(s)
% Inputs:
%   s : scalar coupling strength

    close all
    
    numPts= 2^6+1;
    theta= linspace(-pi/4, 3*pi/4, numPts); % for sampling points on a circle
                                   % in the mu2(0)-mu8(0) plane (going from -pi/4
                                   % to 3pi/4 lets us exploit origin symmetry
                                   % in (mu2,mu8) and do half the calculations)
    
    for mu3= 0:pi/16:pi/2
        mu2Vec= [];    
        mu3Vec= [];
        mu8Vec= [];
        for t= theta(1):pi/(length(theta)-1):theta(end)
            if mod(t,pi/4) == 0 % display at every 25% of calculations (for sanity)
                disp(['done with ', num2str((t+pi/4)*4/pi),'/4', ' for mu3=', num2str(mu3)])
            end
            
            mu2= 3*pi*cos(t); mu8= 3*pi*sin(t);
            
            output_IVP= main_IVP_solver([0 mu2 mu3 0 0 0 0 mu8 0 0 0 0],s,false);
            
            if isempty(output_IVP.tconj)
                % display a warning if the solver didn't find any conj pts; you
                % probably didn't integrate far enough along in time (increase
                % tf in main_IVP_solver)
            warning(['Didn''t capture conj pt for (mu2,mu8,mu3)=(' \num2str(mu2) ...
                ',' num2str(mu8) ',' num2str(mu3) ')']);
            else
                % multiply the initial condition by the conj pt time in order
                % to determine the initial condition at which the first conj
                % pt occurs
                mu2Vec(end+1)= mu2*output_IVP.tconj(1);
                mu8Vec(end+1)= mu8*output_IVP.tconj(1);
                mu3Vec(end+1)= mu3*output_IVP.tconj(1);
            end
        end

        figure(1);
        plot3(mu2Vec, mu8Vec, mu3Vec, 'Marker', '.', 'MarkerSize', 15, ...
            'Color', hsv2rgb([0.499*(mu3/mu3MAX+1),1,1]), 'LineStyle','none');
        hold on; grid on
        % exploit planar origin symmetry
        plot3(-mu2Vec, -mu8Vec, mu3Vec, 'Marker', '.', 'MarkerSize', 15, ...
            'Color', hsv2rgb([0.499*(mu3/mu3MAX+1),1,1]),'LineStyle','none');
        title(['Location of First Conjugate Point wrt \mu_2(0), \mu_8(0), \mu_3(0) for \lambda = ', ...
            num2str(s)]);
        xlabel('\mu_2(0)'); ylabel('\mu_8(0)'); zlabel('\mu_3(0)');
        view(3)
    end
end