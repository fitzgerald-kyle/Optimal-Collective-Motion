function output= pSpaceSearch(x0, x1, p0OG, mu0, numPts, dp1, dp2, dp3)
    close all

    numTraj= 10;    

    params= muSearchParameters(200);
    output.optimalp0= [];
    output.nonoptimalp0= [];
    output.optimalMu= [];
    output.nonoptimalMu= [];
    optIdx= 1; nonoptIdx= 1;
    for j= 1:numTraj
        p0 = p0OG;
        mu= mu0;
        for i= 1:ceil(numPts/numTraj)   
            p0= p0 + [2*dp1*rand-dp1, 2*dp2*rand-dp2, 2*dp3*rand-dp3, ...
                2*dp1*rand-dp1, 2*dp2*rand-dp2, 2*dp3*rand-dp3];
            
            output_IVP= solve_IVP(x0,p0,1,params,mu0);
            err= output_IVP.x(end,:)-x1;
            lastx= output_IVP.x(end,:);
            iter= 1; dMu= 1e-12;
            while norm(err) > params.tol && iter <= params.nmax
                if dMu > 1
                    mu= mu+1;
                else
                    mu= mu+dMu;
                end
                output_IVP= solve_IVP(x0,p0,1,params,mu);
                err= output_IVP.x(end,:)-x1;
                muJacobian= (output_IVP.x(end,:)-lastx) / dMu;
                dMu= -muJacobian' \ err';
                
                dMu= line_search(x0,p0,x1,mu,dMu,norm(err),params);
                
                lastx= output_IVP.x(end,:);
                %fprintf('iteration %i: error %.6f\n', iter, norm(err));
                iter= iter+1;
            end
            
            if isnan(dMu)
                continue
            end
                        
            if isempty(output_IVP.tconj)
                output.optimalp0(optIdx,:)= p0;
                output.optimalMu(optIdx)= mu;
                optIdx= optIdx+1;
            else
                output.nonoptimalp0(nonoptIdx,:)= p0;
                output.nonoptimalMu(nonoptIdx)= mu;
                nonoptIdx= nonoptIdx+1;
            end
        end
    end    
    
    figure
    %axis equal
    if ~isempty(output.optimalp0)
        plot3(output.optimalp0(:,1), output.optimalp0(:,2), output.optimalp0(:,3), ...
            '.b', 'MarkerSize', 15);
        hold on
    end
    if ~isempty(output.nonoptimalp0)
        plot3(output.nonoptimalp0(:,1), output.nonoptimalp0(:,2), output.nonoptimalp0(:,3), ...
            '.r', 'MarkerSize', 15);
    end
    xlabel('p1(0)'); ylabel('p2(0)'); zlabel('p3(0)');
    grid on
end


function dMu = line_search(x0,p0,x1,mu,dMu,err,params)

% This function uses a line search to ensure that the error decreases when
% updating p0
% Function Inputs:
%   x0 : initial condition for x(0)
%   p0 : guess for initial condition p(0) 
%   xf : desired boundary condition for x(tf)
%   dp0 : search direction for p0
%   tf : end of time interval
%   err : err at current value of p(0)
%   params : structure containing parameters
% Function Outputs:
%   step :  step size that ensures error decreases in the direction dp0

% Initialize step size and new error
step = 1/params.step_ratio;
newerr = step*err;
prevErr= err;
    
% Line search loop
while newerr > (1-params.decrease_param)*err || newerr > prevErr
    
    % Compute step size and p0 after step
    step = params.step_ratio*step;
    dMu= step*dMu;
    newMu = mu+dMu;
    
    % Compute error after step
    output_IVP = solve_IVP(x0,p0,1,params,newMu);
    
    % Store error between current and desired value of x(tf)
    neweta = [output_IVP.x(end,1)-x1(1),... 
              output_IVP.x(end,2)-x1(2),...
              wrapToPi(output_IVP.x(end,3))-x1(3),...
              output_IVP.x(end,4)-x1(4),... 
              output_IVP.x(end,5)-x1(5),...
              wrapToPi(output_IVP.x(end,6))-x1(6)];
    prevErr= newerr;
    newerr = norm(neweta);
    
    % Check if step is less than minstep
    if norm(step) < params.minstep
        break
    end 
end
    
end
    
% change solveBVP to have dMu instead of dp0;
% plot a 3D graph of p01-p02-p03 with different colors for different
% numbers of non-optimal solutions (multiple for different p04-p05-p06);
% spreadsheet of solutions?;
% round p0 to hundredth then calculate?