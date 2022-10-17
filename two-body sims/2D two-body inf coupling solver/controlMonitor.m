function solns= controlMonitor(x0, x1, p0, muVec)
    solns= pertSolver(x0, x1, p0, muVec, 0);
    figure(1); hold on; figure(2); hold on; figure(3); hold on;
    for i= 1:length(solns)
        mu= solns{i}.mu;
        if mod(mu, 5) ~= 0
            continue
        end
        
        % plot u1
        figure(1);
        plot(solns{i}.t, 1/(2*mu+1)*((1+mu)*solns{i}.p(:,3)+mu*solns{i}.p(:,6)), ...
            'Color', hsv2rgb([(solns{i}.mu-muVec(1))/(muVec(end)-muVec(1))*0.9, 1, 1]));
        % plot u2
        figure(2);
        plot(solns{i}.t, 1/(2*mu+1)*(mu*solns{i}.p(:,3)+(1+mu)*solns{i}.p(:,6)), ...
            'Color', hsv2rgb([(solns{i}.mu-muVec(1))/(muVec(end)-muVec(1))*0.9, 1, 1]));
        % plot average of p3 and p6
        figure(3);
        plot(solns{i}.t, (solns{i}.p(:,3)+solns{i}.p(:,6))/2, ...
            'Color', hsv2rgb([(solns{i}.mu-muVec(1))/(muVec(end)-muVec(1))*0.9, 1, 1]));
    end
end