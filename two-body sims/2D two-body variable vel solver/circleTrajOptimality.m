function circleTrajOptimality(x0)
    muVec= [1 10 100];
    for i= 1:3
        suboptp03= [];
        suboptp06= [];
        optp03= [];
        optp06= [];
        for theta= 0:.1:6.2
            fprintf('solving theta=%.1f\n', theta)
            for n= 0:8
                L= 2^n;
                output_IVP= solve_IVP(x0,[0 0 cos(theta) 0 0 sin(theta)], ...
                    1,parameters(1000),muVec(i),2^n);
                if ~isempty(output_IVP.tconj)
                    %optimal= 0;
                    suboptp03(end+1)= L*cos(theta);
                    suboptp06(end+1)= L*sin(theta);
                else
                    optp03(end+1)= L*cos(theta);
                    optp06(end+1)= L*sin(theta);
                end
            end
        end
        figure(i);
        plot(suboptp03, suboptp06, '.r');
        hold on
        plot(optp03, optp06, '.b');
        title(['Circle Optimality for mu=' num2str(muVec(i))]);
        xlabel('p_3(0)'); ylabel('p_6(0)');
        drawnow
    end
end