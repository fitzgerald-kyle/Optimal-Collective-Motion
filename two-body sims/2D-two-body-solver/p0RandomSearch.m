function p0RandomSearch(x0, x1, mu)
    t= 0; tPrev= 0;
    mu= mu - 0.05;
    tic;
    while t < 25000 && mu <= 2
        mu= mu + 0.05;
        p0= zeros(1,6);
        p0Range= [500,500,6,500,500,6];
        while t-tPrev < 600
            for i= 1:6
                p0(i)= p0Range(i)*rand - 0.5*p0Range(i);
            end
            fprintf('Trying (%.2f %.2f %.2f %.2f %.2f %.2f)... \n', ...
                p0(1), p0(2), p0(3), p0(4), p0(5), p0(6));
            pertSolver(x0, x1, p0, mu:-.01:0, 0);
            t= toc;
        end
        tPrev= t;
    end
end