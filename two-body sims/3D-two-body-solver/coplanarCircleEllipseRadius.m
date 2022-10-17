close all

sMax = 10000;
numPts = 100;

s = logspace(-log10(sMax), log10(sMax), numPts)';
r = zeros(numPts, 1);
for i = 1:numPts
    [u2Vec, u5Vec] = circleOptimalityCoplanar(1, s(i));
    if ~isempty(u2Vec)
        r(i) = sqrt(u2Vec^2+u5Vec^2);
    end
    if mod(i, round(numPts/10)) == 0
        disp(['Done with s=' num2str(s(i))]);
    end
end

plot(s, r);
set(gca,'XScale','log')

%{
fo = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [0; 4.5], ...
'Lower', [0; -inf]);
ft = fittype('1/(s+b) + d', 'independent', 's', 'dependent', 'r', 'options', fo);
[curveFit,~] = fit(s, r, ft);
%}

%
fo = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [1; 0; 4.5], ...
    'Lower', [-inf; 0; -inf]);
ft = fittype('a/(s+b) + d', 'independent', 's', 'dependent', 'r', 'options', fo);
[curveFit,~] = fit(s, r, ft);
%}

%{
fo = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [1; 0; .1; 4.5], ...
    'Lower', [-inf; 0; -inf; -inf]);
ft = fittype('a/(s+b)^c + d', 'independent', 's', 'dependent', 'r', 'options', fo);
[curveFit,~] = fit(s, r, ft);
%}
hold on
plot(s, feval(curveFit, s));
xlabel('s'); ylabel('Semi-minor axis');