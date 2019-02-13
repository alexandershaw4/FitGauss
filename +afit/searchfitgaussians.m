function [P,F,N,Y,AllF] = searchfitgaussians(t,x,routine,npeaks)
% A wrapper on fitgaussians that tries to optimise the number of gaussians
% that are required to fit the data x at locations t by performing a greedy
% search. Fitgaussians.m uses fminsearch to fit a multi-gaussian model of
% order n. Fitgaussiansbayes uses a Bayesian EM routine to minimse a 
% complexity-adjusted error term (free-energy).
%
% This routine recursively calls fitgaussians(bayes) with different n's
% and evaluates them. It can be considered a model compairson routine of
% different dimensionality gaussian mixture models.
%
%
% Generate some example data and call the fit / search:
%
% w = 1:100 ; % time / freq / unit vector
% S = afit.makef(w,[10 20 50 70],[10 8 5 3],[1 2 5 3]); % some data
% 
% plot the input:
% figure, plot(w,S)
%
% If I know I want to search for 4 Gaussian components:
% [P,f] = afit.fitgaussiansbayes(w,S,4);
%
% returns:
% P = 
%   f: [11 21 51 69]
%   a: [8.7244 7.9161 4.8467 2.5846]
%   w: [1.4209 1.6610 5.3906 2.8153]%
%
% If I don't know how many components there are:
% [P,F,N,Y,AllF] = afit.searchfitgaussians(t,x); % fit it
%
% Or specifya first guess number (e.g. 5):
% [P,F,N,Y,AllF] = afit.searchfitgaussians(t,x,[],5); % fit it
%
% AS2019


if nargin < 4 || isempty(npeaks)
    PKS  = findpeaks(x,t,'SortStr','descend');
    MAXC = length(PKS);
    nc   = fliplr(1:MAXC);
else
    MAXC = npeaks;
    nc   = fliplr(1:npeaks);
end


if nargin < 3 || isempty(routine)
    routine = 'bayes';
end


for i = 1:MAXC
    
    fprintf('Searching iteration %d / %d\n',i,MAXC);
    
    switch routine
        case 'fminsearch'
            [P(i),F(i)] = afit.fitgaussians(t,x, nc(i) );
        case 'bayes'
            [P(i),F(i)] = afit.fitgaussiansbayes(t,x, nc(i) );
    end
    close
end

AllF    = F;
[F,ID]  = max(F); % winning number of components
N       = nc(ID);
P       = P(ID);

fprintf('Model %d won with %d components... \n',ID, N);

for i = 1:N
    Y0(i,:) = afit.makef(t,P.f(i),P.a(i),P.w(i));
end

Y.Model      = sum(Y0,1);
Y.Components = Y0;
Y.Input      = x;
Y.Units      = t;

