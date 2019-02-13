function [P,F] = fitgaussiansbayes(t,x,n)
global dx dt p0 n0
%
%
% AS2019

dx = x;
dt = t;
n0 = n;

% We'll try getting initial (prior) points using 'findpeaks'
[PKS,LOCS,W] = findpeaks(x,t,'NPeaks',n);

% If find peaks hasn't found good initial points, just even space the
% number of guestimated bumps along t
if length(LOCS) ~= n
    P0 = PKS;
    L0 = LOCS;
    W0 = W;
    
    ddt  = floor( length(t)/(n+1) );
    LOCS = t(ddt:ddt:ddt*n);
    PKS  = repmat( mean(x)/n , [1 n]);
    W    = repmat( 1         , [1 n]);
    
    % use whatever initial points find-peaks managed to find...
    for i = 1:length(L0)
        LOCS(i) = L0(i);
        PKS(i)  = W0(i);
        W(i)    = W0(i);
    end
    
end


% Initial parameters and points in state space
for i    = 1:n
    f(i) = (LOCS(i))-1;
    a(i) = PKS(i)/n;
    w(i) = W(i)/(n*n);
end

% Parameter set
p.f = f;
p.a = a;
p.w = w;
X   = spm_vec(p);
p0  = p;

% Set it up as per a dynamic model
M    = [];
M.IS  = @obj;
M.pE = X;
M.pC = diag(~~X);

xU    = [];
xU.u  = 0;

xY.y  = x;
xY.Hz = t;
M.x   = [];
M.m   = 1;
M.l   = 1;

% Call the external Bayesian EM routine
[Qp,Cp,Eh,F] = spm_nlsi_GN(M,xU,xY);

% return parameters in state space
P = spm_unvec(Qp,p);


function e = obj(X,z0,z1,varargin)
global dx dt p0 n0

% n bump
XP = spm_unvec(X,p0);
e  = afit.makef(dt(:)',XP.f,XP.a,XP.w);




% % plot?
%plot(1:length(dx),dx); hold on
%plot(1:length(dx),Y); hold off;  drawnow;

% error
%e = sum( dx(:) - Y(:) ).^2;
    
    