function [P,F] = fitgaussiansbayespriors(t,x,priors)
global dx dt p0 n0
%
%
% AS2019

dx = x;
dt = t;
n0 = length(priors);

% Initial parameters and points in state space
n        = length(priors);
for i    = 1:n
    f(i) = priors(i)-1;
    a(i) = mean(x)/n;
    w(i) = 4      /(n);
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
    
    