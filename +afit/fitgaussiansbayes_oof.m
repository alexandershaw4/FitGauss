function [P,F] = fitgaussiansbayes(t,x,n,w)
global dx dt p0 n0
%
%
% AS2019

dx = x;
dt = t;
n0 = n;

if n == 1
    [PKS,LOCS]=max(x);
    W = length(t)/8;
    LOCS = t(LOCS);
    try 
        W = w;
    end
else

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
end


% Initial parameters and points in state space
for i    = 1:n
    f(i) = (LOCS(i));%-1;
    a(i) = PKS(i)/n;
    w(i) = W(i)/(n*n);
end

% Noise shape (one-over-f squared [oof]) parameters:
% a * (1./f.^b);
noi.a = 1e-2;
noi.b = 2;

% Parameter set
p.f = f;
p.a = a;
p.w = w;
p.n = noi;
X   = spm_vec(p);
p0  = p;

% Set it up as per a dynamic model
M    = [];
M.IS = @obj;
M.pE = X;
M.pC = diag(~~X)*4;

xU    = [];
xU.u  = 0;

xY.y  = x;
xY.Hz = t;
M.x   = [];
M.m   = 1;
M.l   = 1;

% Call the external Bayesian EM routine
[Qp,Cp,Eh,F] = spm_nlsi_GN(M,xU,xY);

%V = (~~diag(M.pC));
%V(1:3) = 1/16; % f
%V(4:6) = 1/64; % a
%V(7:9) = 1/8 ;   % w
%V(10:11) = 1/32;  %n
%[Qp,F,CP,PP,H] = AO(@obj,X,V,x,128,[],[],-400,[],2,[],'fe',1,1,1);

% return parameters in state space
P = spm_unvec(Qp,p);


function e = obj(X,z0,z1,varargin)
global dx dt p0 n0

% n bump
XP = spm_unvec(X,p0);
e  = afit.makef(dt(:)',XP.f,XP.a,XP.w);

noise = XP.n.a * (1./dt.^XP.n.b);

e = e + noise;




% % plot?
%plot(1:length(dx),dx); hold on
%plot(1:length(dx),Y); hold off;  drawnow;

% error
%e = sum( dx(:) - Y(:) ).^2;
    
    