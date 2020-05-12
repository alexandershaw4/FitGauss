function [P,F] = fitgaussians(t,x,n,optfun)
global dx dt p0 n0

dx = x;
dt = t;
n0 = n;

[PKS,LOCS,W] = findpeaks(x,t,'NPeaks',n);

% initialparameters and point
for i    = 1:n
    %f(i) = t(LOCS(i));
    f(i) = t(LOCS(i));
    a(i) = PKS(i)/n;
    w(i) = W(i)/n;
end

% reversable parameter set
p.f = f;
p.a = a;
p.w = w;
X   = spm_vec(p);
p0  = p;

%[X1,F] = fminsearch(@obj,X);
if nargin < 4 || isempty(optfun)
    optfun = @fminsearch;
end

if ~strcmp(char(optfun),'AO') && ~strcmp(char(optfun),'AObayes')
    [X1,F] = optfun(@obj,X);
else
    V = ones(size(X))/128;
    [X1,F] = optfun(@obj,X,V,[],[],[],[],1e-60);
end



% return parameters in state space
P = spm_unvec(X1,p);


function e = obj(X)
global dx dt p0 n0

% n bump
XP = spm_unvec(X,p0);

%for i = 1:n0
    % compute
    %Y(i,:) = afit.makef(dt(:)',XP.f(i),XP.a(i),XP.w(i));
    
%end

Y = afit.makef(dt(:)',XP.f,XP.a,XP.w);

%Y = sum(Y,1);


% % plot?
plot(1:length(dx),dx); hold on
plot(1:length(dx),Y); hold off;  drawnow;

% error
e = sum(sum( dx(:) - Y(:)' )).^2;
    
    