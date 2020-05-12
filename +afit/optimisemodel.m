function [Qp,Cp,Eh,F,stuff] = optimisemodel(DCM)
global M U P0

pE = DCM.M.pE;
M  = DCM.M;
U  = DCM.xU;
P0 = pE;

[y,w,s,g,t,pst] = feval(M.IS,pE,M,U);

% % optimsation
DCM2      = DCM;
DCM2.M.IS = @evaluate;
DCM2.xY.y = {targets(pst,DCM.xY)};

% check 'integration'
yint = feval(DCM2.M.IS,pE,M,U);

% remove zero-contirbuting states
J = exp(DCM.M.pE.J (1:8) );
DCM2.xY.y{1} = DCM2.xY.y{1} .* repmat(~~J,[size(yint,1),1]);

% fit
%DCM2.M.nograph = 0;
%[Qp,Cp,Eh,F]   = spm_nlsi_GN(DCM2.M,DCM2.xU,DCM2.xY);

[Qp,F,Cp] = AO_DCM(pE,DCM2,36,0,'fe');
Eh = [];

try
    stuff.targets    = {targets(pst,DCM.xY)};
    stuff.prediction = evaluate(Qp,M,U);
end

end

function e = objective(pE)
global M U P0 Y

pE = spm_unvec(pE,P0);
[p,pst] = evaluate(pE,M,U);
t       = targets(pst,Y);

% normalise (make scale free error)
e       = sum( sum( t - p ) ).^2;

% plots?
plot(pst,t,':',pst,p); drawnow;

end

function t = targets(pst,Y)

dt = pst(2) - pst(1);
dt = dt/1000;

fq    = [10 50 50 20 20 15 8 8];
[P,F] = afit.fitgaussiansbayespriors(Y.Hz,Y.y{1},fq); close;

for i = 1:length(fq)
    t(i,:) = afit.makef(Y.Hz(:)',P.f(i),P.a(i),P.w(i));
end

t=abs(t)';

end

function [PopTS,pst] = evaluate(pE,M0,U)
global M
[y,w,s,g,t,pst,l] = feval(M.IS,pE,M,U);

lay = squeeze( l{1}.weighted );

lay0 = zeros(8,size(lay,2));
Ji = find(exp(pE.J));
lay0(Ji,:) = lay;
pf = lay0;

% PopTS = squeeze(s{1}(1,:,1,:));
% 
% J = exp(pE.J (1:8) );
% PopTS = PopTS .* repmat(J', [1 size(PopTS,2)] );
% 
% %ssr = findthenearest(pst,50);
% 
% for i = 1:size(PopTS,1)
%     %pf(i,:) = Afft(PopTS(i,ssr:end),1/M0.dt,4:80);
%     pf(i,:) = Afft(PopTS(i,:),1/M0.dt,4:80);
% 
% end

PopTS = {pf'};
pst = w;
end


