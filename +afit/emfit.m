function [Qp,Cp,Eh,F] = emfit(M,U,Y,n)

% integration function
y = Y.y;

% estimate spectral peaks
[P,F] = afit.fitgaussiansbayes(Y.Hz,y{1},n); close

% mock up a spectrum to fit, only featuring narrowband peaks at determined
X = afit.makef(M.Hz,P.f,P.a,P.w/4);

% fit this narrowband spectrum first
Y0      = Y;
Y0.y{1} = X;

% vB EM routine: fit the peaks
[Qp,Cp,Eh,F] = spm_nlsi_GN(M,U,Y0);

% or sampling routine
%M.FS         = 'spm_fs_csd';
%[Qp,qC,qh,F] = spm_nlsi_LS(M,U,Y0);

% now refit the 'full' data using the posteriors of the narrowband fit
M.pE         = Qp;
[Qp,Cp,Eh,F] = spm_nlsi_GN(M,U,Y);
