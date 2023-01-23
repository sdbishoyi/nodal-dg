function [dlogAbsF_dlogt] = ComputeTailDecayRate(T,F)
%
logT   =log10(T);
logAbsF=log10(abs(F));
%
DlogAbsF = logAbsF(5:end) - logAbsF(1:end-4);
DlogT    = logT(5:end)    - logT(1:end-4);

% analytically the log time spacings should be this when T is equally space
%DtTest = transpose((1+(1:length(T)-1))./ (1:length(T)-1));
%logDtTest=log10(DtTest);
%max(abs(DlogT-logDtTest));
dlogAbsF_dlogt = DlogAbsF./DlogT;