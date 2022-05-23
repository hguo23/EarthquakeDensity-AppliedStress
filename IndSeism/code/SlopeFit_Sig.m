clearvars p p_low p_high dfe
I = ~isnan(sig) & sig>0;
sig = sig(I);
eqd = eqd(I);

N = max([round(length(sig)/l) 3e1]);
dN = N;
[binp,binpneg,binppos,p, p_low, p_high,dfe] = calP_window(sig,eqd,N,dN);
% [binp,binpneg,binppos,p, p_low, p_high,dfe] = calP_window(sig,eqd,N,dN,'display');

if length(unique(binp(p>0)))>2
    [~,~,~,~,bvec] = FitSlope_error(binp,p,p_low,p_high,dfe,conf,Nbstrp);
else
    bvec = nan(Nbstrp,1); 
end

