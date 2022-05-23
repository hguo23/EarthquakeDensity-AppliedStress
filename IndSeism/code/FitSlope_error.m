function [a,b,b_low,b_high,bvec] = FitSlope_error(x,y,ymin,ymax,dfe,conf,Nbstrp)
% This file fits y = 10^a * 10^(bx)
% inputs must be from gaussian-fit

% INPUT
% x: x
% y: y
% ymin: minimum y
% ymax: maximum y
% dfe: degrees of freedom (used in tinv)
% conf: confidence interval

% OUTPUT
% a,b: averaged parameters
% b_low, b_high: lower and upper boundary of b at conf
% bvec: vectors of b in bootsttrapping

% Author: Huiyun Guo at UCSC (hguo23@ucsc.edu)
% Date: 2022 Jan 27
%%
I = y <= 0 | isnan(y) | isinf(y) | x <= 0 | isnan(x) | isinf(x) | isnan(ymin) | isinf(ymin) | isnan(ymax) | isinf(ymax);
x = x(~I); y = y(~I); ymin = ymin(~I); ymax = ymax(~I); dfe = dfe(~I);
dy = (log10(ymax)-log10(ymin))/2; % amplitude range

G = [x'*0+1 log10(x')]; % least-square fit
N = size(G,1); 

if N>1     
    cstore = zeros(Nbstrp, 2); 
    for ib = 1:Nbstrp
        p = rand(1,N)*0.9+0.05; % random number between 0.05 to 0.95, based on previous conf in gaussian fitting
        t = tinv(p,dfe);
        f = log10(y)+t.*dy; % a random point between ymin and ymax based on the gaussian distribution
        cb = (G'*G)\(f*G)';
        cstore(ib,:) = cb;
    end
    a = mean(cstore(:,1));
    b = mean(cstore(:,2));
    [b_low, b_high] = ConfidenceInter(conf, cstore(:,2));
    % ConfidenceInter: a confidence interval function
    bvec = cstore(:,2);
else
    a = nan;
    b = nan;
    b_low = nan;
    b_high = nan;
    bvec = nan;
end
end