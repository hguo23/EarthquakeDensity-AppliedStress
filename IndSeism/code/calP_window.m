function [svec,svecneg,svecpos,P, P_low, P_high,dfe] = calP_window(stress,p,l,sl,varargin)
% This file calculates expectation of gaussian fit
% INPUT
% stress: applied stresses
% p: event density
% bin: edges of bins
% l, sl: smoothing parameters
% OUTPUT
% P: expectation
% P_low: minimum expectation at specific confidence interval
% P_high: maximum expectation at specific confidence interval
% dfe: degrees of freedom (used in tinv)

% Author: Huiyun Guo at UCSC (hguo23@ucsc.edu)
% Date: 2022 Jan 27

%%
conf = 90; % confidence interval

[stress,I] = sort(stress);
p = p(I);
svec(1) = nan;
svecneg(1) = nan;
svecpos(1) = nan;
P(1) = nan;
P_low(1) = nan;
P_high(1) = nan;
dfe(1) = nan;

count0 = 0;
% figure(100);clf;
if ~isempty(p)
    for i = 1:sl:length(p)-l
        count0 = count0+1;
        I = i:i+l;
        svec(count0) = 10.^(mean(log10(stress(I))));
        svecneg(count0) = min(stress(I));
        svecpos(count0) = max(stress(I));
        P(count0) = nan;
        P_low(count0) = nan;
        P_high(count0) = nan;
        dfe(count0) = nan;

        n = min([10 round(length(p(I))/3)]); % bin numbers for histogram
        [count,binedge] = hist(log10(p(I)),n);
        I = count>0; count = count(I); binedge = binedge(I);
        if sum(I)>7 % enough bins for fitting
            fo = fitoptions('gauss1',...
                            'Lower',[0,min(binedge),0],...
                            'Upper',[max(count)*1.5,max(binedge),Inf]...
                            );
            f = fit(binedge.',count.','gauss1',fo);
            fitp = f.a1*exp(-((binedge-f.b1)/f.c1).^2);
            c = sqrt(sum(((count-fitp)/max(count)).^2))/n;
            matrx = confint(f,conf/1e2);
            if c < 0.05 & (matrx(2,2)-matrx(1,2))<1.5
                P(count0) = 10^(f.b1);
                P_low(count0) = 10^(matrx(1,2));
                P_high(count0) = 10^(matrx(2,2));
                dfe(count0) = sum(I)-3;
            end   
            %% display
            if ~isempty(varargin)
                subplot(4,5,count0);hold on;
                if c<0.05
                    plot(f, binedge, count);
                else
                    scatter(binedge, count,'.');
                end
                xlabel('log_{10}\rho');
                ylabel('Number of Earthquakes');
                legend(num2str(c));
                title(['Bin No.',num2str(i)]);
%                 title([num2str(bin(i)),'< \sigma <=',num2str(bin(i+l))]);
                box on;hold off;
            end
            
        else
            if ~isempty(varargin)
                subplot(4,5,count0);hold on;
                scatter(binedge, count,'.');
                xlabel('log_{10}\rho');
                ylabel('Number of Earthquakes');
                title(['Bin No.',num2str(i)]);
%                 title([num2str(bin(i)),'< \sigma <=',num2str(bin(i+l))]);
                box on;hold off;
            end
        end           
    end
end
end
