function [LowerBound, HigherBound] = ConfidenceInter(ConfidenceInterval,v)
LowerBound = zeros(length(ConfidenceInterval), 1);
HigherBound = zeros(length(ConfidenceInterval), 1);
for q = 1: length(ConfidenceInterval)
    confidence_interval = ConfidenceInterval(q);
    q1 = prctile(v, (100-confidence_interval)/2);
    q2 = prctile(v, 100-(100-confidence_interval)/2);
    LowerBound(q) = q1;
    HigherBound(q) = q2;
end
