% Compute steady state probabilities given transition probabilities
% PaGq
% June 2018
% Input: P stochastic matrix with columns adding to 1
% Output: pie vector of ss probabilities

function pie = ss_prob(P)

mm = size(P,1) ;
AA = [eye(mm) - P;ones(1,mm)] ;
TT = (AA'*AA)\AA' ;

pie = TT(:,end) ;