function P = SelectOperator(combination,ngaus)
% construction of the gappy selection algorithm 
if nargin == 0 
    combination = [ 3 4]'; 
    ngaus = 4;
end

P = zeros(length(combination),ngaus); 

iind = 1:length(combination);  iind = iind';
subindices = [iind combination];
combination = double(combination);
IndicesLinear = sub2ind(size(P),iind,combination) ;

P(IndicesLinear) = 1 ;

 