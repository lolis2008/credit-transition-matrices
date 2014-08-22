function nok = check_monotonic_rows(M)
% For a given transition matrix M, this function checks whether there are 
% still violoation on the rows,
% that is if values on each side (left and right) of the diagonal entry do
% not decrease monotonically
%
% Input : M transition matrix (NxP) including default column
% Output: nok: true if matrix does not satisfy monoticity row condition
% Author: Fabian Ojeda

% omit default column
Mt = M(:,1:end-1);

% check upper triangle - for non-mononotically decreasing elements
nok1 = any(any(triu(diff(triu(Mt),1,2),1)  > 0));
% check lower triangle - for non-mononotically decreasing elements
nok2 = any(any(triu(diff(rot90(tril(Mt),2),1,2),1)  > 0));
nok = nok1 | nok2;

end