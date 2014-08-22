function nok = check_monotonic_columns(M)
% Check for a given transition matrix M, whether the column values 
% above and below the diagonal decrease monotonically.
% Return true if there are still violoation on the columns.
%
% Input : M transition matrix (NxP) including default column
% Output: nok: true if matrix does not satisfy monoticity column conditions
% Author: Fabian Ojeda

% omit default column
Mt = M(:,1:end-1);

% check lower triangle - for non-mononotically decreasing elements
nok1 = any(any(tril(diff(tril(Mt),1,1),-1)  > 0));
% check upper triangle - for non-mononotically decreasing elements
nok2 = any(any(tril(diff(rot90(triu(Mt),2),1,1),-1)  > 0));
nok = nok1 | nok2

end