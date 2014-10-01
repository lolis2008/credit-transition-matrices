function M = MMWeiSmoothing(M,per,AR,NR,normdiag)
% Pre-processing for Migration Matrices as described in
% J. Wei (2000) "A Multi-Factor, Markov Chain Model for Credit Migrations and Credit
% Spreads", pg. 12.
%
%% Inputs
% M        : migration matrix of N x P
% per      : 1 if values are given as percentages [0-100] range 
%            (otherwise are probabilities in [0-1] range (default 0)
% AR       : 1 if last row contains the absorbing default row (default 0)
% NR       : 1 if last column contains Not Rated values (default 0)
% normdiag : include diagonal elements when normalizing rows (default 0)
%% Output
% M : migration matrix of N x P
%
%% Written by
% Fabian Ojeda (KBC Risk Modelling Department)
%

%% input parameters defaults
if nargin < 2
    per = 0;
end
   
if nargin < 3
    AR = 0;
end

if nargin < 4 
    NR = 0;
end

if nargin < 5
    normdiag = 0;
end

if per
    per = 100;
else
    per = 1;
end

if AR
    AR = M(end,:);
    M(end,:) = [];
else
    AR = [];
end


%% default column assumed to be last column
D = size(M,2);


%% STEP 1: redistribution of not rated (NR) to other ratings except default (D)
if NR
    % default column second from right
    D = size(M,2)-1; 
    % weighting based on NR column
    w = (per+M(:,end))./per;
    % remove NR column
    M(:,end) = [];
    % apply weighting (not on default)
    T = bsxfun(@times,M(:,1:(D-1)),w);
    % insert back ignoring default column
    M(:,1:(D-1)) = T;
end


%% STEP 2: row values should decrease on each side of diagonal
Mt = monotonic_rows(M); % operates first on upper triangle
Mt = monotonic_rows(Mt,1); % operates on lower triangle


%% STEP 3: column values should decrease on each side of diagonal
% do this step until matrix converges 
max_iter = 100;
iter = 0;
while check_columns(Mt) && iter < max_iter
    Mt = monotonic_columns(Mt,per);% operates first on lower triangle
    Mt = monotonic_columns(Mt,per,1); % operates on upper triangle
    iter = iter + 1;
end

%% rescales rows to sum up to (1 or 100)
ndf = 1:D-1; % index excluding default column
n = size(Mt,1);
diagM = zeros(n,1);
if ~normdiag
    diagM = diag(Mt(:,ndf)); % diagonal elements
end
row_per = per - diagM - Mt(:,D); % row percentage excluding default (and diagonal)
sum_row = sum(Mt(:,ndf),2) - diagM; % row sum excluding default (and diagonal)

% avoid zero off diagonal rows 
inz = sum_row > 0;
Mt(inz,ndf) = bsxfun(@times,Mt(inz,ndf),row_per(inz)./sum_row(inz));

% diagonal index
idx = 1:n+1:n^2; 
if ~normdiag
    Mt(idx) = diagM; % restore diagonal
end
M = [Mt;AR]; % reinsert absorbing state if given

end
    
function M = monotonic_rows(M,lower)
% 
% Probability should decline monotonically on each side of the diagonal entry. 
% Whenever there is a violation, the entry is set equal to the previous
% rating’s entry and the difference is equally distributed among the entries 
% between the diagonal entry
% and the entry in question
% 
% M     : migration matrix N x P 
% lower : 1 process the upper triangle (default 0) 

if nargin < 2
    lower = 0;
end

% omit default column
D = M(:,end);
M(:,end) = [];

if lower
    % rotates the lower triangle such matrix can be process left-right
    % top-bottom
    M= rot90(M,2); 
end

% loop upper part matrix
for i=1:size(M,1)
    for j = i:size(M,2)-1
        % compare values
        trans = M(i,j+1) - M(i,j);
        if trans > 0 %M(i,j+1) > M(i,j)
            % replace with previou's rating
            M(i,j+1) = M(i,j);
            % distribute difference 
            % among the entries between the diagonal entry and the entry in question
            idx = i+1:j; 
            if numel(idx) > 0 
                M(i,idx) = M(i,idx) + trans/numel(idx);
            end            
        end
    end
end

if lower
    % revert rotation 
    M = rot90(M,2);
end

% reinsert default column
M = [M D];
end


function M = monotonic_columns(M,per,upper)
% within each column, the entries on each side of the diagonal entry
% should also monotonically decline. To minimize excessive arbitrary adjustments, 
% whenever there is a violoation, the entry in question is swapt with the previous entry, 
% and the two row’s diagonal entries are adjusted to ensure a row sum of 1.0
%
% M      :  migration matrix
% per    : 100 or 1 (depending on the matrix values range)
% upper  :  1 process the upper triangle (default 0) 

if nargin < 3
    upper = 0;
end

% extract default column
D = M(:,end);
M(:,end) = [];

if upper
    % rotates the upper triangle such matrix can be process 
    % top-bottom left-right
    M= rot90(M,2);
    D = rot90(D,2);
end

% loop through rows then by columns
for j=1:size(M,2)
    for i=j:size(M,1)-1
        % compare values
        if M(i+1,j) > M(i,j)
            % swap values
            temp = M(i,j);
            M(i,j) = M(i+1,j);
            M(i+1,j) = temp;
            % rescale 2 rows diagonal entries 
            % such that these rows sum up to (1 or 100)
            rw = [i i+1];            
            digpos = rw + (rw-1) * size(M,1);
            sum_rows = sum(M(rw,:),2);
            M(digpos) = M(digpos) + (per-sum_rows-D(rw))';
        end      
    end
end

if upper
    % revert rotation
    M= rot90(M,2);
    D = rot90(D,2);
end

% reinsert default column
M = [M D];

end

function nok = check_columns(M)
% this function  checks whether there are still violoation on the columns,
% that is if values on each side of the diagonal entry do not decrease
% monotonically
%
% Input : M transition matrix (NxP) including default column
% Output: nok: true if matrix does not satisfy monoticity column condition
% Author: Fabian Ojeda

% omit default column
Mt = M(:,1:end-1);

% check lower triangle - for non-mononotically decreasing elements
nok1 = any(any(tril(diff(tril(Mt),1,1),-1)  > 0));
% check upper triangle - for non-mononotically decreasing elements
nok2 = any(any(tril(diff(rot90(triu(Mt),2),1,1),-1)  > 0));
nok = nok1 | nok2;

end

function nok = check_rows(M)
% this function  checks whether there are still violoation on the rows,
% that is if values on each side of the diagonal entry do not decrease
% monotonically
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




% % % loop upper part matrix
% % for i=1:size(M,1)
% %     for j = i:D-2
% %         % compare values
% %         if M(i,j+1) > M(i,j)
% %             M(i,j+1) = M(i,j);
% %         end
% %     end
% % end
% % 
% % % loop lower matrix (starting from last row)
% % k=1;
% % for i=size(M,1):-1:1
% %     for j=D-k:-1:2
% %         %  compare values
% %         if M(i,j-1) > M(i,j)
% %             M(i,j-1) = M(i,j);
% %         end
% %     end
% %     k=k+1;
% % end
% % 
% % %---------------------------------------------------------
% % %3: column values should decrease on each side of diagonal
% % %---------------------------------------------------------
% % 
% % % loop lower part
% % for i=1:size(M,1)-1
% %     for j=1:i
% %         % compare values
% %         if M(i+1,j) > M(i,j)
% %             % swap values
% %             temp = M(i,j);
% %             M(i,j) = M(i+1,j);
% %             M(i+1,j) = temp;
% %             % rescale 2 rows to sum up to per (1 or 100)
% %             rw = [i i+1];
% %             sum_rows = sum(M(rw,ndf),2);
% %             M(rw,ndf) = bsxfun(@times,M(rw,ndf),((per-M(rw,D))./sum_rows));
% %         end      
% %     end
% % end
% % 
% % % loop upper part starting from
% % for i=1:size(M,1)-1
% %     for j=1:i
% %         % compare values
% %         if M(i+1,j) > M(i,j)
% %             % swap values
% %             temp = M(i,j);
% %             M(i,j) = M(i+1,j);
% %             M(i+1,j) = temp;
% %             % rescale 2 rows to sum up to per (1 or 100)
% %             rw = [i i+1];
% %             sum_rows = sum(M(rw,ndf),2);
% %             M(rw,ndf) = bsxfun(@times,M(rw,ndf),((per-M(rw,D))./sum_rows));
% %         end      
% %     end
% % end
% % 
% % 
% % 
% % mask = sign(bsxfun(@minus,(1:size(M,2)),(1:size(M,1))'));
% % 
% % for i=1:size(M,1)
% %     for j=1:D-1
% %         offset = j-mask(i,j);
% %         trans = M(i,j)-M(i,offset);
% %         if trans > 0
% %             % replace with last rating
% %             M(i,j) = M(i,offset);
% %         end
% %     end
% % end
% % 
% % %-------------------------------------------------------
% % %3: column values should decrease on each side of diagonal
% % %-------------------------------------------------------
% % for j=1:D-1
% %     for i=1:size(M,1)
% %         offset = i+mask(i,j);
% %         trans = M(i,j)-M(offset,j);
% %         if trans > 0
% %             % swap values
% %             temp = M(offset,j);
% %             M(offset,j) = M(i,j);
% %             M(i,j) = temp;
% %             % rescale the two rows to 1
% %             M(i,1:D-1) = (M(i,1:D-1)./sum(M(i,1:D-1))).*(per-M(i,D));
% %             M(offset,1:D-1) = (M(offset,1:D-1)./sum(M(offset,1:D-1))).*(per-M(offset,D));
% %         end
% %     end
% % end
% % 
% % %% rescale all rows to sumup to 1
% % sum_row = sum(M(:,1:D-1),2);
% % M(:,1:D-1) = bsxfun(@times,M(:,1:D-1),((per-M(:,D))./sum_row));
% % 
% % end