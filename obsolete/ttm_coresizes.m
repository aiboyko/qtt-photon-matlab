function sz=ttm_coresizes(ttm_A,varargin)
%this function takes a tt_matrix and returns all shapes of all its cores as
% d-by-4 vector
    cA=core2cell(ttm_A);

    for i=1:numel(cA)
        sz(:,i)=[size(cA{i}) ones(1,4-numel(size(cA{i})))];
    end
    
    if nargin>1
        if strcmp('only_sizes',varargin)
            sz=sz([2 3 ],:);
        end
    end
end