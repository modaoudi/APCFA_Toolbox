function [nGrps,gmns,gstds,gmin,gmax]  = get_clusters(data,thrsh,t)

    % This is an internal routine for the APCFA toolbox.
    % The main routine to execute the entire toolbox is APCFA_toolbox

    [A,~] = sort(data);                                                     % Sorting all mass defect points at a given nominal mass; The first value is zero
    if isempty(A) == 0
        if t == 0
            nGrps = sum(abs(diff(A)) > thrsh) + 1;                          % Number of groups with the set threshold
            disp(['Number of clusters detected:', num2str(nGrps)])
        else
            nGrps = sum(abs(diff(A)) > thrsh) + 1;                         % Number of groups with the set threshold
        end
        % if A contain only one element consider this
        if nGrps == 1 & size(A,1) == 1
            gmns  = mean(A);        % Compute the mean of the group
            gmin  = min(A);         % Compute the min of the group
            gmax  = max(A);         % Compute the max of the group
            gstds = std(A);         % Compute the std of the group
        elseif nGrps == 1 & size(A,1) > 1
            gmns  = mean(A);
            gmin  = min(A);
            gmax  = max(A); 
            gstds = std(A);
            % If A contains more than one element run cluster recognition
        else
            C = clusterdata(A,'maxclust',nGrps);             % Do the cluster with those groups C is index to grp; to work with this one At least two observation must be present
            cnts = accumarray(C,1);                          % The counts of each group
            gmns = accumarray(C,A,[],@mean);                 % Compute means of each group
            gmin  = accumarray(C,A,[],@min);                 % Compute the min of each group
            gmax  = accumarray(C,A,[],@max);                 % Compute the max of each group
            gstds = accumarray(C,A,[],@std);                 % Compute stds of each group

            % Sorting the groups
            [gmns, I] = sort(gmns);
            gstds = gstds(I);
            gmin = gmin(I);
            gmax = gmax(I);
        end
    else
        nGrps = 0;
        gmns  = 0;
        gmin  = 0;
        gmax  = 0;
        gstds = 0;
    end
end