function [d_gospa, x_to_y_assignment, decomposed_cost] = ...
    GOSPA_extended(x_mat, y_mat, p, c, alpha)
% AUTHOR: Abu Sajana Rahmathullah
%Modified by Angel Garcia-Fernandez to account for extended targets

n_ouput_arg=nargout;

checkInput();

nx = size(x_mat.x, 2); % no of points in x_mat
ny = size(y_mat.x, 2); % no of points in y_mat

% compute cost matrix
cost_mat = zeros(nx, ny);
for ix = 1:nx
    for iy = 1:ny
        x_pos=x_mat.x(1:2, ix);
        X=x_mat.X(:,:, ix);
        
        y_pos=y_mat.x(1:2, iy);
        Y=y_mat.X(:,:, iy);
        
       
        cost_mat(ix, iy) ...
             = min(GaussianWassersteinDistance(x_pos,X,y_pos,Y), c);       
        
%         cost_mat(ix, iy) ...
%             = min(computeBaseDistance(x_mat(:, ix), y_mat(:, iy)), c);
    end
end

% intialise output values
decomposed_cost     = struct( ...
    'localisation', 0, ...
    'missed',       0, ...
    'false',        0);


x_to_y_assignment   = [];
opt_cost            = 0;

dummy_cost = (c^p) / alpha; % penalty for the cardinality mismatch

% below, cost is negated to make it compatible with auction algorithm
if nx == 0 % when x_mat is empty, all entries in y_mat are false
    opt_cost              = -ny * dummy_cost;
    decomposed_cost.false = opt_cost;
else
    if ny == 0 % when y_mat is empty, all entries in x_mat are missed
        opt_cost               = -nx * dummy_cost;
        
        if(alpha==2)
            decomposed_cost.missed = opt_cost;
        end
    else % when both x_mat and y_mat are non-empty, use auction algorithm
        cost_mat = -(cost_mat.^p);
        [x_to_y_assignment, y_to_x_assignment, ~] ...
            = auctionAlgorithm(cost_mat, 20*(nx * ny));
        % use the assignments to compute the cost
        for ind = 1:nx
            if x_to_y_assignment(ind) ~= 0
                opt_cost = opt_cost + cost_mat(ind,x_to_y_assignment(ind));
                
                if(alpha==2)
                    
                    decomposed_cost.localisation = ...
                        decomposed_cost.localisation ...
                        + cost_mat(ind,x_to_y_assignment(ind)) ...
                        .* double(cost_mat(ind,x_to_y_assignment(ind)) > -c^p);
                    
                    decomposed_cost.missed      = decomposed_cost.missed ...
                        - dummy_cost ...
                        .* double(cost_mat(ind,x_to_y_assignment(ind)) == -c^p);
                    
                    decomposed_cost.false       = ...
                        decomposed_cost.false ...
                        - dummy_cost ...
                        .* double(cost_mat(ind,x_to_y_assignment(ind)) == -c^p);
                end
            else
                opt_cost               = opt_cost - dummy_cost;
                if(alpha==2)
                    decomposed_cost.missed = decomposed_cost.missed - dummy_cost;
                end
            end
        end
        opt_cost = opt_cost - sum(y_to_x_assignment == 0) * dummy_cost;
        if(alpha==2)
            decomposed_cost.false = decomposed_cost.false ...
                - sum(y_to_x_assignment == 0) * dummy_cost;
        end
    end
end

% final output
d_gospa                      = (-opt_cost)^(1/p);
decomposed_cost.localisation = (-decomposed_cost.localisation);
decomposed_cost.missed       = (-decomposed_cost.missed);
decomposed_cost.false        = (-decomposed_cost.false);

 function checkInput()
        if size(x_mat, 1) ~= size(y_mat, 1)
            error('The number of rows in x_mat & y_mat should be equal.');
        end
        if ~((p >= 1) && (p < inf))
            error('The value of exponent p should be within [1,inf).');
        end
        if ~(c>0)
            error('The value of base distance c should be larger than 0.');
        end
        
        if ~((alpha > 0) && (alpha <= 2))
            error('The value of alpha should be within (0,2].');
        end
        if alpha ~= 2 && n_ouput_arg==3
            warning(['decomposed_cost is not valid for alpha = ' ...
                num2str(alpha)]);
        end
    end
end
