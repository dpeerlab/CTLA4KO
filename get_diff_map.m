function GD = get_diff_map(data, varargin)
% the data must be cells on the columns and features on the rows

% number of cells
N = size(data, 2);

if size(data, 1) < size(data, 2)
    data = data';
end


% set up default parameters
k = 15;
ka = 5;
epsilon = 1;
n_eigs = 10;
eigs_method = 'normal';
knn_method = 'linear';
num_iters = 10;
knn_graph = 0;


% get the input parameters
if ~isempty(varargin)
    for j = 1:length(varargin)
        % k nearest neighbor 
        if strcmp(varargin{j}, 'ka')
            ka = varargin{j+1};
        end
        % for knn-autotune
        if strcmp(varargin{j}, 'k')
            k = varargin{j+1};
        end
        % epsilon
        if strcmp(varargin{j}, 'epsilon')
            epsilon = varargin{j+1};
        end
        % number of eigen-vectors
        if strcmp(varargin{j}, 'n_eigs')
            n_eigs = varargin{j+1};            
        end
        % method to compute eigen-vectors
        if strcmp(varargin{j}, 'method')
            eigs_method = varargin{j+1};
        end                
        % method to compute eigen-vectors
        if strcmp(varargin{j}, 'knn_method')
            knn_method = varargin{j+1};            
        end
        % number of iterations for parallel
        if strcmp(varargin{j}, 'num_iters')
            num_iters = varargin{j+1};            
        end
        % check if a knn-graph is already provided
        if strcmp(varargin{j}, 'knn_graph')
            disp('Using provided knn-graph')
            idx = varargin{j+1};
            dist = varargin{j+2};
            knn_graph = 1;
        end                
        
    end
end

disp('Using parameters:')
disp(['k = ', num2str(k)])
disp(['ka = ', num2str(ka)])
disp(['epsilon = ', num2str(epsilon)])
disp(['n_eigs =	', num2str(n_eigs)])

s1 = size(data, 1);
s2 = size(data, 2);


disp(['Size of data is ', num2str(s1) , ' rows and ', num2str(s2), ' columns.'])

if knn_graph == 0
    disp(['Doing ' knn_method, ' knnsearch'])
    if strcmp(knn_method, 'linear')
        [idx, dist] = knnsearch(data, data, 'k', k);
    else
        [idx, dist] = parfor_knnsearch(data, k, 'num_iters', num_iters);
    end
end

if ka ~= 0
    disp('Adapting sigma')
    dist = bsxfun(@rdivide, dist, dist(:,ka));
end

disp('Constructing graph')
i = repmat((1:N)',1,size(idx,2));
i = i(:);
j = idx(:);
s = dist(:).^2;
if epsilon > 0
    W = sparse(i, j, s);
else
    W = sparse(i, j, ones(size(s))); % unweighted kNN graph
end

distance_graph = W;

disp 'Symmetrize distances'
W = W + W';

if epsilon > 0
    disp 'Computing kernel'
    [i,j,s] = find(W);
    i = [i; (1:N)'];
    j = [j; (1:N)'];
    s = [s./(epsilon^2); zeros(N,1)];
    s = exp(-s);
    W = sparse(i,j,s);
end

disp 'Markov normalization'
%P = bsxfun(@rdivide, W, sum(W,2)); % Markov normalization
Dinv = sparse(1:size(W, 1), 1:size(W, 1), 1./sum(W, 2));
P = Dinv * W;

disp 'Computing Eigen-Values and Eigen-Vectors'
if strcmp(eigs_method, 'normal')
    [V, D] = eigs(P, n_eigs, 'LM', struct('disp', 0,'tol',1e-4,'maxit',1000));
    % sort them
    v = diag(D);
    [~, idxs] = sort(v, 'descend');

    D = D(idxs, idxs);
    V = V(:, idxs);
    V = P*V;
    for k = 1:size(V,2)
        V(:,k) = V(:,k)/norm(V(:,k));
    end
    v = diag(D);
    

% elseif strcmp(eigs_method, 'svd')
%     [U, ~, ~] = svd(P, 0);
% 	U = bsxfun(@rdivide, U, U(:,1));
% 	mappedX = U(:,2:no_dims + 1);
% elseif strcmp(eigs_method, 'random')
%     [U, ~, ~] = randPCA(P, n_eigs);
%     U = bsxfun(@rdivide, U, U(:,1));
%     mappedX = U(:,2:no_dims + 1);
end

% collect all the output
GD.W = W;
GD.T = P;
GD.EigenVals = v;
GD.EigenVecs = V;
GD.distance_graph = distance_graph;
GD.parameters.knn = k;
GD.parameters.ka = ka;
GD.parameters.epsilon = epsilon;
GD.parameters.n_eigs = n_eigs;
GD.parameters.eigs_method = eigs_method;

disp('Done')
end

% function [idx, dist] = parfor_knnsearch(data, k, varargin)
%     num_iters = 10;
%     chunk_size = ceil(size(data, 1) / num_iters);
%     
%     if ~isempty(varargin)
%         for j = 1:length(varargin) 
%             if strcmp(varargin{j}, 'num_iters')
%                 chunk_size = ceil(size(data, 1) / num_iters);
%             end
%         end
%     end            
%     
%     temp_idx = cell(num_iters, 1);
%     temp_dist = cell(num_iters, 1);
%     start = 1;    
%     samp_ids = cell(num_iters, 1);
%     for ct = 1:num_iters
%         terminate = min(ct*chunk_size, size(data, 1));
%         samp_ids{ct} = start:terminate;
%         start = terminate + 1;
%     end
%     
%     parfor ct = 1:num_iters                
%         [temp_idx{ct}, temp_dist{ct}] = knnsearch(data, data(samp_ids{ct}, :), 'k', k);        
%     end
%     
%     idx = nan(size(data, 1), k);
%     dist = nan(size(data, 1), k);
%     start = 1;    
%     for ct = 1:num_iters
%         terminate = min(ct*chunk_size, size(data, 1));
%         idx(start:terminate, :) = temp_idx{ct};
%         dist(start:terminate, :) = temp_dist{ct};
%         start = terminate + 1;
%     end
% end

function [idx_final, dist_final] = parfor_knnsearch(data, k, varargin)
    num_iters = 10;
    chunk_size = ceil(size(data, 1) / num_iters);
    n = size(data, 1);
    if ~isempty(varargin)
        for j = 1:length(varargin) 
            if strcmp(varargin{j}, 'num_iters')
                chunk_size = ceil(size(data, 1) / num_iters);
            end
        end
    end            
    
    total_chunks = ceil( n / chunk_size );
    idx_final    = []; %nan(size(data, 1), k);
    dist_final   = []; %nan(size(data, 1), k);
    
    tic
    % iterate over submatrices
    parfor iter = 1:total_chunks
        
        from = 1+chunk_size*(iter-1);
		to = min( from + chunk_size - 1, n );
		rx = from:to;
		
		[ idx, d ] = knnsearch( data, data( rx, : ), 'k', k );		

		idx_final = [idx_final; idx];
        dist_final = [dist_final; d];
    end
    

end