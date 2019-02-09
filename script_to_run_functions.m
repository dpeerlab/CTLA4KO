% This script demonstrates our use of several computational tools

% Note on usage:
% To be able to run all the functions, please make sure that they are on
% your MATLAB path.

% The PCHA code was taken from: http://www.mortenmorup.dk/MMhomepageUpdated_files/Page327.htm
% The code to find knee point is from: https://www.mathworks.com/matlabcentral/fileexchange/35094-knee-point

% Please contact Spencer Wei (spencer.wei@gmail.com) or Roshan Sharma
% (rs3380@columbia.edu) for any questions or additional scripts.

% We will assume that the data on which analysis is to be done as "data"

% Clustering using phenograph:
phenograph_clusters = phenograph(data, 30);

% Identifying archetypes: We will run the method to asking for number of
% archetypes from 2 to 25
% XC are the archetypes: Note the matrix may need to be transposed
% S, C are the convex coefficients
% SSE are the sum of sqaured errors

XC = cell(24, 1);
S = cell(24, 1);
C = cell(24, 1);
SSE = nan(24, 1);
for number_of_archetypes = 2:25
    j = number_of_archetypes - 1;
    [XC{j}, S{j}, C{j}, SSE(j)] = PCHA(data', number_of_archetypes, 1:size(data, 1), 1:size(data, 1), 0);
end

% Identifying the optimal number of archetypes
[res_x, idx_of_result] = knee_pt(SSE/sum(sum(data.^2)), 2:25, 1);

% Running diffusion components
G = get_diff_map(data', 'k', 30, 'ka', 10, 'n_eigs', 50);

% The output of G is a MATLAB struct type
% It has eigenvectors stored in G.EigenVecs, eigenvalues in G.EigenVals,
% the transition or Markov matrix in G.T. Please type G on MATLAB command
% window to see everything it stores.

% Diffusion distance
% t_opt is the optimal time
eigenVecs_new = G.EigenVecs;
eigenVals_new = G.EigenVals;
    
DM = @(t)(eigenVecs_new(:,2:end)*(sparse(diag(eigenVals_new(2:end).^t))));
DMt = DM(t_opt);

% This defines the diffusion distance matrix
