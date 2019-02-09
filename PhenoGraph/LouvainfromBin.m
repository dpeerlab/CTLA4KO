function [c,Q,bestpartition] = LouvainfromBin( filename, numiters )
ps = '~/Documents/MATLAB/PhenoGraph/Louvain';
if isempty(ps)
    error('Please enter the path to Louvain C++ code')
end
filename = strrep( filename, '.bin', '' );
% begin
fprintf(1, 'MATLAB: calling convert:\n');
command = [ps '/convert -i ' filename '.bin -o ' filename '_graph.bin -w ' filename '_graph.weights' ];
fprintf(1,'%s\n', command );
system( command );
% run community detection
for iter = 1:numiters
    
    fprintf(1,'MATLAB: running community detection, ITERATION %i\n', iter );
    command = [ps '/community ' filename '_graph.bin -l -1 -v -w ' filename '_graph.weights > ' filename '.tree'];
    fprintf( 1, '%s\n', command );
    [s,r] = system( command );
    fprintf(1, '\n');
    % find each iteration's modularity
    fprintf( 1, 'MATLAB: modularity scores:\n' );
    q = find_modularity( r );
    
    % find number of lvevls
    command = [ps '/hierarchy ' filename '.tree' ];
    fprintf(1, '%s\n', command );
    [s,r] = system( command );
    fprintf(1, '\n' );
    
    r = strtok(r,10);
    r = regexprep( r, 'Number of levels: ', '' );
    num_levels = str2num( r )-1;
    
    fprintf( 1, 'MATLAB: max level is %d\n', num_levels );
    
    % import levels
    for level = 1:num_levels
        fprintf( 1, 'MATLAB: importing level %d\n', level )
        command = [ps '/hierarchy ' filename '.tree -l ' num2str( level ) ' > ' filename '.tmp' ];
        fprintf(1, '%s\n', command );
        system( command );
        hierarchy_output = load( [filename '.tmp'] );
        c{iter,level} = hierarchy_output(:,2) + 1;
        Q{iter,level} = q(level);
    end
    
end
% find best partition
maxmod = 0;
for i = 1:numel(Q)
    if Q{i} > maxmod
        maxmod = Q{i};
        [I,J] = ind2sub( size(Q), i );
    end
end
bestpartition = c{I,J};
% clean up
files{1} = dir([filename '*.tmp']);
files{2} = dir([filename '*.tree']);
files{3} = dir([filename '*.weights']);
files{4} = dir([filename '*.bin']);
for i = 1:4
    for j = 1:length(files{i})
        delete( files{i}(j).name )
    end
end
%----------------------
function Q = find_modularity( r )
% Q = find_modularity( r )
% convert the text output into modularity score of each iteration
signature = '  modularity increased from %f to %f';
idx = 0;
while( ~isempty( r ) )
    % read a line and match it to the signature
    [token, r] = strtok( r, char( 10 ) );
    a = sscanf( token, signature );
    
    if( ~isempty( a ) )
        % signature matched copy new modularity
        idx = idx + 1;
        Q( idx ) = a( 2 );
    end
end
