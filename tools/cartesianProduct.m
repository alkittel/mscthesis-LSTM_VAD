function result = cartesianProduct(sets)
%%cartesianProduct
% TABLE = cartesianProduct({A, B, ...}) returns an array with all combinations
% of the vectors A, B, ..., etc.

c = cell(1, numel(sets));
[c{:}] = ndgrid( sets{:} );
result = cell2mat( cellfun(@(v)v(:), c, 'UniformOutput',false) );