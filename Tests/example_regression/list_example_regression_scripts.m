function scriptTable = list_example_regression_scripts()
%LIST_EXAMPLE_REGRESSION_SCRIPTS List selectable simulation scripts.
%
%   list_example_regression_scripts
%   scripts = list_example_regression_scripts

catalog = example_regression_catalog();
scriptTable = struct2table(rmfield(catalog, 'OutputNames'), 'AsArray', true);
scriptTable.OutputNames = string(cellfun(@(names) strjoin(names, '; '), ...
    {catalog.OutputNames}, 'UniformOutput', false)).';

if nargout == 0
    disp(scriptTable);
    clear scriptTable;
end
end
