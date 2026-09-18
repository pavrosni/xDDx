function absolutePath = resolve_script_relative_path(scriptDirectory, inputPath)
%RESOLVE_SCRIPT_RELATIVE_PATH Resolve a user path to an absolute canonical path.

scriptDirectory = char(scriptDirectory);
inputPath = char(inputPath);

if isempty(inputPath)
    error('resolve_script_relative_path:EmptyPath', 'Path value cannot be empty.');
end

pathObject = java.io.File(inputPath);
if pathObject.isAbsolute()
    absolutePath = char(pathObject.getCanonicalPath());
else
    absolutePath = char(java.io.File(fullfile(scriptDirectory, inputPath)).getCanonicalPath());
end

end
