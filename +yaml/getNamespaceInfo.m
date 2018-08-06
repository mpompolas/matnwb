%returns sources given namespace path
function [filelist, name, dependencies] = getNamespaceInfo(namespacePath)
fid = fopen(namespacePath);
namespaceText = fread(fid, '*char')';
namecell = regexp(namespaceText, '\s+name:\s*(\S+)', 'tokens', 'once');
name = namecell{1};
filelist = misc.flattenTokens(regexp(namespaceText, '\s+-\s*source:\s*(\S+)', 'tokens'));
dependencies = misc.flattenTokens(regexp(namespaceText, '\s+-\s*namespace:\s*(\S+)', 'tokens'));
fclose(fid);
end