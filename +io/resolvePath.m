function o = resolvePath(nwb, path)
dotTok = split(path, '.');
tokens = split(dotTok{1}, '/');
%skip first `/` if it exists
if isempty(tokens{1})
    tokens(1) = [];
end

%process slash tokens
o = nwb;
prefix = '';
for i=1:length(tokens)
    tok = tokens{i};
    errmsg = 'Could not resolve path `%s`.  Could not find `%s`.';
    if isa(o, 'types.untyped.Set')
        if any(strcmp(keys(o), tok))
            o = o.get(tok);
        else
            error(errmsg, path, tok);
        end
        continue;
    end
    %is class
    props = properties(o);
    if any(strcmp(props, tok))
        o = o.(tok);
    elseif any(strcmp(props, [prefix '_' tok]))
        o = o.([prefix '_' tok]);
        prefix = '';
    elseif any(startsWith(props, tok))
        if isempty(prefix)
            prefix = tok;
        else
            prefix = [prefix '_' tok];
        end
    else
        %dig one level into untyped sets because we don't know
        %if the untyped set is extraneous to the type or not.
        found = false;
        for j=1:length(props)
            set = o.(props{j});
            if isa(set, 'types.untyped.Set') &&...
                    any(strcmp(keys(set), tok))
                o = set.get(tok);
                found = true;
                break;
            end
        end
        if ~found
            error(errmsg, path, tok);
        end
    end
        [o, tokens] = resolveSet(o, tokens);
    else
        [o, tokens] = resolveObj(o, tokens);
    end
    if isempty(o)
        error(errmsg, path);
    end
end
end

function [o, remainder] = resolveSet(obj, tokens)
tok = tokens{1};
if any(strcmp(keys(obj), tok))
    o = obj.get(tok);
    remainder = tokens(2:end);
else
    o = [];
    remainder = tokens;
end
end

function [o, remainder] = resolveObj(obj, tokens)
props = properties(obj);
toklen = length(tokens);
eagerlist = cell(toklen,1);
for i=1:toklen
    eagerlist{i} = strjoin(tokens(1:i), '_');
end
% stable in this case preserves ordering with eagerlist bias
[eagers, ei, ~] = intersect(eagerlist, props, 'stable');
if isempty(eagers)
    % go one level down and check for sets
    proplen = length(props);
    issetprops = false(proplen, 1);
    for i=1:proplen
        issetprops(i) = isa(obj.(props{i}), 'types.untyped.Set');
    end
    setprops = props(issetprops);
    setpropslen = length(setprops);
    minlen = length(tokens) + 1;
    for i=1:setpropslen
        [new_o, new_tokens] = resolveSet(obj.(setprops{i}), tokens);
        new_toklen = length(new_tokens);
        if new_toklen < minlen
            o = new_o;
            remainder = new_tokens;
            minlen = new_toklen;
        end
    end
else
    o = obj.(eagers{end});
    remainder = tokens(ei(end)+1:end);
end
end