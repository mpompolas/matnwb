classdef SpecFile < types.untyped.MetaClass
% SpecFile Dataset for storing contents of a specification file for either the core format or an extension.  Name should match name of file.`


% PROPERTIES
properties
    data; % property of type char
    help; % A help statement
    namespaces; % Namespaces defined in the file
end

methods
    function obj = SpecFile(varargin)
        % SPECFILE Constructor for SpecFile
        %     obj = SPECFILE(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % data = char
        % namespaces = char
        % help = char
        varargin = [{'help' 'Contents of format specification file.'} varargin];
        obj = obj@types.untyped.MetaClass(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'data',[]);
        addParameter(p, 'namespaces',[]);
        addParameter(p, 'help',[]);
        parse(p, varargin{:});
        obj.data = p.Results.data;
        obj.namespaces = p.Results.namespaces;
        obj.help = p.Results.help;
        if strcmp(class(obj), 'types.core.SpecFile')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.data(obj, val)
        obj.data = obj.validate_data(val);
    end
    function obj = set.help(obj, val)
        obj.help = obj.validate_help(val);
    end
    function obj = set.namespaces(obj, val)
        obj.namespaces = obj.validate_namespaces(val);
    end
    %% VALIDATORS
    
    function val = validate_data(obj, val)
        val = types.util.checkDtype('data', 'char', val);
    end
    function val = validate_help(obj, val)
        val = types.util.checkDtype('help', 'char', val);
    end
    function val = validate_namespaces(obj, val)
        val = types.util.checkDtype('namespaces', 'char', val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.untyped.MetaClass(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.help)
            io.writeAttribute(fid, [fullpath '/help'], class(obj.help), obj.help, false);
        end
        if ~isempty(obj.namespaces)
            io.writeAttribute(fid, [fullpath '/namespaces'], class(obj.namespaces), obj.namespaces, true);
        else
            error('Property `namespaces` is required.');
        end
    end
end

end