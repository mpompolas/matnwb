classdef NWBData < types.untyped.MetaClass
% NWBData The attributes specified here are included in all interfaces.


% PROPERTIES
properties
    data; % property of type any
    help; % Short description of what this type of Interface contains.
end

methods
    function obj = NWBData(varargin)
        % NWBDATA Constructor for NWBData
        %     obj = NWBDATA(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % data = any
        % help = char
        obj = obj@types.untyped.MetaClass(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'data',[]);
        addParameter(p, 'help',[]);
        parse(p, varargin{:});
        obj.data = p.Results.data;
        obj.help = p.Results.help;
        if strcmp(class(obj), 'types.core.NWBData')
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
    %% VALIDATORS
    
    function val = validate_data(obj, val)
    end
    function val = validate_help(obj, val)
        val = types.util.checkDtype('help', 'char', val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.untyped.MetaClass(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.help)
            io.writeAttribute(fid, [fullpath '/help'], class(obj.help), obj.help, false);
        else
            error('Property `help` is required.');
        end
    end
end

end