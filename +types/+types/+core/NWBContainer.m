classdef NWBContainer < types.untyped.MetaClass
% NWBContainer The attributes specified here are included in all interfaces.


% PROPERTIES
properties
    help; % Short description of what this type of NWBContainer contains.
end

methods
    function obj = NWBContainer(varargin)
        % NWBCONTAINER Constructor for NWBContainer
        %     obj = NWBCONTAINER(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % help = char
        obj = obj@types.untyped.MetaClass(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'help',[]);
        parse(p, varargin{:});
        obj.help = p.Results.help;
        if strcmp(class(obj), 'types.core.NWBContainer')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.help(obj, val)
        obj.help = obj.validate_help(val);
    end
    %% VALIDATORS
    
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