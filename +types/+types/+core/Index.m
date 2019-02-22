classdef Index < types.core.NWBData
% Index Pointers that index data values


% PROPERTIES
properties
    target; % the target dataset that this index applies to
end

methods
    function obj = Index(varargin)
        % INDEX Constructor for Index
        %     obj = INDEX(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % target = ref to NWBData object
        varargin = [{'help' 'indexes into a list of values for a list of elements'} varargin];
        obj = obj@types.core.NWBData(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'target',[]);
        parse(p, varargin{:});
        obj.target = p.Results.target;
        if strcmp(class(obj), 'types.core.Index')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.target(obj, val)
        obj.target = obj.validate_target(val);
    end
    %% VALIDATORS
    
    function val = validate_data(obj, val)
    end
    function val = validate_target(obj, val)
        % Reference to type `NWBData`
        val = types.util.checkDtype('target', 'types.untyped.ObjectView', val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.NWBData(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.target)
            io.writeAttribute(fid, [fullpath '/target'], class(obj.target), obj.target, false);
        else
            error('Property `target` is required.');
        end
    end
end

end