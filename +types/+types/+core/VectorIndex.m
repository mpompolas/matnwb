classdef VectorIndex < types.core.Index
% VectorIndex Pointers that index data values



methods
    function obj = VectorIndex(varargin)
        % VECTORINDEX Constructor for VectorIndex
        %     obj = VECTORINDEX(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        obj = obj@types.core.Index(varargin{:});
        if strcmp(class(obj), 'types.core.VectorIndex')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    
    %% VALIDATORS
    
    function val = validate_data(obj, val)
    end
    function val = validate_target(obj, val)
        % Reference to type `VectorData`
        val = types.util.checkDtype('target', 'types.untyped.ObjectView', val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.Index(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
    end
end

end