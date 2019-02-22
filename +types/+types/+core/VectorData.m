classdef VectorData < types.core.NWBData
% VectorData Data values indexed by pointer


% PROPERTIES
properties
    description; % A short description of what these vectors are
end

methods
    function obj = VectorData(varargin)
        % VECTORDATA Constructor for VectorData
        %     obj = VECTORDATA(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % description = char
        varargin = [{'help' 'Values for a list of elements'} varargin];
        obj = obj@types.core.NWBData(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'description',[]);
        parse(p, varargin{:});
        obj.description = p.Results.description;
        if strcmp(class(obj), 'types.core.VectorData')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.description(obj, val)
        obj.description = obj.validate_description(val);
    end
    %% VALIDATORS
    
    function val = validate_data(obj, val)
    end
    function val = validate_description(obj, val)
        val = types.util.checkDtype('description', 'char', val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.NWBData(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.description)
            io.writeAttribute(fid, [fullpath '/description'], class(obj.description), obj.description, false);
        else
            error('Property `description` is required.');
        end
    end
end

end