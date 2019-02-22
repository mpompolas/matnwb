classdef Image < types.core.NWBData
% Image Image base type.


% PROPERTIES
properties
    description; % description of image
    resolution; % pixels / cm
end

methods
    function obj = Image(varargin)
        % IMAGE Constructor for Image
        %     obj = IMAGE(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % description = char
        % resolution = float
        obj = obj@types.core.NWBData(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'description',[]);
        addParameter(p, 'resolution',[]);
        parse(p, varargin{:});
        obj.description = p.Results.description;
        obj.resolution = p.Results.resolution;
        if strcmp(class(obj), 'types.core.Image')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.description(obj, val)
        obj.description = obj.validate_description(val);
    end
    function obj = set.resolution(obj, val)
        obj.resolution = obj.validate_resolution(val);
    end
    %% VALIDATORS
    
    function val = validate_data(obj, val)
    end
    function val = validate_description(obj, val)
        val = types.util.checkDtype('description', 'char', val);
    end
    function val = validate_help(obj, val)
        val = types.util.checkDtype('help', 'char', val);
    end
    function val = validate_resolution(obj, val)
        val = types.util.checkDtype('resolution', 'float', val);
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
        if ~isempty(obj.resolution)
            io.writeAttribute(fid, [fullpath '/resolution'], class(obj.resolution), obj.resolution, false);
        else
            error('Property `resolution` is required.');
        end
    end
end

end