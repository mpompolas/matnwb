classdef Images < types.core.NWBDataInterface
% Images A NWBDataInterface for storing images that have some relationship


% PROPERTIES
properties
    description; % Description of images in this container
    image; % Images stored in this NWBDataInterface
end

methods
    function obj = Images(varargin)
        % IMAGES Constructor for Images
        %     obj = IMAGES(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % description = char
        % image = Image
        varargin = [{'help' 'A collection of images that have some meaningful relationship'} varargin];
        obj = obj@types.core.NWBDataInterface(varargin{:});
        [obj.image, ivarargin] = types.util.parseConstrained(obj,'image', 'types.core.Image', varargin{:});
        varargin(ivarargin) = [];
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'description',[]);
        parse(p, varargin{:});
        obj.description = p.Results.description;
        if strcmp(class(obj), 'types.core.Images')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.description(obj, val)
        obj.description = obj.validate_description(val);
    end
    function obj = set.image(obj, val)
        obj.image = obj.validate_image(val);
    end
    %% VALIDATORS
    
    function val = validate_description(obj, val)
        val = types.util.checkDtype('description', 'char', val);
    end
    function val = validate_image(obj, val)
        constrained = { 'types.core.Image' };
        types.util.checkSet('image', struct(), constrained, val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.NWBDataInterface(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.description)
            io.writeAttribute(fid, [fullpath '/description'], class(obj.description), obj.description, false);
        else
            error('Property `description` is required.');
        end
        if ~isempty(obj.image)
            refs = obj.image.export(fid, fullpath, refs);
        else
            error('Property `image` is required.');
        end
    end
end

end