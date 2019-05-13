classdef OpticalChannel < types.core.NWBContainer
% OpticalChannel One of possibly many groups storing channel-specific data COMMENT: Name is arbitrary but should be meaningful


% PROPERTIES
properties
    description; % Any notes or comments about the channel
    emission_lambda; % Emission wavelength for channel in nm
end

methods
    function obj = OpticalChannel(varargin)
        % OPTICALCHANNEL Constructor for OpticalChannel
        %     obj = OPTICALCHANNEL(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % description = char
        % emission_lambda = float
        varargin = [{'help' 'Metadata about an optical channel used to record from an imaging plane'} varargin];
        obj = obj@types.core.NWBContainer(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'description',[]);
        addParameter(p, 'emission_lambda',[]);
        parse(p, varargin{:});
        obj.description = p.Results.description;
        obj.emission_lambda = p.Results.emission_lambda;
        if strcmp(class(obj), 'types.core.OpticalChannel')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.description(obj, val)
        obj.description = obj.validate_description(val);
    end
    function obj = set.emission_lambda(obj, val)
        obj.emission_lambda = obj.validate_emission_lambda(val);
    end
    %% VALIDATORS
    
    function val = validate_description(obj, val)
        val = types.util.checkDtype('description', 'char', val);
    end
    function val = validate_emission_lambda(obj, val)
        val = types.util.checkDtype('emission_lambda', 'float', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.NWBContainer(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.description)
            if startsWith(class(obj.description), 'types.untyped.')
                refs = obj.description.export(fid, [fullpath '/description'], refs);
            elseif ~isempty(obj.description)
                io.writeDataset(fid, [fullpath '/description'], class(obj.description), obj.description, false);
            end
        else
            error('Property `description` is required.');
        end
        if ~isempty(obj.emission_lambda)
            if startsWith(class(obj.emission_lambda), 'types.untyped.')
                refs = obj.emission_lambda.export(fid, [fullpath '/emission_lambda'], refs);
            elseif ~isempty(obj.emission_lambda)
                io.writeDataset(fid, [fullpath '/emission_lambda'], class(obj.emission_lambda), obj.emission_lambda, false);
            end
        else
            error('Property `emission_lambda` is required.');
        end
    end
end

end