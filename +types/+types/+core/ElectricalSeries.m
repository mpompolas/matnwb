classdef ElectricalSeries < types.core.TimeSeries
% ElectricalSeries Stores acquired voltage data from extracellular recordings. The data field of an ElectricalSeries is an int or float array storing data in Volts. TimeSeries::data array structure: :blue:`[num times] [num channels] (or [num_times] for single electrode).`


% PROPERTIES
properties
    electrodes; % the electrodes that this series was generated from
end

methods
    function obj = ElectricalSeries(varargin)
        % ELECTRICALSERIES Constructor for ElectricalSeries
        %     obj = ELECTRICALSERIES(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % electrodes = DynamicTableRegion
        varargin = [{'data_unit' 'volt' 'help' 'Stores acquired voltage data from extracellular recordings'} varargin];
        obj = obj@types.core.TimeSeries(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'electrodes',[]);
        parse(p, varargin{:});
        obj.electrodes = p.Results.electrodes;
        if strcmp(class(obj), 'types.core.ElectricalSeries')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.electrodes(obj, val)
        obj.electrodes = obj.validate_electrodes(val);
    end
    %% VALIDATORS
    
    function val = validate_data(obj, val)
        val = types.util.checkDtype('data', 'numeric', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[Inf], [Inf Inf]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_data_unit(obj, val)
        val = types.util.checkDtype('data_unit', 'char', val);
    end
    function val = validate_electrodes(obj, val)
        val = types.util.checkDtype('electrodes', 'types.core.DynamicTableRegion', val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.TimeSeries(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.electrodes)
            refs = obj.electrodes.export(fid, [fullpath '/electrodes'], refs);
        else
            error('Property `electrodes` is required.');
        end
    end
end

end