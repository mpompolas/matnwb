classdef VoltageClampSeries < types.core.PatchClampSeries
% VoltageClampSeries Stores current data recorded from intracellular voltage-clamp recordings. A corresponding VoltageClampStimulusSeries (stored separately as a stimulus) is used to store the voltage injected.


% PROPERTIES
properties
    capacitance_fast; % Unit: Farad
    capacitance_fast_unit; % The base unit of measure used to store data. This should be in the SI unit. COMMENT: This is the SI unit (when appropriate) of the stored data, such as Volts. If the actual data is stored in millivolts, the field 'conversion' below describes how to convert the data to the specified SI unit.
    capacitance_slow; % Unit: Farad
    capacitance_slow_unit; % The base unit of measure used to store data. This should be in the SI unit. COMMENT: This is the SI unit (when appropriate) of the stored data, such as Volts. If the actual data is stored in millivolts, the field 'conversion' below describes how to convert the data to the specified SI unit.
    resistance_comp_bandwidth; % Unit: Hz
    resistance_comp_bandwidth_unit; % The base unit of measure used to store data. This should be in the SI unit. COMMENT: This is the SI unit (when appropriate) of the stored data, such as Volts. If the actual data is stored in millivolts, the field 'conversion' below describes how to convert the data to the specified SI unit.
    resistance_comp_correction; % Unit: %
    resistance_comp_correction_unit; % The base unit of measure used to store data. This should be in the SI unit. COMMENT: This is the SI unit (when appropriate) of the stored data, such as Volts. If the actual data is stored in millivolts, the field 'conversion' below describes how to convert the data to the specified SI unit.
    resistance_comp_prediction; % Unit: %
    resistance_comp_prediction_unit; % The base unit of measure used to store data. This should be in the SI unit. COMMENT: This is the SI unit (when appropriate) of the stored data, such as Volts. If the actual data is stored in millivolts, the field 'conversion' below describes how to convert the data to the specified SI unit.
    whole_cell_capacitance_comp; % Unit: Farad
    whole_cell_capacitance_comp_unit; % The base unit of measure used to store data. This should be in the SI unit. COMMENT: This is the SI unit (when appropriate) of the stored data, such as Volts. If the actual data is stored in millivolts, the field 'conversion' below describes how to convert the data to the specified SI unit.
    whole_cell_series_resistance_comp; % Unit: Ohm
    whole_cell_series_resistance_comp_unit; % The base unit of measure used to store data. This should be in the SI unit. COMMENT: This is the SI unit (when appropriate) of the stored data, such as Volts. If the actual data is stored in millivolts, the field 'conversion' below describes how to convert the data to the specified SI unit.
end

methods
    function obj = VoltageClampSeries(varargin)
        % VOLTAGECLAMPSERIES Constructor for VoltageClampSeries
        %     obj = VOLTAGECLAMPSERIES(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % capacitance_fast = float32
        % capacitance_fast_unit = char
        % capacitance_slow = float32
        % capacitance_slow_unit = char
        % resistance_comp_bandwidth = float32
        % resistance_comp_bandwidth_unit = char
        % resistance_comp_correction = float32
        % resistance_comp_correction_unit = char
        % resistance_comp_prediction = float32
        % resistance_comp_prediction_unit = char
        % whole_cell_capacitance_comp = float32
        % whole_cell_capacitance_comp_unit = char
        % whole_cell_series_resistance_comp = float32
        % whole_cell_series_resistance_comp_unit = char
        varargin = [{'capacitance_fast_unit' 'Farad' 'capacitance_slow_unit' 'Farad' 'help' 'Current recorded from cell during voltage-clamp recording' 'resistance_comp_bandwidth_unit' 'Hz' 'resistance_comp_correction_unit' 'pecent' 'resistance_comp_prediction_unit' 'pecent' 'whole_cell_capacitance_comp_unit' 'Farad' 'whole_cell_series_resistance_comp_unit' 'Ohm'} varargin];
        obj = obj@types.core.PatchClampSeries(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'capacitance_fast',[]);
        addParameter(p, 'capacitance_fast_unit',[]);
        addParameter(p, 'capacitance_slow',[]);
        addParameter(p, 'capacitance_slow_unit',[]);
        addParameter(p, 'resistance_comp_bandwidth',[]);
        addParameter(p, 'resistance_comp_bandwidth_unit',[]);
        addParameter(p, 'resistance_comp_correction',[]);
        addParameter(p, 'resistance_comp_correction_unit',[]);
        addParameter(p, 'resistance_comp_prediction',[]);
        addParameter(p, 'resistance_comp_prediction_unit',[]);
        addParameter(p, 'whole_cell_capacitance_comp',[]);
        addParameter(p, 'whole_cell_capacitance_comp_unit',[]);
        addParameter(p, 'whole_cell_series_resistance_comp',[]);
        addParameter(p, 'whole_cell_series_resistance_comp_unit',[]);
        parse(p, varargin{:});
        obj.capacitance_fast = p.Results.capacitance_fast;
        obj.capacitance_fast_unit = p.Results.capacitance_fast_unit;
        obj.capacitance_slow = p.Results.capacitance_slow;
        obj.capacitance_slow_unit = p.Results.capacitance_slow_unit;
        obj.resistance_comp_bandwidth = p.Results.resistance_comp_bandwidth;
        obj.resistance_comp_bandwidth_unit = p.Results.resistance_comp_bandwidth_unit;
        obj.resistance_comp_correction = p.Results.resistance_comp_correction;
        obj.resistance_comp_correction_unit = p.Results.resistance_comp_correction_unit;
        obj.resistance_comp_prediction = p.Results.resistance_comp_prediction;
        obj.resistance_comp_prediction_unit = p.Results.resistance_comp_prediction_unit;
        obj.whole_cell_capacitance_comp = p.Results.whole_cell_capacitance_comp;
        obj.whole_cell_capacitance_comp_unit = p.Results.whole_cell_capacitance_comp_unit;
        obj.whole_cell_series_resistance_comp = p.Results.whole_cell_series_resistance_comp;
        obj.whole_cell_series_resistance_comp_unit = p.Results.whole_cell_series_resistance_comp_unit;
        if strcmp(class(obj), 'types.core.VoltageClampSeries')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.capacitance_fast(obj, val)
        obj.capacitance_fast = obj.validate_capacitance_fast(val);
    end
    function obj = set.capacitance_fast_unit(obj, val)
        obj.capacitance_fast_unit = obj.validate_capacitance_fast_unit(val);
    end
    function obj = set.capacitance_slow(obj, val)
        obj.capacitance_slow = obj.validate_capacitance_slow(val);
    end
    function obj = set.capacitance_slow_unit(obj, val)
        obj.capacitance_slow_unit = obj.validate_capacitance_slow_unit(val);
    end
    function obj = set.resistance_comp_bandwidth(obj, val)
        obj.resistance_comp_bandwidth = obj.validate_resistance_comp_bandwidth(val);
    end
    function obj = set.resistance_comp_bandwidth_unit(obj, val)
        obj.resistance_comp_bandwidth_unit = obj.validate_resistance_comp_bandwidth_unit(val);
    end
    function obj = set.resistance_comp_correction(obj, val)
        obj.resistance_comp_correction = obj.validate_resistance_comp_correction(val);
    end
    function obj = set.resistance_comp_correction_unit(obj, val)
        obj.resistance_comp_correction_unit = obj.validate_resistance_comp_correction_unit(val);
    end
    function obj = set.resistance_comp_prediction(obj, val)
        obj.resistance_comp_prediction = obj.validate_resistance_comp_prediction(val);
    end
    function obj = set.resistance_comp_prediction_unit(obj, val)
        obj.resistance_comp_prediction_unit = obj.validate_resistance_comp_prediction_unit(val);
    end
    function obj = set.whole_cell_capacitance_comp(obj, val)
        obj.whole_cell_capacitance_comp = obj.validate_whole_cell_capacitance_comp(val);
    end
    function obj = set.whole_cell_capacitance_comp_unit(obj, val)
        obj.whole_cell_capacitance_comp_unit = obj.validate_whole_cell_capacitance_comp_unit(val);
    end
    function obj = set.whole_cell_series_resistance_comp(obj, val)
        obj.whole_cell_series_resistance_comp = obj.validate_whole_cell_series_resistance_comp(val);
    end
    function obj = set.whole_cell_series_resistance_comp_unit(obj, val)
        obj.whole_cell_series_resistance_comp_unit = obj.validate_whole_cell_series_resistance_comp_unit(val);
    end
    %% VALIDATORS
    
    function val = validate_capacitance_fast(obj, val)
        val = types.util.checkDtype('capacitance_fast', 'float32', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_capacitance_fast_unit(obj, val)
        val = types.util.checkDtype('capacitance_fast_unit', 'char', val);
    end
    function val = validate_capacitance_slow(obj, val)
        val = types.util.checkDtype('capacitance_slow', 'float32', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_capacitance_slow_unit(obj, val)
        val = types.util.checkDtype('capacitance_slow_unit', 'char', val);
    end
    function val = validate_resistance_comp_bandwidth(obj, val)
        val = types.util.checkDtype('resistance_comp_bandwidth', 'float32', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_resistance_comp_bandwidth_unit(obj, val)
        val = types.util.checkDtype('resistance_comp_bandwidth_unit', 'char', val);
    end
    function val = validate_resistance_comp_correction(obj, val)
        val = types.util.checkDtype('resistance_comp_correction', 'float32', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_resistance_comp_correction_unit(obj, val)
        val = types.util.checkDtype('resistance_comp_correction_unit', 'char', val);
    end
    function val = validate_resistance_comp_prediction(obj, val)
        val = types.util.checkDtype('resistance_comp_prediction', 'float32', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_resistance_comp_prediction_unit(obj, val)
        val = types.util.checkDtype('resistance_comp_prediction_unit', 'char', val);
    end
    function val = validate_whole_cell_capacitance_comp(obj, val)
        val = types.util.checkDtype('whole_cell_capacitance_comp', 'float32', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_whole_cell_capacitance_comp_unit(obj, val)
        val = types.util.checkDtype('whole_cell_capacitance_comp_unit', 'char', val);
    end
    function val = validate_whole_cell_series_resistance_comp(obj, val)
        val = types.util.checkDtype('whole_cell_series_resistance_comp', 'float32', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = fliplr(val.dims);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_whole_cell_series_resistance_comp_unit(obj, val)
        val = types.util.checkDtype('whole_cell_series_resistance_comp_unit', 'char', val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.PatchClampSeries(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.capacitance_fast)
            if startsWith(class(obj.capacitance_fast), 'types.untyped.')
                refs = obj.capacitance_fast.export(fid, [fullpath '/capacitance_fast'], refs);
            elseif ~isempty(obj.capacitance_fast)
                io.writeDataset(fid, [fullpath '/capacitance_fast'], class(obj.capacitance_fast), obj.capacitance_fast, false);
            end
        end
        if ~isempty(obj.capacitance_fast_unit) && ~isempty(obj.capacitance_fast)
            io.writeAttribute(fid, [fullpath '/capacitance_fast/unit'], class(obj.capacitance_fast_unit), obj.capacitance_fast_unit, false);
        end
        if ~isempty(obj.capacitance_slow)
            if startsWith(class(obj.capacitance_slow), 'types.untyped.')
                refs = obj.capacitance_slow.export(fid, [fullpath '/capacitance_slow'], refs);
            elseif ~isempty(obj.capacitance_slow)
                io.writeDataset(fid, [fullpath '/capacitance_slow'], class(obj.capacitance_slow), obj.capacitance_slow, false);
            end
        end
        if ~isempty(obj.capacitance_slow_unit) && ~isempty(obj.capacitance_slow)
            io.writeAttribute(fid, [fullpath '/capacitance_slow/unit'], class(obj.capacitance_slow_unit), obj.capacitance_slow_unit, false);
        end
        if ~isempty(obj.resistance_comp_bandwidth)
            if startsWith(class(obj.resistance_comp_bandwidth), 'types.untyped.')
                refs = obj.resistance_comp_bandwidth.export(fid, [fullpath '/resistance_comp_bandwidth'], refs);
            elseif ~isempty(obj.resistance_comp_bandwidth)
                io.writeDataset(fid, [fullpath '/resistance_comp_bandwidth'], class(obj.resistance_comp_bandwidth), obj.resistance_comp_bandwidth, false);
            end
        end
        if ~isempty(obj.resistance_comp_bandwidth_unit) && ~isempty(obj.resistance_comp_bandwidth)
            io.writeAttribute(fid, [fullpath '/resistance_comp_bandwidth/unit'], class(obj.resistance_comp_bandwidth_unit), obj.resistance_comp_bandwidth_unit, false);
        end
        if ~isempty(obj.resistance_comp_correction)
            if startsWith(class(obj.resistance_comp_correction), 'types.untyped.')
                refs = obj.resistance_comp_correction.export(fid, [fullpath '/resistance_comp_correction'], refs);
            elseif ~isempty(obj.resistance_comp_correction)
                io.writeDataset(fid, [fullpath '/resistance_comp_correction'], class(obj.resistance_comp_correction), obj.resistance_comp_correction, false);
            end
        end
        if ~isempty(obj.resistance_comp_correction_unit) && ~isempty(obj.resistance_comp_correction)
            io.writeAttribute(fid, [fullpath '/resistance_comp_correction/unit'], class(obj.resistance_comp_correction_unit), obj.resistance_comp_correction_unit, false);
        end
        if ~isempty(obj.resistance_comp_prediction)
            if startsWith(class(obj.resistance_comp_prediction), 'types.untyped.')
                refs = obj.resistance_comp_prediction.export(fid, [fullpath '/resistance_comp_prediction'], refs);
            elseif ~isempty(obj.resistance_comp_prediction)
                io.writeDataset(fid, [fullpath '/resistance_comp_prediction'], class(obj.resistance_comp_prediction), obj.resistance_comp_prediction, false);
            end
        end
        if ~isempty(obj.resistance_comp_prediction_unit) && ~isempty(obj.resistance_comp_prediction)
            io.writeAttribute(fid, [fullpath '/resistance_comp_prediction/unit'], class(obj.resistance_comp_prediction_unit), obj.resistance_comp_prediction_unit, false);
        end
        if ~isempty(obj.whole_cell_capacitance_comp)
            if startsWith(class(obj.whole_cell_capacitance_comp), 'types.untyped.')
                refs = obj.whole_cell_capacitance_comp.export(fid, [fullpath '/whole_cell_capacitance_comp'], refs);
            elseif ~isempty(obj.whole_cell_capacitance_comp)
                io.writeDataset(fid, [fullpath '/whole_cell_capacitance_comp'], class(obj.whole_cell_capacitance_comp), obj.whole_cell_capacitance_comp, false);
            end
        end
        if ~isempty(obj.whole_cell_capacitance_comp_unit) && ~isempty(obj.whole_cell_capacitance_comp)
            io.writeAttribute(fid, [fullpath '/whole_cell_capacitance_comp/unit'], class(obj.whole_cell_capacitance_comp_unit), obj.whole_cell_capacitance_comp_unit, false);
        end
        if ~isempty(obj.whole_cell_series_resistance_comp)
            if startsWith(class(obj.whole_cell_series_resistance_comp), 'types.untyped.')
                refs = obj.whole_cell_series_resistance_comp.export(fid, [fullpath '/whole_cell_series_resistance_comp'], refs);
            elseif ~isempty(obj.whole_cell_series_resistance_comp)
                io.writeDataset(fid, [fullpath '/whole_cell_series_resistance_comp'], class(obj.whole_cell_series_resistance_comp), obj.whole_cell_series_resistance_comp, false);
            end
        end
        if ~isempty(obj.whole_cell_series_resistance_comp_unit) && ~isempty(obj.whole_cell_series_resistance_comp)
            io.writeAttribute(fid, [fullpath '/whole_cell_series_resistance_comp/unit'], class(obj.whole_cell_series_resistance_comp_unit), obj.whole_cell_series_resistance_comp_unit, false);
        end
    end
end

end