classdef EventWaveform < types.core.NWBDataInterface
% EventWaveform Represents either the waveforms of detected events, as extracted from a raw data trace in /acquisition, or the event waveforms that were stored during experiment acquisition.


% PROPERTIES
properties
    spikeeventseries; % SpikeEventSeries object containing detected spike event waveforms
end

methods
    function obj = EventWaveform(varargin)
        % EVENTWAVEFORM Constructor for EventWaveform
        %     obj = EVENTWAVEFORM(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % spikeeventseries = SpikeEventSeries
        varargin = [{'help' 'Waveform of detected extracellularly recorded spike events'} varargin];
        obj = obj@types.core.NWBDataInterface(varargin{:});
        [obj.spikeeventseries, ivarargin] = types.util.parseConstrained(obj,'spikeeventseries', 'types.core.SpikeEventSeries', varargin{:});
        varargin(ivarargin) = [];
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        parse(p, varargin{:});
        if strcmp(class(obj), 'types.core.EventWaveform')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.spikeeventseries(obj, val)
        obj.spikeeventseries = obj.validate_spikeeventseries(val);
    end
    %% VALIDATORS
    
    function val = validate_spikeeventseries(obj, val)
        constrained = {'types.core.SpikeEventSeries'};
        types.util.checkSet('spikeeventseries', struct(), constrained, val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.NWBDataInterface(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.spikeeventseries)
            refs = obj.spikeeventseries.export(fid, fullpath, refs);
        end
    end
end

end