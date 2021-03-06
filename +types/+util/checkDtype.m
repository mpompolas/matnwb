function val = checkDtype(name, type, val)
%ref
%any, double, int/uint, char
persistent WHITELIST;
if isempty(WHITELIST)
    WHITELIST = {...
        'types.untyped.ExternalLink'...
        'types.untyped.SoftLink'...
        };
end
if isstruct(type)
    names = fieldnames(type);
    assert(isstruct(val) || istable(val) || isa(val, 'containers.Map'), ...
        'types.untyped.checkDtype: Compound Type must be a struct, table, or a containers.Map');
    if (isstruct(val) && isscalar(val)) || isa(val, 'containers.Map')
        %check for correct array shape
        sizes = zeros(length(names),1);
        for i=1:length(names)
            if isstruct(val)
                subv = val.(names{i});
            else
                subv = val(names{i});
            end
            assert(isvector(subv),...
                ['types.util.checkDtype: struct of arrays as a compound type ',...
                'cannot have multidimensional data in their fields.  Field data ',...
                'shape must be scalar or vector to be valid.']);
            sizes(i) = length(subv);
        end
        sizes = unique(sizes);
        assert(isscalar(sizes),...
            ['struct of arrays as a compound type ',...
            'contains mismatched number of elements with unique sizes: [%s].  ',...
            'Number of elements for each struct field must match to be valid.'], ...
            num2str(sizes));
    end
    for i=1:length(names)
        pnm = names{i};
        subnm = [name '.' pnm];
        typenm = type.(pnm);
        
        if (isstruct(val) && isscalar(val)) || istable(val)
            val.(pnm) = types.util.checkDtype(subnm,typenm,val.(pnm));
        elseif isstruct(val)
            for j=1:length(val)
                elem = val(j).(pnm);
                assert(~iscell(elem) && ...
                    (isempty(elem) || ...
                    (isscalar(elem) || (ischar(elem) && isvector(elem)))),...
                    ['Fields for an array of structs for '...
                    'compound types should have non-cell scalar values or char arrays.']);
                val(j).(pnm) = types.util.checkDtype(subnm, typenm, elem);
            end
        else
            val(names{i}) = types.util.checkDtype(subnm,typenm,val(names{i}));
        end
    end
else
    errid = 'MATNWB:INVALIDTYPE';
    errmsg = ['Property `' name '` must be a ' type '.'];
    if isempty(val)
        return;
    end
    if isa(val, 'types.untyped.DataStub')
        %grab first element and check
        truval = val;
        if any(val.dims == 0)
            val = [];
        else
            val = val.load(1);
        end
    elseif isa(val, 'types.untyped.Anon')
        truval = val;
        val = val.value;
    else
        truval = [];
    end
    
    if any(strcmpi(type, {'single' 'double' 'logical' 'numeric'})) ||...
            startsWith(type, {'int' 'uint' 'float'})
        %all numeric types
        try
            val = types.util.correctType(val, type);
        catch ME
            error('MATNWB:CASTERROR', 'Could not cast type `%s` to `%s` for property `%s`',...
                class(val), type, name);
        end
    elseif strcmp(type, 'isodatetime')
        addpath(fullfile(fileparts(which('nwbfile')), 'external_packages', 'datenum8601'));
        assert(ischar(val) || iscellstr(val) || isdatetime(val) ||...
            (iscell(val) && all(cellfun('isclass', val, 'datetime'))), errid, errmsg);
        if ischar(val) || iscellstr(val)
            if ischar(val)
                val = {val};
            end
            
            datevals = cell(size(val));
            % one of:
            % +-hh:mm
            % +-hhmm
            % +-hh
            % Z
            tzre_pattern = '(?:[+-]\d{2}(?::?\d{2})?|Z)$';
            for i = 1:length(val)
                dnum = datenum8601(val{i});
                
                tzre_match = regexp(val{i}, tzre_pattern, 'once');
                if isempty(tzre_match)
                    tz = 'local';
                else
                    tz = val{i}(tzre_match:end);
                    if strcmp(tz, 'Z')
                        tz = 'UTC';
                    end
                end
                datevals{i} = ...
                    datetime(dnum(1), 'TimeZone', tz, 'ConvertFrom', 'datenum');
            end
            val = datevals;
        end
        
        if isdatetime(val)
            val = {val};
        end
        
        for i=1:length(val)
            if isempty(val{i}.TimeZone)
                val{i}.TimeZone = 'local';
            end
            val{i}.Format = 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSZZZZZ';
        end
        
        if isscalar(val)
            val = val{1};
        end
    elseif strcmp(type, 'char')
        assert(ischar(val) || iscellstr(val), errid, errmsg);
    else%class, ref, or link
        
        noncell = false;
        if ~iscell(val)
            val = {val};
            noncell = true;
        end
        for i=1:length(val)
            subval = val{i};
            if isempty(subval)
                continue;
            end
            
            if ~isa(subval, type) && ~any(strcmp(class(subval), WHITELIST))
                error(errid, errmsg);
            end
        end
        if noncell
            val = val{1};
        end
    end
    
    %reset to datastub/anon
    if ~isempty(truval)
        val = truval;
    end
end
end