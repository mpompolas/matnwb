classdef DynamicTable < types.core.NWBDataInterface
% DynamicTable A group containing multiple datasets that are aligned on the first dimension (Currently, this requirement if left up to APIs to check and enforce). Apart from a column that contains unique identifiers for each row there are no other required datasets. Users are free to add any number of VectorData objects here. Table functionality is already supported through compound types, which is analogous to storing an array-of-structs. DynamicTable can be thought of as a struct-of-arrays. This provides an alternative structure to choose from when optimizing storage for anticipated access patterns. Additionally, this type provides a way of creating a table without having to define a compound type up front. Although this convenience may be attractive, users should think carefully about how data will be accessed. DynamicTable is more appropriate for column-centric access, whereas a dataset with a compound type would be more appropriate for row-centric access. Finally, data size should also be taken into account. For small tables, performance loss may be an acceptable trade-off for the flexibility of a DynamicTable. For example, DynamicTable was originally developed for storing trial data and spike unit metadata. Both of these use cases are expected to produce relatively small tables, so the spatial locality of multiple datasets present in a DynamicTable is not expected to have a significant performance impact. Additionally, requirements of trial and unit metadata tables are sufficiently diverse that performance implications can be overlooked in favor of usability.


% PROPERTIES
properties
    colnames; % The names of the columns in this table. This should be used to specifying an order to the columns
    description; % Description of what is in this dynamic table
    id; % The unique identifier for the rows in this dynamic table
    vectordata; % The vector columns in this dynamic table
    vectorindex; % The indices for the vector columns in this dynamic table
end

methods
    function obj = DynamicTable(varargin)
        % DYNAMICTABLE Constructor for DynamicTable
        %     obj = DYNAMICTABLE(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % colnames = char
        % description = char
        % id = ElementIdentifiers
        % vectordata = VectorData
        % vectorindex = VectorIndex
        varargin = [{'help' 'A column-centric table'} varargin];
        obj = obj@types.core.NWBDataInterface(varargin{:});
        [obj.vectordata, ivarargin] = types.util.parseConstrained(obj,'vectordata', 'types.core.VectorData', varargin{:});
        varargin(ivarargin) = [];
        [obj.vectorindex, ivarargin] = types.util.parseConstrained(obj,'vectorindex', 'types.core.VectorIndex', varargin{:});
        varargin(ivarargin) = [];
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'colnames',[]);
        addParameter(p, 'description',[]);
        addParameter(p, 'id',[]);
        parse(p, varargin{:});
        obj.colnames = p.Results.colnames;
        obj.description = p.Results.description;
        obj.id = p.Results.id;
        if strcmp(class(obj), 'types.core.DynamicTable')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.colnames(obj, val)
        obj.colnames = obj.validate_colnames(val);
    end
    function obj = set.description(obj, val)
        obj.description = obj.validate_description(val);
    end
    function obj = set.id(obj, val)
        obj.id = obj.validate_id(val);
    end
    function obj = set.vectordata(obj, val)
        obj.vectordata = obj.validate_vectordata(val);
    end
    function obj = set.vectorindex(obj, val)
        obj.vectorindex = obj.validate_vectorindex(val);
    end
    %% VALIDATORS
    
    function val = validate_colnames(obj, val)
        val = types.util.checkDtype('colnames', 'char', val);
    end
    function val = validate_description(obj, val)
        val = types.util.checkDtype('description', 'char', val);
    end
    function val = validate_id(obj, val)
        val = types.util.checkDtype('id', 'types.core.ElementIdentifiers', val);
    end
    function val = validate_vectordata(obj, val)
        constrained = { 'types.core.VectorData' };
        types.util.checkSet('vectordata', struct(), constrained, val);
    end
    function val = validate_vectorindex(obj, val)
        constrained = { 'types.core.VectorIndex' };
        types.util.checkSet('vectorindex', struct(), constrained, val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.NWBDataInterface(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.colnames)
            io.writeAttribute(fid, [fullpath '/colnames'], class(obj.colnames), obj.colnames, true);
        else
            error('Property `colnames` is required.');
        end
        if ~isempty(obj.description)
            io.writeAttribute(fid, [fullpath '/description'], class(obj.description), obj.description, false);
        else
            error('Property `description` is required.');
        end
        if ~isempty(obj.id)
            refs = obj.id.export(fid, [fullpath '/id'], refs);
        else
            error('Property `id` is required.');
        end
        if ~isempty(obj.vectordata)
            refs = obj.vectordata.export(fid, fullpath, refs);
        end
        if ~isempty(obj.vectorindex)
            refs = obj.vectorindex.export(fid, fullpath, refs);
        end
    end
end

end