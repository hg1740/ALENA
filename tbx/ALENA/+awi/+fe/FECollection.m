classdef FECollection < awi.fe.FEBaseClass
    %FECollection Describes a collection of 'awi.fe' objects without having
    %to define the individual object handles. This saves on overheads when
    %thousands of nodes/elements are required to describe a model.
    %
    % Detailed Description
    %   - As a general rule, any 'awi.fe.FECollection' object will not
    %     aggregate objects that themselves are not subclasses of
    %     'awi.fe.FECollection', i.e. the object will not store any handles
    %     to other 'awi.fe' objects. This is to save on overheads when
    %     defining large meshes, however, this is at a loss of robustness
    %     and generality.
    
    %IDNumbers
    properties (Hidden)
        %Collections of 'aw.fe' objects can have multiple ID numbers 
        IDNumbers
        %Groups within the collection - to allow easy access and indexing
        SubGroups
    end
    
    methods % set / get
        function set.IDNumbers(obj, val) %set.IDNumbers
            validateIDvector(obj, val, 'IDNumbers');
            obj.IDNumbers = val;                
        end
        function set.SubGroups(obj, val) %set.SubGroups
            if isempty(val) %User trying to reset
                obj.SubGroups = [];
                return
            end
            validateattributes(val, {'cell'}, {'vector'}, class(obj), 'SubGroups');
            assert(mod(numel(val), 2) == 0, ['Expected ''SubGroups'' ', ...
                'to have an even number of elements.']);
            param = val(1 : 2 : end);
            value = val(2 : 2 : end);
            assert(iscellstr(param), ['Expected ''SubGroups'' to be ', ...
                'formatted as Param/Value cell-array.']);
            assert(all(cellfun(@isnumeric, value)), ['Expected ''SubGroups'' ', ...
                'values to be numeric, positive and integer-valued, ', ...
                'i.e. index numbers.'])
            idx = cellfun(@iscolumn, value);
            value(idx) = cellfun(@(x) x', value(idx), 'Unif', false);
            validateattributes(horzcat(value{:}), {'numeric'}, ...
                {'nonnegative', 'integer'}, class(obj), 'SubGroup index numbers');
            val(2 : 2 : end) = value;
            obj.SubGroups = val;
        end
    end
    
    methods (Access = protected) % validation
        function validateIDvector(obj, val, prp)
            if isempty(val)
                return
            end
            validateattributes(val, {'numeric'}, {'integer', 'row', ...
                'real', 'nonnan', 'nonnegative'}, class(obj), prp);
        end
    end
    
end

