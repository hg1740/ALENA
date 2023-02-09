classdef (ConstructOnLoad) Searchable < matlab.mixin.SetGet
    %
    % Searchable provides functionality allowing object arrays to be searched for user-specified content.
    %
    % Examples: (where fcn is a constructor of a class that implements Searchable)
    %
    %  obj = fcn('Name', 'Peter', 'Type', 'Employee', 'Description', {'Peter is an employee'});
    %  obj(2) = fcn('Name', 'Paul', 'Type', 'Manager', 'Description', {'Paul is a manager'});
    %
    %  res = find(obj, 'Type', 'Manager')
    %
    % Various modifiers can be specified in varargin to control the type of comparisons:
    %  '-nocase' makes string comparisons case-insensitive, default is case-sensitive
    %  '-partial' makes string comparisons based on the shorter of pairs of strings, default is full compare
    %  '-noignorenan' makes numeric comparisons treat nans as different, default is to treat them as same
    
    properties
        
        %None
        
    end
    
    methods % construction / destruction
        
        function obj = Searchable(varargin)
        
            %Nothing to do
            
        end
        
    end
    
    methods (Sealed) % Helper functions acting on a heterogrenous array of Searchables
        
        function varargout = findall(obj, varargin)
            %
            % Extends generic find facility to work with a hierarchical array
            %  of Collectable objects
            
            %Valid startring point ?
            assert(isa(obj, 'mvc.mixin.Collectable'), 'findall only applicable to a Collection class');
            
            %Flatten the collection, and pass it on
            [varargout{1:nargout}] = find(flatlist(obj), varargin{:});
            
        end
            
        function ret = find(obj, varargin)
            %
            % Generic find facility, looks for user specified content in object array
            %
            % TODO: this will need to grow considerably (e.g. regexp, case, etc etc)
            
            %Look for special tokens
            [return_option, logical_option, case_option, partial_option, nan_option, varargin] = i_specialTokens(varargin{:});
            
            %Ensure that anything left is in pairs
            assert(mod(numel(varargin), 2) == 0, 'search criteria must be specified in [param, value] pairs');
            
            %Start from this mask
            [bMask, varargin] = i_initialiseMask(obj, logical_option, varargin{:});
            
            %The function used to decide if two values are equal
            fcn = @(x,y)i_isequal(x, y, case_option, partial_option, nan_option);
            
            %For each property specified in inputs
            for i = 1:2:numel(varargin)
                
                %What logcial combination ?
                if strcmp(logical_option, '-and')
                    
                    %Only need to check things that we've not yet eliminated
                    % (conversion to logical makes masking the mask easier)
                    idx = find(bMask);
                    
                    %Which objects have the specified property ?
                    bPrp = arrayfun(@(x)isprop(x, varargin{i}), obj(idx));
                    
                    %Which of those objects have the specified value ?
                    bVal = arrayfun(@(x)fcn(x.(varargin{i}), varargin{i + 1}), obj(idx(bPrp)));
                    
                    %Hence update mask (carefully)
                    bMask(idx(bPrp)) = bVal;
                    
                    %Only continue if it's worth it
                    if ~any(bMask(:))
                        
                        %We've already eliminated everything, so we're done
                        break;
                    
                    end
                    
                else
                    
                    %Only need to check things that we've not yet included
                    % (conversion to logical makes masking the mask easier)
                    idx = find(~bMask);
                    
                    %Which objects have the specified property ?
                    bPrp = arrayfun(@(x)isprop(x, varargin{i}), obj(idx));
                    
                    %Which of those objects have the specified value ?
                    bVal = arrayfun(@(x)fcn(x.(varargin{i}), varargin{i + 1}), obj(idx(bPrp)));
                    
                    %Hence update mask (carefully)
                    bMask(idx(bPrp)) = bVal;
                    
                    if all(bMask(:))
                        
                        %We've already included everything, so we're done
                        break;
                        
                    end
                    
                end
                
            end
            
            %Send back what ?
            switch return_option
                
                case '-object'
                    
                    %The object(s) that match
                    ret = obj(bMask);
                    
                case '-mask'
                    
                    %The mask that matches
                    ret = bMask;
                    
                case '-index'
                    
                    %The indices that match
                    ret = find(bMask);
                    
                otherwise
                    error('bad return option');
            end
            
            function [bMask, remainder] = i_initialiseMask(obj, logical_option, varargin)
                
                %Initialise remainder
                remainder = varargin;
                
                %Start here
                switch logical_option
                    
                    case '-and'
                        
                        %Start by assuming everything is IN
                        bMask = true(size(obj));
                        
                    case '-or'
                        
                        %Start by assuming everything is OUT
                        bMask = false(size(obj));
                        
                    otherwise
                        error('bad option');
                end                        
                
                %Look for keywords that can be used to initialise mask
                bIsa = cellfun(@(x)ischar(x) && strcmp('isa', x), varargin);
                for j = find(bIsa)
                   
                    %Look for match
                    bMaskJ = arrayfun(@(x)isa(x, varargin{j + 1}), obj);

                    %Set mask accordingly
                    if strcmp(logical_option, '-and')                            
                        bMask = bMask & bMaskJ;
                    else
                        bMask = bMask | bMaskJ;
                    end
                                        
                end
            
                %Strip from remainder
                remainder(sort([find(bIsa), find(bIsa) + 1])) = [];
                
                %Look for keywords that can be used to initialise mask
                bClass = cellfun(@(x)ischar(x) && strcmp('class', x), varargin);
                for j = find(bClass)
                   
                    %Look for match
                    bMaskJ = arrayfun(@(x)strcmp(class(x), varargin{j + 1}), obj);

                    %Set mask accordingly
                    if strcmp(logical_option, '-and')                            
                        bMask = bMask & bMaskJ;
                    else
                        bMask = bMask | bMaskJ;
                    end
                                        
                end
            
                %Strip from remainder
                remainder(sort([find(bClass), find(bClass) + 1])) = [];
                
            end
            
            function [return_option, logical_option, case_option, partial_option, nan_option, remainder] = i_specialTokens(varargin)
                
                %Get return option
                [return_option, remainder] = i_option({'-object', '-mask', '-index'}, varargin{:});
                
                %Get logical option
                [logical_option, remainder] = i_option({'-and', '-or'}, remainder{:});    
                
                %Get case option
                [case_option, remainder] = i_option({'-case', '-nocase'}, remainder{:});    
                
                %Get partial option
                [partial_option, remainder] = i_option({'-nopartial', '-partial'}, remainder{:});    
                
                %Get nan option
                [nan_option, remainder] = i_option({'-ignorenan', '-noignorenan'}, remainder{:});    
                
                function [option, remainder] = i_option(options, varargin)
                    
                    %Initialise remainder
                    remainder = varargin;
                    
                    %Only interested in char within varargin
                    bdx = find(cellfun(@ischar, varargin));
                    
                    %Look for match
                    [b, jdx] = ismember(options, varargin(bdx));
                    
                    %Anything ?
                    if any(b)
                        
                        %Get whatever is chosen
                        option = options{b};
                        
                        %Eliminate from remainder
                        remainder(bdx(jdx(b))) = [];
                        
                    else
                        
                        %Go with first option as default
                        option = options{1};
                        
                    end
                    
                end
                
            end
            
            function b = i_isequal(x, y, case_option, partial_option, nan_option)
            
                %Case option requires some preprocessing
                if strcmpi(case_option, '-nocase')
                    x = i_lower(x);
                    y = i_lower(y);
                end
            
                %Partial option requires some preprocessing
                if strcmpi(partial_option, '-partial')
                    [x, y] = i_partial(x, y);
                end
                
                %Pass it on according to nan option
                if strcmpi(nan_option, '-ignorenan')
                    b = isequaln(x, y);
                else
                    b = isequal(x, y);
                end

                function [x,y] = i_partial(x,y)
                
                    %What have we got ?
                    if ischar(x) && ischar(y)

                        %Take the shorter of the two strings
                        n = min([numel(x), numel(y)]);
                        x = x(1:n);
                        y = y(1:n);
                       
                    else
                        %TODO: one string, one cell
                    end
                    
                end
                
                function x = i_lower(x)
                    if ischar(x)
                        x = lower(x);
                    elseif iscell(x)
                        b = cellfun(@ischar, x);
                        x(b) = cellfun(@i_lower, x(b), 'UniformOutput', false);
                    end
                end
                
            end
            
        end
        
        function [b, idx] = ismember(obj, list)
        
            %Initialise
            b = false(size(obj));
            idx = zeros(size(obj));
            
            %Check, carefully (avoids problems with builtin ismember, which throws error because eq is not sealed
            for i = 1:numel(obj)
                for j = 1:numel(list)
                    if obj(i) == list(j)
                        b(i) = true;
                        idx(i) = j;
                        break;
                    end
                end
            end
            
        end
        
    end
    
end
