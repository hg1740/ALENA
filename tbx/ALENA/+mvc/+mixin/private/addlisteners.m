function varargout = addlisteners(obj, varargin)
%
% addlisteners(obj, prp, fcn)
%
% Adds listener(s) to object(s) OBJ, such that changes (postset) to property(ies) PRP triggers function FCN.
%
% addlisteners(obj, fcn)
%
% Adds listener(s) to object(s) OBJ, such that changes (postset) to all (observable) properties triggers function FCN.

%Do it for each object individually (because addlistener is not sealed)
L = arrayfun(@(x)i_addlisteners(x, varargin{:}), obj, 'UniformOutput', false);

%What to send back ?
if nargout == 0
    
    %Nothing
    
elseif nargout == 1
    
    %Expand everything
    varargout{1} = horzcat(L{:});
    
elseif nargout == numel(obj)
    
    %Deal them out
    [varargout{1:nargout}] = deal(L{:});
    
else
    error('invalid return arguments');
end

    function L = i_addlisteners(obj, prp, fcn, nam)
        
        %If list of properties not explicitly supplied
        if nargin < 3
            
            %Then go with all
            fcn = prp;
            prp = [];
            
        end

        %If event name not specified
        if nargin < 4
            nam = 'PostSet';
        end
        
        %Limit to observables only
        mc = metaclass(obj);
        
        %Nothing specific ?
        if isempty(prp)
            
            %Listen to all
            prp = {mc.PropertyList([mc.PropertyList.SetObservable]).Name};
            
            %Check for dynamic properties 
            if isa(obj, 'mvc.mixin.Dynamicable')
                prp = [prp, obj.ObservableDynamicPropertyNames];
            end
            
        else
              
            names = {mc.PropertyList([mc.PropertyList.SetObservable]).Name};
            
            if isa(obj, 'mvc.mixin.Dynamicable')
                names = [names, obj.ObservableDynamicPropertyNames];
            end
            
            %Down-select
            prp = intersect(prp, names);

        end
        
        %Create listeners
        L = cellfun(@(x)addlistener(obj, x, nam, fcn), prp, 'UniformOutput', false);
        
        %Concatenate
        L = horzcat(L{:});
        
    end

end