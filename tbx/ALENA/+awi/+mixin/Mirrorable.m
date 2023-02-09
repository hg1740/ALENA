classdef (ConstructOnLoad) Mirrorable < matlab.mixin.SetGet
    %Mirrorable Allows a part to be mirrored in a global axis plane.
    %
    % To construct a mirrored copy of another object this class must be
    % initiatied with the parameter '-Mirror' immediately followed by a
    % handle to the object that will be the 'mirrored-parent'.
    %   e.g. myLeftWing = awi.model.LiftingSurface('-Mirror', myRightWing)
    % 
    % 
    %
    % Useful for going from a half-aircraft model to a full-aircraft model
    % in one step. It would be extra nice if the mirrored object could be
    % dependent on its parent so that changes to the parent are
    % automatically propogated across.
    
    %Inernal references to mirrror properties and associated function
    %handles as well as any listener objects.
    properties (SetAccess = private, Hidden = true)
        %Store an internal list of properties that can be mirrored.
        %Listeners will be automatically generated for these properties
        %when a mirrored object is constructed.
        MirroredProperties = {};
        %Store an internal list of functions that will be invoked by the
        %listeners for the properties in 'MirroredProperties'.
        MirrorFunctions = {};
        %Store an internal list of any 'proplistener' objects 
        MirrorListeners = [];
    end
    
    %Mirror object and mirror plane
    properties (SetAccess = private)
        %The object that this object will mirror
        MirrorParent
        %The global plane that the object is being mirrored in.
        MirrorPlane = '';
    end
    
    properties (Constant, Hidden = true)
        %Default function that will be assigned to a mirror property
        DefaultMirrorFunc = @assignMirroredValue;
    end
    
    methods % set / get
        function set.MirrorParent(obj, val) %set.MirrorParent 
            %set.MirrorParent Set method for the property 'MirrorParent'.
            
            thisClass = class(obj);
            
            %Can only mirror objects of the same class            
            assert(strcmp(thisClass, class(val)), sprintf(['Can only '  , ...
                'mirror objects of the same class. Expected the object ', ...
                'to be of class ''%s''.'], thisClass));
            
            obj.MirrorParent = val;
            
        end
        function set.MirrorPlane(obj, val)  %set.MirrorPlane  
            %set.MirrorPlane Set method for the property 'MirrorPlane'.
            %
            % 'MirrorPlane' must be one of the following strings: 'XY',
            % 'XZ', 'YZ'.
            %
            % The value will automatically be converted to upper case.
            
            if ~isempty(val) %Always okay
                validatestring(upper(val), {'XY', 'XZ', 'YZ'}, ...
                    class(obj), 'MirrorPlane');
            end
            obj.MirrorPlane = upper(val);
            
        end
    end
    
    methods % construction
        function obj = Mirrorable(varargin)
            %Mirrorable Class constructor for the 'awi.mixin.Mirrorable'
            %class.                       
            
        end
    end
    
    methods (Access = protected) % mirroring the object
        function setMirrorBehaviour(obj, varargin)
            %setMirrorBehaviour Creates the listener objects for the
            %various mirror properties.
            
            if nargin < 2 %Escape route
                return
            end
            
            assert(isa(obj, 'awi.model.Component'), ['Can only mirror ', ...
                'objects that are of class ''awi.model.Component''.']);
            
            %Parse inputs
            [props2Listen, ListenFunc] = parseInputs(obj, varargin{:});   
           
            if isempty(props2Listen) || isempty(ListenFunc) %Escape route
                return
            end
                        
            %Are we already listnening to some of these properties? 
            %   - e.g. Because we have already called this method in a
            %          superclass constructor.
            if ~isempty(obj.MirrorListeners)
                disp('Wait');                
            end
            
            %Make the listeners
            lh = arrayfun(@(i) addlistener(obj.MirrorParent, ...
                props2Listen{i}, 'PostSet', ListenFunc{i}), ...
                1 : numel(props2Listen), 'Unif', false);
            obj.MirrorListeners = [obj.MirrorListeners, horzcat(lh{:})];
            
            %Mirror the object
            obj.mirrorElement;
            
            function [props2Listen, ListenFunc] = parseInputs(obj, varargin)
                %parseInputs Parses the 'varargin' cell array and pulls out
                %pertinent arguments.
                
                %Sensible defaults
                props2Listen = [];
                ListenFunc   = [];
                
                %Index 'varargin' to return just the parameters
                param = varargin(1 : 2 : end);
                
                if ~iscellstr(param) %Escape route
                    return
                end
                
                %Search 'varargin' for the special token '-Mirror'
                ind = find(ismember(param, '-Mirror'), 1);
                
                if isempty(ind) %Escape route
                    return
                end
                
                %If we get this far then the object is being initiated as a
                %mirrored copy of another object.
                
                %Pull the additonal arguments out of 'varargin'.
                MirrorObj = varargin{2 * ind};
                
                %Check for additional arguments...
                
                %Has the user defined a list of properties to listen to?
                ind = find(ismember(param, 'MirroredProperties'), 1);
                if isempty(ind)
                    props2Listen = obj.MirroredProperties;
                else
                    props2Listen = varargin{2 * ind};
                end                
                if isempty(props2Listen)  %Escape route
                    return
                end
                
                %Has the user defined a list of function handles associated
                %with the mirrored properties?
                ind = find(ismember(param, 'MirrorFunctions'), 1);
                if isempty(ind)
                    ListenFunc = obj.MirrorFunctions;
                else
                    ListenFunc = varargin{2 * ind};
                end                
                if isempty(ListenFunc) %Escape route
                    return
                end
                
                %Store a reference to the mirror parent object.
                obj.MirrorParent = MirrorObj;
                
                %Has the mirror-plane been defined?
                ind = find(ismember(param, 'MirrorPlane'), 1);
                if isempty(ind)
                    %If not, default to 'XZ'.
                    obj.MirrorPlane = 'XZ';
                else
                    obj.MirrorPlane = varargin{2 * ind};
                end
                                
            end
            
        end        
        function mirrorElement(obj)
            %mirrorThis Mirrors this object about the plane defined by
            %'obj.MirrorPlane'.
            %
            % At this level we cannot deal with the specifics of the
            % geometry and the 'mirror-plane'. Instead, we can make sure we
            % respect any Parent/Child relationships that the mirror parent
            % has established.
                       
            %Any information about the hierachy?
            if ~isa(obj, 'mvc.mixin.Collectable')
                return
            end
            
            %Parent?
            if isempty(obj.MirrorParent.Parent)
                return
            end
            
            %Mirrored objects have the same parent as their 'MirrorParent'
            obj.MirrorParent.Parent.add(obj);
            
            %Children?
            if isempty(obj.MirrorParent.Children)
                return
            end
            
        end
    end
    
    methods (Access = protected) % default callbacks for mirrored properties
        function assignMirroredValue(obj, src, evt)
            %assignMirroredValue Passes the mirrored value straight across.
        end
        function assignNegativeMirroredValue(obj, src, evt)
            %assignNegativeMirroredValue Assigns the negative of the
            %property value.
            
            
        end
    end
    
    methods (Access =  protected) % manipulating the 'mirrorable' properties
        function addMirrorProps(obj, propNames, funcHandles)
            %addMirrorProps Updates the list of mirrored properties and
            %mirror functions.
            
            if nargin ~= 3 %Escape route
                return
            end            
            
            %Check inputs
            assert(iscellstr(propNames), ['Expected ''propNames'' ', ...
                'to be a cell array of strings.']);
            assert(iscell(funcHandles), ['Expected ''funcHandles'' to ', ...
                'be a cell array of function handles']);
            assert(numel(propNames) == numel(funcHandles), ['Expected ' , ...
                'the number of mirrored properties to match the number ', ...
                'of function handles.']);
            
            %Replace any empty cells with the default function handles
            idx = cellfun(@isempty, funcHandles);
            funcHandles(idx) = {obj.DefaultMirrorFunc};
            idx = cellfun(@(x) isa(x, 'function_handle'), funcHandles);
            assert(all(idx), ['Expected ''funcHandles'' to be a cell ', ...
                'array of function handles.']);          
            
            %Append
            obj.MirroredProperties = [obj.MirroredProperties, propNames];
            obj.MirrorFunctions    = [obj.MirrorFunctions, funcHandles];
            
        end
    end
    
end

    