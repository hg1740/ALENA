classdef (ConstructOnLoad) Entity < ...
        mvc.model.DrawableThing & awi.mixin.Buildable & awi.mixin.DefaultBuild
    %Entity Base-class of all application specific objects comprising of an
    %array of generic objects that are:
    %
    % * Nameable: Each object in the object array has metadata that are
    %   manageable.
    %
    % * Searchable: The contents of the object array are searchable.
    %
    % * Contextable: The contents of the object array can present
    %   context-sensitive options to the user.
    %
    % * Collectable: Objects of this type can be arranged in a hierarchical
    %   fashion (i.e. with parent/child relationships).
    %
    % * Drawable : So that an object of class Entity can be drawn using a
    %   object-specific draw method.
    %
    % * UiTools: Instances of the object array can interact with the user
    %
    % * Debugable: To provide the user with quick and easy access to debug
    %   functions
    %
    % * Buildable: Allows the object to be defined using multiple
    %   parameterisation methods.
    %
    % * DefaultBuild: Allows a set of methods to be defined which return
    %   default objects. Useful for testing.
    %
    % The Entity class describes an object that does not have a position in
    % 3D space. It can be considered as an abstract quantity in the
    % conceptual sense, but it is not an abstract class.
    %
    %   For example: An Entity can be used to represent an element of the
    %   framework that does not need to be plotted, such as the
    %   awi.model.Framework class and other objects that are 'Collectable'.
                
    %Identifying the object
    properties (SetAccess = private)
        UniqueIdentifier
    end
    
    methods % construction / destruction
        
        function obj = Entity(varargin)
            
            %Pass it on
            obj@mvc.model.DrawableThing(varargin{:});
            obj@awi.mixin.Buildable;
            
            %Configure collectables - An Entity can collect other Entities
            obj.addCollectionSpec(@awi.model.Entity, [], [], [], true); % but the ability to do so is hidden
            
        end
        
    end
    
    methods (Hidden)
        function setUniqueIdentifier(obj, UUID, bReplace)
            %setUniqueIdentifier Assigns the character 
            
            if nargin < 2
                UUID = [];
            end
            if nargin < 3
                bReplace = false;
            end
            
            nObj = numel(obj);
            
            if isempty(UUID)
                error('Update code');
                %Make default UUIDs
                
                if bReplace
                    UUID = makeUUIDs(nObj);
                else
                    currentUUID = get(obj, {'UniqueIdentifier'});
                    bEmpty = cellfun(@isempty, currentUUID);
                    UUID = makeUUIDs(nnz(bEmpty));
                    obj = obj(~bEmpty);
                end
            end
            
            assert(and(iscellstr(UUID), nObj == numel(UUID)), ...
                ['Expected UUID to be a cell array of character ' , ...
                'vectors with the same number of elements as the ', ...
                'number of objects in the object array ''obj''.']);
            %Flatten the array
            if nnz(size(UUID) > 1) > 1
                UUID = reshape(UUID, [numel(UUID), 1]);
            end
            if isrow(UUID)
                UUID = UUID';
            end
            set(obj, {'UniqueIdentifier'}, UUID);
            
            function makeUUIDs
                % See https://uk.mathworks.com/matlabcentral/answers/459210-fast-uuid-generation-in-matlab
                uuid = java.util.UUID.randomUUID;
            end
            
        end
    end
end
