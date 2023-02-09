classdef Grid < matlab.mixin.SetGet
    
    % Properties
    properties
        id    = [];
        cs    = [];
        coord = [];
        cd    = [];
        ps    = [];
        seid  = [];
        part  = [];
        type  = [];
    end
    
    methods
        
        function obj = Grid(~)
            
        end
        
        function obj = set.id(obj,Value)
            obj.id     = Value;
        end
        
        
        function obj = set.cs(obj,Value)
            obj.cs     = Value;
        end
        
        
        function obj = set.coord(obj,Value)
            obj.coord    = Value;
        end
        
        
        function obj = set.cd(obj,Value)
            obj.cd     = Value;
        end
        
        
        function obj = set.ps(obj,Value)
            obj.ps     = Value;
        end
        
        
        function obj = set.seid(obj,Value)
            obj.seid    = Value;
        end
        
        function f = genNastranDeck(obj)
            f = {};
            f = cat(1,f,sprintf('GRID*   %16.0f%16.0f%16.9E%16.9E',obj.id,obj.cs,obj.coord(1),obj.coord(2)));
            f = cat(1,f,sprintf('*       %16.9E%16.0f%16.0f%16.0f',obj.coord(3),obj.cd,obj.ps,obj.seid));
        end
        
        function f = genNeocassDeck(obj)
            
            id_st = sprintf('%-8d', obj.id);
            cs_st = sprintf('%-8d', obj.cs);
            C1_st = sprintf('%-8.4f', obj.coord(1));
            C2_st = sprintf('%-8.4f', obj.coord(2));
            C3_st = sprintf('%-8.4f', obj.coord(3));
            
            f{1} = sprintf('GRID    %s%s%s%s%s', id_st,cs_st,strjust(C1_st,'right'),strjust(C2_st,'right'),strjust(C3_st,'right'));
        end
        
        %% Set field values of object
        % Define a set method for this class to make sure relevant
        % fields are set and whole object not redefined in a similar
        % way to structure behaviour
        function obj = setFields(obj,Value)
            
            if isa(Value,'Grid')
                obj = Value;
            elseif isobject(Value) || isstruct(Value)
                fieldsStruct = fieldnames(Value);
                fieldsObject = fieldnames(obj);
                fieldsObjectReduced = regexprep(fieldsObject,'Obj','');
                
                [~,iVal,iObj] = intersect(fieldsStruct,fieldsObjectReduced);
                
                for i = 1:numel(iObj)
                    obj.(fieldsObject{iObj(i)}) = Value.(fieldsStruct{iVal(i)});
                end
            end
            
        end
        
        function h = plot(obj,fig)
            
            if nargin == 2
                figure(fig);
                hold on;
            else
                figure;
            end
            h = plotobj(obj,'r','o');
        end
        
        function h = plotaero(obj,fig)
            
            if nargin == 2
                figure(fig);
                hold on;
            else
                figure;
            end
            obj1 = findobj(obj,'type','Aerodynamic');
            if ~isempty(obj1)
               h = plotobj(obj1,'b','o');
            end
        end
        
        function h = plotbeam(obj,fig)
            
            if nargin == 2
                figure(fig);
                hold on;
            else
                figure;
            end
            obj1 = findobj(obj,'type','Beam');
            if ~isempty(obj1)
               h = plotobj(obj1,'r','o');
            end
        end
                
    end
    
    methods (Access = private)
       
        function h = plotobj(obj,c,mk)
            
            h = hggroup;
            
            x = struct2mat(obj,'coord',1);
            y = struct2mat(obj,'coord',2);
            z = struct2mat(obj,'coord',3);
            
            plot3(x,y,z,'MarkerFaceColor',c,'Marker',mk,'MarkerEdgeColor','k','LineStyle','none','Parent',h);
            axis equal;
            
        end
        
    end
end
