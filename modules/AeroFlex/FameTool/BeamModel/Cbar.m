classdef Cbar < matlab.mixin.SetGet
    % CBAR Defines a simple beam element.
    %
    
    % Properties
    properties
        id     = [];    % Unique element identification number. (0 < Integer < 100,000,000)
        pid    = [];    % Property identification number of a PBAR element
        conn   = [];    % Grid point identification numbers of connection points (2 unique values)
        barg0  = [];    % Alternate method to supply the orientation vector ? using Grid point G0. The direction of ? is from NA to G0. ? is then translated to NA
        orient = [];    % Components of orientation vector ?, from NA, in the displacement coordinate system at NA
        offset = 'GGG'; % Offset identification (see Table 5.7), see NEOCASS manual
        part   = [];
        R      = [];
        D      = [];
        M      = [];
        colloc = [];
        
    end
    
    methods
        
        function obj = Cbar(~)
            
            % Should evaluate the following upon construction?
            % collocation point     |
            % Rotation matrix       | This is currently in the
            % Stiffness matrix      | BeamModel class definition
            % Mass matrix           |
            
        end
        
        function f = genNastranDeck(obj,ff)
            if nargin <2
                ff = 8;
            end
            
            if ff == 16
                %16 field format definition
                f = {};
                f = cat(1,f,sprintf('CBAR*   %16.0f%16.0f%16.0f%16.0f',obj.id,obj.pid,obj.conn(1),obj.conn(2)));
                f = cat(1,f,sprintf('*       %16.9E%16.9f%16.9f%s',obj.orient(1),obj.orient(2),obj.orient(3),obj.offset));
            end
            
            if ff == 8
                f = {};
                f = cat(1,f,sprintf('CBAR    %-8i%8i%8i%8i%8s%8s%8s%8s', obj.id, obj.pid, obj.conn(1), ...
                    obj.conn(2), num2_8ch(obj.orient(1)),num2_8ch(obj.orient(2)),num2_8ch(obj.orient(3)),obj.offset));
            end
        end
        
        function f = genNeocassDeck(obj)
            f = {};
            f = cat(1,f,sprintf('CBAR    %-8i%8i%8i%8i%8s%8s%8s%8s', obj.id, obj.pid, obj.conn(1), ...
                obj.conn(2), num2_8ch(obj.orient(1)),num2_8ch(obj.orient(2)),num2_8ch(obj.orient(3)),obj.offset));
        end
        
        
        
    end
end
