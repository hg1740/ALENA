classdef Conm2 < matlab.mixin.SetGet
    %CONM2 Defines a concentrated mass at a Grid point.
    %
    %  Concentrated mass similar to NASTRAN definition
    
    % Properties
    properties
        id     = [];    % Element identification number. (0 < Integer < 100,000,000)
        grid   = [];    % Grid point identification number. (Integer > 0)
        cid    = [];    % Coordinate system identification number
        m      = [];    % Mass value. (Real)
        offset = [];    % offset distances from the Grid point to the center of gravity of the mass in the coordinate system defined in cid
        i      = [];    % Mass moments of inertia measured at the mass center of gravity in the coordinate system defined by cid
        type   = [];
        part   = [];
        file   = [];
    end
    
    methods
        
        % Constructor
        function obj = Conm2(~)
            
        end
        
        function f = generateNastranDeck(obj)
            Inertia      = obj.i;
            idx          = abs(Inertia)<1e-6;
            Inertia(idx) = 0;
            
            f = {};
            f = cat(1,f,sprintf('CONM2*  %16.0f%16.0f%16.0f%16.9E',obj.id,obj.grid,obj.cid,obj.m));
            f = cat(1,f,sprintf('*       %16.9E%16.9E%16.9E%s',0,0,0,blanks(16)));
            f = cat(1,f,sprintf('*       %16.9E%16.9E%16.9E%16.9E',Inertia(1,1),Inertia(2,1),Inertia(2,2),Inertia(3,1)));
            f = cat(1,f,sprintf('*       %16.9E%16.9E%s',Inertia(3,2),Inertia(3,3),blanks(32)));
        end
        
        function f = genNeocassDeck(obj)

            id   = sprintf('%-8d', obj.id);
            node = sprintf('%-8d', obj.grid);
            cid  = sprintf('%-8d', obj.cid);
            m    = sprintf('%-.2e', obj.m);        
            
            % Define Offsets in correct format
            Offset       = obj.offset;

            Offx = sprintf('%-8.4f',Offset(1));
            Offy = sprintf('%-8.4f',Offset(2));
            Offz = sprintf('%-8.4f',Offset(3));
                    
            % Define inertias in correct format
            Inertia      = obj.i;
            
            I11 = sprintf('%-8s',num2_8ch(Inertia(1,1)));
            I21 = sprintf('%-8s',num2_8ch(Inertia(2,1)));
            I22 = sprintf('%-8s',num2_8ch(Inertia(2,2)));
            I31 = sprintf('%-8s',num2_8ch(Inertia(3,1)));
            I32 = sprintf('%-8s',num2_8ch(Inertia(3,2)));
            I33 = sprintf('%-8s',num2_8ch(Inertia(3,3)));

            f{1} = sprintf('CONM2   %s%s%s%s%s%s%s',id,node,cid,m,Offx,Offy,Offz);
            f{2} = sprintf('        %s%s%s%s%s%s',I11,I21,I22,I31,I32,I33);
          
        end 
        
    end
end
