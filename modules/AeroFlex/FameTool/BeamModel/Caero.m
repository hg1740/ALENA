classdef Caero < matlab.mixin.SetGet
    
    % CAERO aerodynamic panel definition
    %
    % The CAERO card (see Table 5.30) is used to define all the geometric (see
    % Table 5.28) and mesh parameters (see Table 5.29 ) for each box of the
    % aerodynamic Model which are required by the meshing tool coming from
    % TORNADO Vortex Lattice code.
    %
    % See NEOCASS manual p.98 for detailed definition
    
    % Properties
    properties
        id          = [];   % Identification number
        startX      = [];   % X-position of starting corner point
        startY      = [];   % Y-position of starting corner point
        startZ      = [];   % Z-position of starting corner point
        twist1      = [];   % Twist at root chord side of panel
        twist2      = [];   % Twist at outboard side of panel
        c           = [];   % Root chord
        t           = [];   % Taper ratio of panel
        b           = [];   % Span of panel
        sw          = [];   % Sweep angle
        rootAirfoil = [];   % Name of root profile default={}
        tipAirfoil  = [];   % Name of tip profile default={}
        dih         = [];   % Dihedral angle
        meshType    = [];   % Mesh type
        cid         = [];   % Coordinate system identifier
        ny          = [];   % Span wise number of panels
        nx          = [];   % Chord wise number of panels
        flapped     = [];   % Flapped or not?
        fnx         = [];   % Number of chord wise number of control panels
        fc          = [];   % Portion of the chord taken up by the control surface
        csType      = [];   % Control surface type
        csId        = [];   % Control surface ID
        csData      = [];   % [1x4] array
        part        = [];
        csLink      = [];
        symetric    = [];
        %    csData{1} : 0 - no control surface, 1 - a control surface
        %    csData{2} : root chord
        %    csData{3} : tip chord
        %    csData{4} : number of chord wise panels on control surface
    end
    
    properties(Dependent = true)
        outline        = [];
        controloutline = [];
        sw_le          = [];
    end
    
    properties(Hidden = true)
        
        le      = [];
        te      = [];
    end
    
    methods
        
        function obj = Caero(~)
            
        end
        
        %% Set field values of object
        % Define a set method for this class to make sure relevant
        % fields are set and whole object not redefined in a similar
        % way to structure behaviour
        function obj = setFields(obj,Value)
            
            if isa(Value,'Caero')
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
        
        function value = get.outline(obj)
            
            %1--------2
            %|        |
            %4--------3
            
            %First point
            X1 = obj.startX;
            Y1 = obj.startY;
            Z1 = obj.startZ;
            
            %Second point
            X2 = X1 + 0.25 * obj.c + tan(pi*obj.sw/180) * obj.b - 0.25*obj.t * obj.c;
            if obj.dih == 90
                Y2 = Y1;
                Z2 = Z1 + obj.b;
            else
                Y2 = Y1 + obj.b;
                Z2 = Z1 + obj.b * tan(pi*obj.dih/180);
            end
            
            %Third point
            X3 = X2 + obj.t * obj.c;
            Y3 = Y2;
            Z3 = Z2;
            
            %Fourth point
            X4 = X1 + obj.c;
            Y4 = Y1;
            Z4 = Z1;
            
            value = [X1,Y1,Z1;X2,Y2,Z2;X3,Y3,Z3;X4,Y4,Z4];
            
        end
        
        function value = get.controloutline(obj)
            
            %1--------2
            %|        |
            %4--------3
            
            value = [];
            
            if ~isempty(obj.csData)
                
                %First point
                X10 = obj.startX;
                X1 = obj.startX + (1-obj.csData(2)) * obj.c;
                Y1 = obj.startY;
                Z1 = obj.startZ;
                
                %Second point
                X20 = X10 + 0.25 * obj.c + tan(pi*obj.sw/180) * obj.b - 0.25*obj.t * obj.c;
                X2 =  X20 + (1-obj.csData(3)) * obj.t * obj.c;
                if obj.dih == 90
                    Y2 = Y1;
                    Z2 = Z1 + obj.b;
                else
                    Y2 = Y1 + obj.b;
                    Z2 = Z1 + obj.b * tan(pi*obj.dih/180);
                end
                
                %Third point
                X3 = X20 + obj.t * obj.c;
                Y3 = Y2;
                Z3 = Z2;
                
                %Fourth point
                X4 = X10 + obj.c;
                Y4 = Y1;
                Z4 = Z1;
                
                value = [X1,Y1,Z1;X2,Y2,Z2;X3,Y3,Z3;X4,Y4,Z4];
            end
            
        end
        
        function value = get.sw_le(obj)
            %1--------2
            %|        |
            %4--------3
            
            %First point
            X1 = obj.startX;
            
            %Second point
            X2 = X1 + 0.25 * obj.c + tan(pi*obj.sw/180) * obj.b - 0.25*obj.t * obj.c;
            
            % Extract the leading edge sweep
            value = 180*atan((X2-X1)/obj.b)/pi;
        end
        
        function f = genNeocassDeck(obj)
            id_str     = sprintf('%-8d', obj.id);
            dih_str    = sprintf('%-8.5f', obj.dih);
            cid_str    = sprintf('%-8d', obj.cid);
            ny_str     = sprintf('%-8d', obj.ny);
            nx_str     = sprintf('%-8d', obj.nx);
            Raf_str    = sprintf('%-8s', obj.rootAirfoil);
            Taf_str    = sprintf('%-8s', obj.tipAirfoil);
            Mt_str     = sprintf('%-8g', obj.meshType);
            cx_str     = sprintf('%-8.4f', obj.startX);
            cy_str     = sprintf('%-8.4f', obj.startY);
            cz_str     = sprintf('%-8.4f', obj.startZ);
            Chord_str  = sprintf('%-8g', obj.c);
            Span_str   = sprintf('%-8.4f', obj.b);
            Taper_str  = sprintf('%-8.5f', obj.t);
            Sweep_str  = sprintf('%-8.4f', obj.sw);
            twist1_str = sprintf('%-8.5f', obj.twist1);
            twist2_str = sprintf('%-8.5f', obj.twist2);
            
            f{1,1} = sprintf('CAERO1  %s%s%s%s%s%s%s%s%s', id_str, dih_str, ...
                cid_str, ny_str, nx_str, Raf_str, Taf_str, Mt_str);
            f{2,1} = sprintf('        %s%s%s%s%s%s%s%s%s%s', cx_str, cy_str, ...
                cz_str, Chord_str, Span_str, Taper_str, Sweep_str, twist1_str, twist2_str);
            
            if ~isempty(obj.csType)
                csDum = sprintf('%-8g', 1);
                csBeg = sprintf('%-8.5f', obj.csData(2));
                csEnd = sprintf('%-8.5f', obj.csData(3));
                csNx  = sprintf('%-8d', obj.csData(4));
                ail   = sprintf('%-8s', obj.csId);
                f     = cat(1,f,sprintf('        %s%s%s%s%s',csDum,csEnd,csBeg,csNx,ail));
            end
            
        end
        
        function f = genNastranDeck(obj)
            error('Not Coded yet');
        end
        
        function [h,g] = plot(obj,fig,value)
            
            if nargin < 2
                figure;
                value1 = [0 1 0];
                value2 = [1 1 0];
            elseif nargin < 3
                figure(fig);
                hold on;
                value1 = [0 1 0];
                value2 = [1 1 0];
            else 
                figure
                value1 = value;
                value2 = value;
            end
            
            if nargin == 2
                figure(fig);
                hold on;
            else
                figure;
            end
            
            %Plot the wing outline
            NewOutline = [obj.outline];
            h = patch(NewOutline(:,1:3:end),NewOutline(:,2:3:end),NewOutline(:,3:3:end),value1);
            
            %Plot the control surface outline
            NewControlOutline = [obj.controloutline];
            g = patch(NewControlOutline(:,1:3:end),NewControlOutline(:,2:3:end),NewControlOutline(:,3:3:end),value2);
            
            axis equal
            
        end
        
        function [mac,mac_LE,mac_AC] = recoverMAC(obj)
            
            % Identify where the unique spanwise locations are
            
            %[panels,ia,ic] = unique([obj.startY],'stable');
            [~,ia,ic] = unique([obj.startY]+[obj.b],'stable');
            
            np    = numel(ia);
            MACP  = zeros(np,1);
            AREAP = zeros(np,1);
            XLE   = zeros(np,1);
            XAC   = zeros(np,1);
            
            for i = 1:np
                
                Tempobj = (obj(ic == i)); 
                [MACP(i), AREAP(i), XLE(i), XAC(i)] = MAC_patch([Tempobj.outline]);
            end
            
            
            % area weighted results
            AREATOT = sum(abs(AREAP));
            mac     = dot(MACP, abs(AREAP))/ AREATOT;
            mac_LE  = dot(XLE, abs(AREAP)) / AREATOT;
            mac_AC  = dot(XAC, abs(AREAP)) / AREATOT;
            
%             fprintf('\nMean Aerodynamic Chord: %8.4f',mac);
%             fprintf('\nMean Aerodynamic Chord: %8.4f',mac_LE);
%             fprintf('\nMean Aerodynamic Chord: %8.4f',mac_AC);
        end
        
    end
end

function [MAC, A, xLE, xAC] = MAC_patch(TempOutline)

Outline = TempOutline(1:4,1:3);
Outline(1,1) = min(TempOutline(1,1:3:end));
Outline(2,1) = min(TempOutline(2,1:3:end));
Outline(3,1) = max(TempOutline(3,1:3:end));
Outline(4,1) = max(TempOutline(4,1:3:end));

% Take into account that more than one panel maybe input here!
c = Outline(4,1)-Outline(1,1);

tp = (Outline(3,1)-Outline(2,1))/c;

sp = abs(Outline(2,2)-Outline(1,2));

sw = 180*atan(((Outline(2,1) + 0.25*c*tp)-(Outline(1,1) + 0.25*c))/sp)/pi;

% patch MAC
MAC = (2/3)*(1+tp+tp^2).*c/(1+tp);

% patch Area
A = c*(1+tp)*sp/2;

try
    eta = interp1([c,tp*c],[0,1],MAC);
catch
    eta = 0.5;
end

xAC = tan(pi*sw/180).*eta.*sp + 0.25.*c + Outline(1,1);
xLE = xAC - 0.25*MAC;

end

