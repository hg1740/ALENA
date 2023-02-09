classdef PanelPropCollection < awi.fe.FECollection
    %PanelPropCollection Describes a collection of 'awi.fe.PanelProp'
    %objects.
    %
    % See also awi.fe.PanelProp
    
    %Primary properties
    properties
        %Property identification number
        PID
        %Material identification number for the membrane
        MID1
        %Material identification number for the bending
        MID2
        %Material identification number for transverse shear
        MID3
        %Material identification number for membrane-bending coupling
        MID4
        %Membrane thickness
        Thickness
        %Bending moment of inertia ratio (12 * I / T^3)
        BendingInertiaRatio = 1;
        %Ratio of the transverse shear thickness to the membrane thickness
        TransverseShearThicknessRatio = 0.83;
        %Non-structural mass per unit area
        NSM
        %Fibre distances for stress calculations
        FibreDistance
    end
    
    %Handle to 'awi.fe' objects
    properties
        %Handle to the 'awi.fe.Material' object that describes the membrane
        %material properties
        MembraneMaterial
        %Handle to the 'awi.fe.Material' object that describes the bending
        %material properties
        BendingMaterial
        %Handle to the 'awi.fe.Material' object that describes the shear
        %material properties
        ShearMaterial
        %Handle to the 'awi.fe.Material' object that describes the 
        %membrane-bending coupling material properties
        MBCouplingMaterial
    end
    
    methods % set / get
        function set.PID(obj, val)                           %set.PID
            obj.IDNumbers = val;
        end
        function set.MID1(obj, val)                          %set.MID1
            validateIDvector(obj, val, 'MID1');
            obj.MID = val;
        end
        function set.MID2(obj, val)                          %set.MID2
            validateIDvector(obj, val, 'MID2');
            obj.MID = val;
        end
        function set.MID3(obj, val)                          %set.MID3
            validateIDvector(obj, val, 'MID3');
            obj.MID = val;
        end
        function set.MID4(obj, val)                          %set.MID4
            validateIDvector(obj, val, 'MID4');
            obj.MID = val;
        end
        function set.Thickness(obj, val)                     %set.Thickness
            validateattributes(val, {'numeric'}, {'row', 'positive', ...
                'nonnan', 'finite', 'real'}, class(obj), 'Thickness');
            obj.Thickness = val;
        end
        function set.BendingInertiaRatio(obj, val)           %set.BendingInertiaRatio
            validateattributes(val, {'numeric'}, {'row', 'positive', ...
                'nonnan', 'finite', 'real'}, class(obj), 'BendingInertiaRatio');
            obj.BendingInertiaRatio = val;
        end
        function set.TransverseShearThicknessRatio(obj, val) %set.TransverseShearThicknessRatio
            validateattributes(val, {'numeric'}, {'row', 'positive', ...
                'nonnan', 'finite', 'real'}, class(obj), 'TransverseShearThicknessRatio');
            obj.TransverseShearThicknessRatio = val;
        end
        function set.NSM(obj, val)                           %set.NSM
            validateattributes(val, {'numeric'}, {'row', 'positive', ...
                'nonnan', 'finite', 'real'}, class(obj), 'NSM');
            obj.NSM = val;
        end
        function set.FibreDistance(obj, val)                 %set.FibreDistance
            validateattributes(val, {'numeric'}, {'2d', 'nrows', 2, ...
                'nonnan', 'finite', 'real'}, class(obj), 'FibreDistance');
            obj.FibreDistance = val;
        end
        function set.MembraneMaterial(obj, val)              %set.MembraneMaterial
            valdiateattributes(val, {'awi.fe.Material'}, {'scalar'}, ...
                class(obj), 'MembraneMaterial');
            obj.MembraneMaterial = val;
        end
        function set.BendingMaterial(obj, val)               %set.MembraneMaterial
            valdiateattributes(val, {'awi.fe.Material'}, {'scalar'}, ...
                class(obj), 'BendingMaterial');
            obj.BendingMaterial = val;
        end
        function set.ShearMaterial(obj, val)                 %set.ShearMaterial
            valdiateattributes(val, {'awi.fe.Material'}, {'scalar'}, ...
                class(obj), 'ShearMaterial');
            obj.ShearMaterial = val;
        end
        function set.MBCouplingMaterial(obj, val)            %set.MBCouplingMaterial
            valdiateattributes(val, {'awi.fe.Material'}, {'scalar'}, ...
                class(obj), 'MBCouplingMaterial');
            obj.MBCouplingMaterial = val;
        end
        function val = get.MID1(obj)                         %get.MID1
            if isempty(obj.MembraneMaterial)
                val = obj.MID1;
            else
                val = obj.MembraneMaterial.ID;
            end
        end
        function val = get.MID2(obj)                         %get.MID2
            if isempty(obj.BendingMaterial)
                val = obj.MID2;
            else
                val = obj.BendingMaterial.ID;
            end
        end
        function val = get.MID3(obj)                         %get.MID3
            if isempty(obj.ShearMaterial)
                val = obj.MID3;
            else
                val = obj.ShearMaterial.ID;
            end
        end
        function val = get.MID4(obj)                         %get.MID4
            if isempty(obj.MBCouplingMaterial)
                val = obj.MID4;
            else
                val = obj.MBCouplingMaterial.ID;
            end
        end
    end
    
    methods % construction
        function obj = PanelPropCollection
            
            %Make a note of the property names
            addFEProp(obj, 'PID', 'MID1',  'MID2',  'MID3','MID4', ...
                'Thickness', 'BendingInertaRatio', ...
                'TransverseShearThicknessRatio', 'NSM', 'FibreDistance');
            
        end
    end
    
end