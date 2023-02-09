classdef (ConstructOnLoad) Material < awi.model.Entity & awi.mixin.BeamPropertyObject
    %Material Describes the properties of a material during the pre- and
    %post-failure states.
    %
    % The parameterisation of the 'awi.model.Material' class is derived
    % directly from the CPACS implementation.
        
    %Material type
    properties
        %Shortcut for controlling some of the material properties
        MaterialType = 'isotropic';
    end
    
    %Density
    properties
        %Density (assumed homogeneous)
        Rho = 0;
    end
    
    %Modulus of elasticity
    properties
        %Modulus of elasticity in direction 1-1 (Young's Modulus - E)
        K11 = 0;
        %Modulus of elasticity in direction 1-2 (Shear Modulus - G)
        K12 = 0;
        %Modulus of elasticity in direction 1-3
        K13 = 0;
        %Modulus of elasticity in direction 1-4
        K14 = 0;
        %Modulus of elasticity in direction 1-5
        K15 = 0;
        %Modulus of elasticity in direction 1-6
        K16 = 0;
        %Modulus of elasticity in direction 2-1
        K21 = 0;
        %Modulus of elasticity in direction 2-2
        K22 = 0;
        %Modulus of elasticity in direction 2-3
        K23 = 0;
        %Modulus of elasticity in direction 2-4
        K24 = 0;
        %Modulus of elasticity in direction 2-5
        K25 = 0;
        %Modulus of elasticity in direction 2-6
        K26 = 0;
        %Modulus of elasticity in direction 3-1
        K31 = 0;
        %Modulus of elasticity in direction 3-2
        K32 = 0;
        %Modulus of elasticity in direction 3-3
        K33 = 0;   
        %Modulus of elasticity in direction 3-4
        K34 = 0; 
        %Modulus of elasticity in direction 3-5
        K35 = 0; 
        %Modulus of elasticity in direction 3-6
        K36 = 0;         
    end
    
    %Poisson's Ratio
    properties
        Nu12 = 0;
        Nu21 = 0;
    end
    
    %Direct failure stresses (tensile and compressive)
    properties
        %Yield stress
        SigmaYield = -1;
        %Tensile failure stress in direction 1-1
        Sigma11_T  = -1;
        %Compressive failure stress in direction 1-1
        Sigma11_C  = -1;
        %Tensile failure stress in direction 2-2
        Sigma22_T  = -1;
        %Compressive failure stress in direction 2-2
        Sigma22_C  = -1;
    end
    
    %Shear failure stresses
    properties
        %Shear failure stress in the plane 1-2
        Tau12 = -1;
        %Shear failure stress in the plane 2-3
        Tau23 = -1;
    end
    
    %Post-failure behaviour
    properties
        PlaticEliminationStrain
        TangentModulus
        TrueStress
    end
        
    %K-Matrix
    properties (Dependent)
        %Modulus matrix
        K
    end
    
    %Shortcuts
    properties
        %Young's Modulus
        E
        %Shear Modulus
        G
        %Poisson's Ratio
        Nu
    end
    
    methods % set / get
        function set.MaterialType(obj, val) %set.MaterialType 
            validatestring(val, {'isotropic'}, class(obj), 'MaterialType');
            obj.MaterialType = val;
        end      
        function set.E(obj, val)            %set.E  
           %set.E Set method for the property 'E'.
           
           obj.K11 = val; %#ok<*MCSUP>
        end
        function set.G(obj, val)            %set.G  
           %set.g Set method for the property 'g'.
           
           obj.K12 = val; %#ok<*MCSUP>
        end
        function set.Nu(obj, val)           %set.Nu 
            %set.Nu Set method for the property 'Nu'.
            
            obj.Nu12 = val;
        end
        function val = get.K(obj)           %get.K  
            
            %K-matrix is dependent on material type
            switch obj.MaterialType
                case 'isotropic'
                    %Same modulus everywhere!
                    val = repmat(obj.K11, [6, 6]);
                otherwise
                    val = [];
            end
            
            %Construct the 6x6 matrix
            %val = [];
%             val = [ ...
%                 obj.K11, obj.K12, obj.K13, obj.K14, obj.K15, obj.K16 ; ...
%                 obj.K21, obj.K22, obj.K23, obj.K24, obj.K25, obj.K26 ; ...
%                 obj.K31, obj.K32, obj.K33, obj.K34, obj.K35, obj.K36 ; ...
%                 obj.K41, obj.K42, obj.K43, obj.K44, obj.K45, obj.K46 ; ...
%                 obj.K51, obj.K52, obj.K53, obj.K54, obj.K55, obj.K56 ; ...
%                 obj.K61, obj.K62, obj.K63, obj.K64, obj.K65, obj.K66];
        end
        function val = get.E(obj)           %get.E  
            %get.E Get method for the dependent property 'E'.
            
            val = obj.K11;
        end
        function val = get.G(obj)           %get.G  
            %get.G Get method for the dependent property 'G'.
            
            val = obj.K12;
        end
        function val = get.Nu(obj)          %get.Nu 
            %get.Nu Get method for the dependent property 'Nu'.
            
            val = obj.Nu12;
        end
    end
    
end

