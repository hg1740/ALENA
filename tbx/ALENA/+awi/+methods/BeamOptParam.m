classdef BeamOptParam < matlab.mixin.SetGet
    %BeamOptParam Data class describing the inputs to the beam structural
    %optimisation.
    
    %Structural
    properties
        %Box beam idealisation
        BoxBeamType = 'SymmetricBox';
        %Minimum thickness of box components [m]
        MinThick = 2.5e-3;
        %Maximum thickness of box components [m]
        MaxThick = 5e-2;
        %Initial thickness of box components [m]
        InitialThick = 2.5e-2;
        %Method for calculating box-beam height
        SparHeight = 'average';
        %Portion of out-of-plane load carried by the spar/cover elements
        TransverseLoadFactor = 0.5;
        %Stringer-to-skin area ratio
        As_bT = 0.5;
        %Rib Pitch [m]
        RibPitch = 0.9;
        %Stringer pitch = [m]
        StringerPitch = 0.15;
    end
    
    %Sizing
    properties (Constant = true)
        %Type of sizing study {'optimisation', 'sizing'}
        Method = 'optimisation';
        %Determines whether to use Nastran's built-in cross-sectional beam
        %library or use design equations to generate custom cross-sections
        UseBeamCrossSection = false;
        %Determines whether design variables are generated for the
        %port-side objects
        OptXZSymmetry = true;
        %Determines whether any port-side models are retained in the model
        ModelXZSymmetry = false;
    end
    
    %Constraints
    properties
        %Flag to apply strength constraints
        StrengthConstr   = true;
        %Flag to apply panel/strigner buckling constraints
        LocalStabConstr  = false;
        %Flag to apply the global buckling constraint
        GlobalBuckConstr = false;
        %Flag to apply analytical (Euler) buckling constraints
        EulerBuckConstr  = false;
        %Flag to apply flutter constraints
        FlutterConstr    = true;
        %Flag to apply divergence constraint
        DivergeConstr    = false;
    end
    
    %Materials
    properties
        %Type of material {'metal', 'composite'}
        MatType = 'metal';
    end
    
    %Reserve Factors
    properties
        %Reserve Factor for strength constraints
        StrengthRF   = 1.5;
        %Reserve Factor for stability constraints
        StabilityRF  = 1;
        %Reserve Factor for global buckling loads
        GlobalBuckRF = 1;
        %Reserve Factor for flutter speeds
        FlutterRF    = 1.15;
        %Maximum Reserve Factor for strength/stability/buckling constraints
        MaxRF        = 1e20;
    end
    
    %Flutter
    properties
        %Number of flutter modes to track in SOL 200
        NumFlutterModes = 10;
        %Flag to ignore rigid body modes when tracking flutter modes
        IgnoreRigidBodyModes = true;
    end
    
    %SOL 200 Parameters
    properties
        %Maximum number of design cycles
        DesMax = 100;
        %Fractional change allowed in each property 
        DelP   = 1.25e-2;
        %Minimum move limit on property values
        DPMin  = 6.25e-4;
        %Fractional change allowed in each design variable
        DelX   = 3.125e-2;
        %Minimum change in design variable
        DXMin  = 3.125e-3
        %Maximum allowable change in design variable
        DXMax  = 1;
        %Relative change in objective function for convergence
        DelObj = 0.0001;
        %Absolute change in objective function for convergence
        DelAbsObj = 10;
        %Logical flag indicating whether to print output of initial
        %analysis in .f06 file
        PrintInitialResults = false;
        %Flag to choose whether to retain all constraints
        ScreenConstr = false;
    end
    
    %Flag for sensitivity analysis 
    properties
        %Logical flag that tells Nastran to exit after sensitivity analysis
        ExitAfterSensitivity = false;
    end
    
    %Helper
    properties (Constant = true, Hidden = true)
        %Property/Parameter mapping
        Sol200ParamMap = { ...
            'Desmax', 'PrintInitialResults', 'DelP', 'DPMin', 'DelX', 'DXMin', 'DXMax', 'DelObj' , 'DelAbsObj'; ...
            'DESMAX', 'NASPR0'             , 'DELP', 'DPMIN', 'DELX', 'DXMIN', 'DXMAX', 'DELOBJ' , 'DABOBJ'};
    end
    
    methods % set / get
        function obj = set.BoxBeamType(obj, val)
            validatestring(val, {'SymmetricBox', 'Box'}, class(obj), 'BoxBeamType');
            obj.BoxBeamType = val;
        end
        function obj = set.MinThick(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'finite', 'real', 'nonnan'}, class(obj), 'MinThick');
            obj.MinThick = val;
        end
        function obj = set.MaxThick(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'finite', 'real', 'nonnan'}, class(obj), 'MaxThick');
            obj.MaxThick = val;
        end
        function obj = set.InitialThick(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'finite', 'real', 'nonnan'}, class(obj), 'InitialThick');
            obj.InitialThick = val;
        end
        function obj = set.SparHeight(obj, val)
            validatestring(val, {'average', 'min', 'max'}, class(obj), 'SparHeight');
            obj.SparHeight = val;
        end
        function obj = set.TransverseLoadFactor(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'nonnegative', ...
                'finite', 'real', 'nonnan', '<=', 1}, class(obj), 'TransverseLoadFactor');
            obj.TransverseLoadFactor = val;
        end
        function obj = set.As_bT(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'nonnegative', ...
                'finite', 'real', 'nonnan', '<=', 1}, class(obj), 'As_bT');
            obj.As_bT = val;
        end
        function obj = set.RibPitch(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'finite', 'real', 'nonnan'}, class(obj), 'RibPitch');
            obj.RibPitch = val;
        end
        function obj = set.StringerPitch(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'finite', 'real', 'nonnan'}, class(obj), 'StringerPitch');
            obj.StringerPitch = val;
        end
        function obj = set.StrengthConstr(obj, val)
            validateattributes(val, {'logical'}, {'scalar'}, class(obj), ...
                'StrengthConstr');
            obj.StrengthConstr = val;
        end
        function obj = set.LocalStabConstr(obj, val)
            validateattributes(val, {'logical'}, {'scalar'}, class(obj), ...
                'LocalStabConstr');
            obj.LocalStabConstr = val;
        end
        function obj = set.GlobalBuckConstr(obj, val)
            validateattributes(val, {'logical'}, {'scalar'}, class(obj), ...
                'GlobalBuckConstr');
            obj.GlobalBuckConstr = val;
        end
        function obj = set.EulerBuckConstr(obj, val)
            validateattributes(val, {'logical'}, {'scalar'}, class(obj), ...
                'EulerBuckConstr');
            obj.EulerBuckConstr = val;
         end        
        function obj = set.FlutterConstr(obj, val)
            validateattributes(val, {'logical'}, {'scalar'}, class(obj), ...
                'FlutterConstr');
            obj.FlutterConstr = val;
        end
        function obj = set.DivergeConstr(obj, val)
            validateattributes(val, {'logical'}, {'scalar'}, class(obj), ...
                'DivergeConstr');
            obj.DivergeConstr = val;
        end
        function obj = set.MatType(obj, val)
            validatestring(val, {'metal', 'composite'}, class(obj), 'MatType');
            obj.MatType = val;
        end
        function obj = set.StrengthRF(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'finite', 'real', 'nonnan'}, class(obj), 'StrengthRF');
            obj.StrengthRF = val;
        end
        function obj = set.StabilityRF(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'finite', 'real', 'nonnan'}, class(obj), 'StabilityRF');
            obj.StabilityRF = val;
        end
        function obj = set.GlobalBuckRF(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'finite', 'real', 'nonnan'}, class(obj), 'GlobalBuckRF');
            obj.GlobalBuckRF = val;
        end
        function obj = set.FlutterRF(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'finite', 'real', 'nonnan'}, class(obj), 'FlutterRF');
            obj.FlutterRF = val;
        end
        function obj = set.MaxRF(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'finite', 'real', 'nonnan'}, class(obj), 'MaxRF');
            obj.MaxRF = val;
        end
        function obj = set.NumFlutterModes(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', ...
                'nonnan', 'finite', 'real'}, class(obj), 'NumFlutterModes');
            obj.NumFlutterModes = val;
        end
        function obj = set.IgnoreRigidBodyModes(obj, val)
            validateattributes(val, {'logical'}, {'scalar'}, class(obj), ...
                'IgnoreRigidBodyModes');
            obj.IgnoreRigidBodyModes = val;
        end
        function obj = set.DesMax(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'integer', 'finite', 'real', 'nonnan'}, class(obj), 'DesMax');
            obj.DesMax = val;
        end
        function obj = set.DelP(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'finite', 'real', 'nonnan', '<=', 1}, class(obj), 'DelP');
            obj.DelP = val;
        end
        function obj = set.DPMin(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'finite', 'real', 'nonnan', '<=', 1}, class(obj), 'DPMin');
            obj.DPMin = val;
        end
        function obj = set.DelX(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'finite', 'real', 'nonnan', '<=', 1}, class(obj), 'DelX');
            obj.DelX = val;
        end
        function obj = set.DXMin(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'finite', 'real', 'nonnan', '<=', 1}, class(obj), 'DXMin');
            obj.DXMin = val;
        end
        function obj = set.DelObj(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'finite', 'real', 'nonnan', '<=', 1}, class(obj), 'DelObj');
            obj.DelObj = val;
        end
        function obj = set.ScreenConstr(obj, val)
            validateattributes(val, {'logical'}, {'scalar'}, class(obj), ...
                'ScreenConstr');
            obj.ScreenConstr = val;
        end
        function obj = set.ExitAfterSensitivity(obj, val)
            validateattributes(val, {'logical'}, {'scalar'}, class(obj), ...
                'ExitAfterSensitivity');
            obj.ExitAfterSensitivity = val;
        end
    end
    
end

