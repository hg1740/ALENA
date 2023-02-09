classdef ElementResults < matlab.mixin.SetGet
    %ElementResults Defines a set of stress results quantities 
    
    %Results data
    properties
        %Handle to the 'awi.fe.Beam' objects -> nNode = numel(obj.Nodes)
        Elements    
        %Time vector - > nT = numel(obj.TimeVector)
        TimeVector  
        
        % Results from "H5Data.ResultSets.STRESS.BEAM" (stress data are in a 1 x nElements x 11 array)
        EID % Element ID
        GRID % Grid ID
        SD % Identification number of coordinate system
        XC % Stress at recovery point C
        XD % Stress at recovery point D
        XE % Stress at recovery point E
        XF % Stress at recovery point F
        MAX % Maximum stress from all recovery points
        MIN % Minimum stress from all recovery points
        
        % Results from "H5Data.ResultSets.STRESS.BEAM_NL"
        % NSX: Stress
        % NSE: Equivalent stress
        % TE: Total strain
        % EPE: Eff. strain Plastic/NLElastic
        % ECE: Eff. creep strain
        GRIDA % Grid ID at end A
        NSXCA
        NSECA
        TECA
        EPECA
        ECECA
        NSXDA
        NSEDA
        TEDA
        EPEDA
        ECEDA
        NSXEA
        NSEEA
        TEEA
        EPEEA
        ECEEA
        NSXFA
        NSEFA
        TEFA
        EPEFA
        ECEFA
        GRIDB % Grid ID at end B
        NSXCB
        NSECB
        TECB
        EPECB
        ECECB
        NSXDB
        NSEDB
        TEDB
        EPEDB
        ECEDB
        NSXEB
        NSEEB
        TEEB
        EPEEB
        ECEEB
        NSXFB
        NSEFB
        TEFB
        EPEFB
        ECEFB
        
    end
    
    %Descriptors
    properties
        %Results type
        ResultsType
        %Title describing the results SET
        Title
    end
    
    properties (Constant)
        ValidResultsType = {'Stress'};
    end
    
    methods % set / get
        function set.Elements(obj, val)        %set.Nodes
            validateattributes(val, {'awi.fe.Beam'}, ...
                {'row'}, class(obj), 'Elements');
            obj.Elements = val;
        end
        function set.TimeVector(obj, val)   %set.TimeVector
            validateattributes(val, {'numeric'}, {'row', 'nonnan', ...
                'finite', 'real', 'nonnegative'}, class(obj), 'TimeVector');
            obj.TimeVector = val;
        end
        function set.ResultsType(obj, val)  %set.ResultsTypes
            val = validatestring(val, obj.ValidResultsType, class(obj), 'ValidResultsTypes');
            obj.ResultsType = val;
        end
        function set.Title(obj, val)        %set.Title
            validateattributes(val, {'char'}, {'row'}, class(obj), 'Title');
            obj.Title = val;
        end
        function set.EID(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'EID');
            obj.EID = val;
        end
        function set.GRID(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'nrows', 1, ...
                'finite'}, class(obj), 'GRID');
            obj.GRID = val;
        end
        function set.SD(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'nrows', 1, ...
                'finite'}, class(obj), 'SD');
            obj.SD = val;
        end
        function set.XC(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'nrows', 1, ...
                'finite'}, class(obj), 'XC');
            obj.XC = val;
        end
        function set.XD(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'nrows', 1, ...
                'finite'}, class(obj), 'XD');
            obj.XD = val;
        end
        function set.XE(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'nrows', 1, ...
                'finite'}, class(obj), 'XE');
            obj.XE = val;
        end
        function set.XF(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'nrows', 1, ...
                'finite'}, class(obj), 'XF');
            obj.XF = val;
        end
        function set.MAX(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'nrows', 1, ...
                'finite'}, class(obj), 'MAX');
            obj.MAX = val;
        end
        function set.MIN(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'nrows', 1, ...
                'finite'}, class(obj), 'MIN');
            obj.MIN = val;
        end
        function set.GRIDA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'GRIDA');
            obj.GRIDA = val;
        end
        
        function set.NSXCA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'NSXCA');
            obj.NSXCA = val;
        end
        function set.NSECA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'NSECA');
            obj.NSECA = val;
        end
        function set.TECA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'TECA');
            obj.TECA = val;
        end
        function set.EPECA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'EPECA');
            obj.EPECA = val;
        end
        function set.ECECA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'ECECA');
            obj.ECECA = val;
        end
        
        function set.NSXDA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'NSXDA');
            obj.NSXDA = val;
        end
        function set.NSEDA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'NSEDA');
            obj.NSEDA = val;
        end
        function set.TEDA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'TEDA');
            obj.TEDA = val;
        end
        function set.EPEDA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'EPEDA');
            obj.EPEDA = val;
        end
        function set.ECEDA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'ECEDA');
            obj.ECEDA = val;
        end
        
        function set.NSXEA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'NSXEA');
            obj.NSXEA = val;
        end
        function set.NSEEA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'NSEEA');
            obj.NSEEA = val;
        end
        function set.TEEA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'TEEA');
            obj.TEEA = val;
        end
        function set.EPEEA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'EPEEA');
            obj.EPEEA = val;
        end
        function set.ECEEA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'ECEEA');
            obj.ECEEA = val;
        end
        
        function set.NSXFA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'NSXFA');
            obj.NSXFA = val;
        end
        function set.NSEFA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'NSEFA');
            obj.NSEFA = val;
        end
        function set.TEFA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'TEFA');
            obj.TEFA = val;
        end
        function set.EPEFA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'EPEFA');
            obj.EPEFA = val;
        end
        function set.ECEFA(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'ECEFA');
            obj.ECEFA = val;
        end
        
      	function set.GRIDB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'GRIDA');
            obj.GRIDB = val;
        end
        
        function set.NSXCB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'NSXCB');
            obj.NSXCB = val;
        end
        function set.NSECB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'NSECB');
            obj.NSECB = val;
        end
        function set.TECB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'TECB');
            obj.TECB = val;
        end
        function set.EPECB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'EPECB');
            obj.EPECB = val;
        end
        function set.ECECB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'ECECB');
            obj.ECECB = val;
        end
        
        function set.NSXDB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'NSXDB');
            obj.NSXDB = val;
        end
        function set.NSEDB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'NSEDB');
            obj.NSEDB = val;
        end
        function set.TEDB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'TEDB');
            obj.TEDB = val;
        end
        function set.EPEDB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'EPEDB');
            obj.EPEDB = val;
        end
        function set.ECEDB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'ECEDB');
            obj.ECEDB = val;
        end
        
        function set.NSXEB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'NSXEB');
            obj.NSXEB = val;
        end
        function set.NSEEB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'NSEEB');
            obj.NSEEB = val;
        end
        function set.TEEB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'TEEB');
            obj.TEEB = val;
        end
        function set.EPEEB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'EPEEB');
            obj.EPEEB = val;
        end
        function set.ECEEB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'ECEEB');
            obj.ECEEB = val;
        end
        
        function set.NSXFB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'NSXFB');
            obj.NSXFB = val;
        end
        function set.NSEFB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'NSEFB');
            obj.NSEFB = val;
        end
        function set.TEFB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'TEFB');
            obj.TEFB = val;
        end
        function set.EPEFB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'EPEFB');
            obj.EPEFB = val;
        end
        function set.ECEFB(obj, val)
            validateattributes(val, {'numeric'}, {'3d', 'column', ...
                'finite'}, class(obj), 'ECEFB');
            obj.ECEFB = val;
        end
    end
    
end

