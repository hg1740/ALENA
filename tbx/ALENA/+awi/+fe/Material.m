classdef Material < awi.fe.FEBaseClass
    %Material Describes an istropic material for use in a finite element
    %model.
    %
    % The definition of the 'Material' object matches that of the MAT1 bulk
    % data type from MSC.Nastran.
    
    %Primary Properties
    properties
        %ID number of the material 
        MID
        %Young's Modulus
        E
        %Shear Modulus
        G
        %Poissons Ratio
        Nu
        %Material Density
        Rho = 0;
    end
    
    methods % set / get
        function set.MID(obj, val) %set.MID 
            %set.MID Set method for the property 'MID'.
            % 
            % Passes the value straight to the inherited 'ID' property
            % which will validate and store the ID number.
            obj.ID = val;
        end
        function set.E(obj, val)   %set.E   
            validateattributes(val, {'numeric'}, {'scalar', 'nonnegative', ...
                'finite', 'real'}, class(obj), 'E');
            obj.E = val;
        end
        function set.G(obj, val)   %set.G   
            validateattributes(val, {'numeric'}, {'scalar', 'nonnegative', ...
                'finite', 'real'}, class(obj), 'G');
            obj.G = val;
        end
        function set.Nu(obj, val)  %set.Nu  
            if isempty(val)
                obj.Nu = val;
                return
            end
            validateattributes(val, {'numeric'}, {'scalar', 'nonnegative', ...
                'finite', 'real'}, class(obj), 'Nu');
            obj.Nu = val;
        end
        function set.Rho(obj, val) %set.Rho 
           validateattributes(val, {'numeric'}, {'scalar', 'nonnegative', ...
                'finite', 'real'}, class(obj), 'Rho');
            obj.Rho = val;
        end
    end
    
    methods % construction
        function obj = Material
            
            %Make a note of the property names
            addFEProp(obj, 'MID', 'E', 'G', 'Nu', 'Rho');
            
        end
    end
 
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.Material' object
            %into a text file using the format of the MSC.Nastran 'MAT1'
            %bulk data entry. 
            %
            % The following assumptions are made:
            %   * The thermal expansion coefficient (A) is assumed to be 0.
            %   * The reference temperature is assumed to be 0.
            %   * The structural damping coefficient for the material is
            %   assumed to be 0.
            %   * The stress limits are not specified.
            %   * The material coordinate system identification (MCSID)
            %   number is not specified.
            
            %By default, do not close the file
            bClose = false;
                        
            if nargin < 2 %Ask the user for the file
                fName = awi.fe.FEBaseClass.getBulkDataFile;
                bClose = true;
                fid = fopen(fName, 'w');                
            end
            
            if nargin < 3 %Comments by standard
                bComment = true;
            end
                        
            if bComment %Helpful comments?
                comment = ['MAT1 : Defines the material properties ', ...
                    'for linear isotropic materials.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'large');
                 
            nObj = numel(obj);
            
            %Card name
            nam   = repmat({'MAT1*'}  , [1, nObj]);
            blnks = repmat({['*', blanks(7)]}, [1, nObj]);
            
            %Set up the format for printing
            data = [ ...
                nam   ; {obj.ID} ; {obj.E} ; {obj.G} ; {obj.Nu} ; ...
                blnks ; {obj.Rho}];
            
            %Write in 16-character column width as standard
            format = [ ...
                '%-8s%-16i%#-16.8g%#-16.8g%#-16.8g\r\n', ...
                '%-8s%#-16.8g\r\n'];
            
            %Write the data to the file
            fprintf(fid, format, data{:});
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

