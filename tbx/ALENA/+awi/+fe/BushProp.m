classdef BushProp < awi.fe.FEBaseClass
    %BushProp Describes the stiffness properties of a generic bush element
    %for use in a finite element model.
    %
    % The definition of the 'BushProp' object mathces that of the PBUSH
    % bulk data type from MSC.Nastran.
    
    %Primary properties
    properties
        %Property identification number
        PID
        %Stiffness in direction 1
        K1
        %Stiffness in direction 2
        K2
        %Stiffness in direction 3
        K3
        %Stiffness in direction 4
        K4
        %Stiffness in direction 5
        K5
        %Stiffness in direction 6
        K6
    end    
    
    methods % set / get
        function set.PID(obj, val)  %set.PID 
            %set.PID Set method for the property 'PID'.
            %
            % Passes the value straight to the inherited 'ID' property
            % which will validate and store the ID number.
            obj.ID = val;
        end
        function set.K1(obj, val)   %set.K1  
            val(isnan(val)) = 0;    %Replace any NaN terms with zeros
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty', ...
                'finite', 'real', 'nonnan'}, class(obj), 'K1');
            obj.K1 = val;
        end
        function set.K2(obj, val)   %set.K2  
            val(isnan(val)) = 0;    %Replace any NaN terms with zeros
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty', ...
                'finite', 'real', 'nonnan'}, class(obj), 'K2'); 
            obj.K2 = val;
        end
        function set.K3(obj, val)   %set.K3  
            val(isnan(val)) = 0;    %Replace any NaN terms with zeros
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty', ...
                'finite', 'real', 'nonnan'}, class(obj), 'K3'); 
            obj.K3 = val;
        end
        function set.K4(obj, val)   %set.K4  
            val(isnan(val)) = 0;    %Replace any NaN terms with zeros
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty', ...
                'finite', 'real', 'nonnan'}, class(obj), 'K4'); 
            obj.K4 = val;
        end
        function set.K5(obj, val)   %set.K5  
            val(isnan(val)) = 0;    %Replace any NaN terms with zeros
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty', ...
                'finite', 'real', 'nonnan'}, class(obj), 'K5'); 
            obj.K5 = val;
        end
        function set.K6(obj, val)   %set.K6  
            val(isnan(val)) = 0;    %Replace any NaN terms with zeros
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty', ...
                'finite', 'real', 'nonnan'}, class(obj), 'K6'); 
            obj.K6 = val;
        end
    end
    
    methods % construction
        function obj = BushProp
            
            %Make a note of the property names
            addFEProp(obj, 'PID', 'K1', 'K2', 'K3', 'K4', 'K5', 'K6');
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.BushProp' object
            %into a text file using the format of the MSC.Nastran 'PBUSH'
            %bulk data entry.
            %
            % The following assumptions are made:
            %   * The bush element only has stiffness in all six-DOF. No
            %   damping, mass or strain coefficients are defined.
            
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
                comment = ['PBUSH : Defines the nominal property ', ...
                    'values for a generalized spring-and-damper ' , ...
                    'structural element.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'large');
            
            nObj = numel(obj);
            
            %Card name
            nam   = repmat({'PBUSH*'}  , [1, nObj]);
            blnks = repmat({['*', blanks(7)]}, [1, nObj]);
            str   = repmat({'K'}, [1, nObj]);
            
            %Set up the format for printing
            data = [ ... 
                nam   ; {obj.ID} ; str      ; {obj.K1} ; {obj.K2} ; ...
                blnks ; {obj.K3} ; {obj.K4} ; {obj.K5} ; {obj.K6}];
            
            %Write in 16-character column width as standard
            format = [ ...
                '%-8s%-16i%-16s%#-16.5g%#-16.5g\r\n', ...
                '%-8s%#-16.5g%#-16.5g%#-16.5g%#-16.5g\r\n'];
            
            %Write the data to the file
            fprintf(fid, format, data{:});
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

