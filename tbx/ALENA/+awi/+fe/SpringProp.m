classdef SpringProp < awi.fe.ScalarProp
    %SpringProp Describes the stiffness, damping coefficient and stress
    %coefficient of a generic spring element.
    %
    % The defintion of the 'SpringProp' object matches that of the PELAS
    % bulk data type from MSC.Nastran.
    
    %Primary properties
    properties
        %Spring stiffness (force per unit displacement)
        K
        %Damping coefficient
        GE
        %Stress coefficient
        S
    end

    methods % set / get
        function set.K(obj, val)    %set.K   
            %Pass it on
            obj.PropVal = val;
        end
        function set.GE(obj, val)   %set.GE  
            validateattributes(val, {'numeric'}, {'scalar', 'finite', ...
                'real', 'nonnan'}, class(obj), 'GE');
            obj.GE = val;
        end
        function set.S(obj, val)    %set.S   
            validateattributes(val, {'numeric'}, {'scalar', 'finite', ...
                'real', 'nonnan'}, class(obj), 'S');
            obj.S = val;
        end
        function val = get.K(obj)   %get.K   
           val = obj.PropVal; 
        end
    end
    
    methods % construction
        function obj = SpringProp
            
            addFEProp(obj, 'PID', 'K', 'GE', 'S');
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.SpringProp' object
            %into a text file using the format of the MSC.Nastran PELAS
            %bulk data entry.
                   
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
                comment = ['PELAS : Specifies the stiffness, damping ', ...
                    'coefficient, and stress coefficient of a scalar ', ...
                    'elastic (spring) element'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, '16');
            
            nObj = numel(obj);
            
            %Card name
            nam   = repmat({'PELAS*'}, [1, nObj]);
            
            %Set up the format for printing
            data = [nam ; {obj.ID} ; {obj.K} ; {obj.GE} ; {obj.S}];
            
            %Write in 16-character column width as standard
            format = '%-8s%-16i%#-16.6g%#-16.6g%#-16.6g\r\n';
            
            %Write the data to the file
            fprintf(fid, format, data{:});
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

