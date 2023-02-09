classdef ScalarMassProp < awi.fe.ScalarProp
    %ScalarMassProp Describes the mass of a scalar mass element.
    %
    % The defintion of the 'ScalarMassProp' object matches that of the
    % PMASS bulk data type from MSC.Nastran.
    
    %Primary properties
    properties
        %Mass (force per unit acceleration)
        M
    end
    
    methods % set / get
        function set.M(obj, val)    %set.M   
            obj.PropVal = val;
        end       
        function val = get.M(obj)   %get.M
            val = obj.PropVal;
        end
    end
    
    methods % construction
        function obj = ScalarMassProp
            
            addFEProp(obj, 'PID', 'M');
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.ScalarMassProp' object
            %into a text file using the format of the MSC.Nastran PDAMP
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
                comment = ['PMASS : Specifies the mass value of a ', ...
                    'scalar mass element'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, '16');
            
            nObj = numel(obj);
            
            %Card name
            nam   = repmat({'PMASS*'}, [1, nObj]);
            
            %Set up the format for printing
            data = [nam ; {obj.ID} ; {obj.M}];
            
            %Write in 16-character column width as standard
            format = '%-8s%-16i%#-16.6g\r\n';
            
            %Write the data to the file
            fprintf(fid, format, data{:});
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

