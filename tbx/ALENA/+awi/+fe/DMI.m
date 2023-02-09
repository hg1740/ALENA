classdef DMI < awi.fe.FEBaseClass
    %DMI describe a user defined matrix submitted as part of the bulk data
    % 
    % The definition of the 'DMI' object matches that of the DMI bulk 
    % data type from MSC.Nastran.   
    
    %Primary Properties
    properties
        % Name of the Matrix
        NAME
        % form of the matrix
        FORM = 2
        % Input Matrix type
        TIN = 1
        % Output Matrix type
        TOUT = 0
        % Matrix Data
        DATA
    end
    
    methods
        function obj = DMI()
            
        end
    end
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            %By default, do not close the file
            bClose = false;
            if isempty(obj.DATA)
                return
            end                        
            if nargin < 2 %Ask the user for the file
                fName  = awi.fe.FEBaseClass.getBulkDataFile;
                bClose = true;
                fid    = fopen(fName, 'w');                
            end
            if nargin < 3 %Comments by standard
                bComment = true;
            end        
            if bComment
                mni.printing.bdf.writeComment('DMI card', fid);
            end
            mni.printing.bdf.writeColumnDelimiter(fid, 'normal');
            dmi_card = mni.printing.cards.DMI(obj.NAME,obj.DATA,...
                obj.FORM,obj.TIN,obj.TOUT);
            dmi_card.writeToFile(fid);
        end
    end
end

