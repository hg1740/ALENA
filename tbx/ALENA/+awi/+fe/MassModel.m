classdef MassModel < awi.fe.FEModel
    %MassModel Describes a collection of mass.
    %
    %
    % Much the same as a standard 'awi.fe.FEModel' object except when the
    % object is written to a file all of the 'awi.fe.PointMass' objects are
    % collected into a single MSC.Nastran 'MASSSET' entry.
    
    %ID number
    properties 
        %Initial counter for the Mass Case ID numbers. No MassModelwill be
        %assigned an ID number lower than this when using the method
        %'assignIDnumbers' in the 'awi.fe.FEModel' object.
        MassID0 = 900;  
    end
    
    methods % set / get
        function set.MassID0(obj, val) %set.MassID0
            %set.MassID0 Set method for the property 'MassID0'.
            
            validateID(obj, val, 'ID0');
            assert(val < 999, ['Cannot assign a Mass ID number ', ...
                'greater than 999.']);
            obj.MassID0 = val;
            
        end
    end
    
    methods % exporting/writing FE data to a file
        function fNam = writeToFile(obj, dn, bComment, bCollect, bWriteOpt, bDesModel)
            %writeToFile Writes the model to a file in the folder specified
            %by 'dn'.
            
            %If we have more than one 'MassCase' rename the first one as
            %that is the name that will be taken for the file.
            if numel(obj) > 1
                nam = obj(1).Name;
                obj(1).Name = 'MassCases';
            end
            
            %Pass it on - Write the mass ID token!
            bMassCase = true;
            fNam = writeToFile@awi.fe.FEModel(obj, dn, bComment, bCollect, bWriteOpt, bDesModel, bMassCase);
                        
            %Return original name to the MassCase
            if numel(obj) > 1
                obj(1).Name = nam;
            end
            
        end        
    end
    
end

