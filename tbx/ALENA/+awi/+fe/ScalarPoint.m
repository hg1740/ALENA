classdef ScalarPoint < awi.fe.Node
    %ScalarPoint Describes a point in the mathermatical formulation of the
    %FE model which has no physical location in 3D-space.
    %
    % The definition of the 'ScalarPoint' object matches that of the SPOINT
    % bulk data type from MSC.Nastran.
    %
    % For visualisation purposes a scalar point can be given a point in 3D
    % space however this information is not transferred to the finite
    % element formulation.
        
    methods
        function obj = ScalarPoint
           
            %Make a note of the property names
            addFEProp(obj, 'ID');
            
        end
    end
    
    methods 
        function writeToFile(obj, fid, bComment)
            
            for i = 1 : numel(obj)
                fprintf(fid, '%-8s%-8i\r\n', 'SPOINT', obj(i).ID);
            end            
            
        end
    end
    
end

