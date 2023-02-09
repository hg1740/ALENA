classdef Aerofoil < awi.model.Component
    %Aerofoil Defines a generic aerofoil cross-section located on a
    %'LiftingSurface' object.
    %
    % An Aerofoil can be defined using a normalised cross-section which is
    % defined in an external file, or they can be defined using the NACA
    % aerofoil series.
    
    
    properties
        %Name of the cross-section.
        Section
        %Non-dimensional position along the parent stick.
        Eta
        %Defines the dimension used for interpolation.
        Eta_Flag = 'R';
    end
    
    methods
        function set.Section(obj, val)
            %set.Name Set method for the property 'Section'.
            %
            % Rules:
            %   - Type : 'char'|'string'
            %   - Attr : vector
            %
            
            validateattributes(val, {'string'}, {'vector', 'nonempty', ...
                'real'}, class(obj), 'Section');
            
            %Preallocate cell array
            C = cell(size(val));
            
            %Extract data from external files
            extID = obj.isExternalFile(val);             
            C(extID) = arrayfun(@(filename) ...
                extractCrossSectionFromFile(obj, filename), val(extID), ...
                'Unif', false);
            
            %Look for special tokens
            %   - Tokens : 'NACA', 'Circle'
            nacaID    = contains(val(~extID), 'NACA');
            C(nacaID) = getNACAProfile(val(nacaID));
        end
    end
    
    methods 
       function data = extractCrossSectionFromFile(obj, filename)
            %extractCrossSectionFromFile
            
            %Do nothing for now
            data = {};
       end 
       
       function data = calculateNacaCoords(obj, NACA, nPoints)
           %calculateNacaCoords
           data = {};
       end
    end
    
    methods (Static)
        function b = isExternalFile(stringOrChar) 
            %isExternalFile Checks if the string input denotes input from 
            %an external file.
            %
            % An external file is denoted by the token '@' as the first
            % character of the string.
    
            validateattributes(stringOrChar, {'char', 'string'}, { ...
                'vector', 'nonempty', 'real'}, 'awi.model.CrossSection');
            
            sz = size(stringOrChar);
            
            %Convert to character array if string is provided.
            if isstring(stringOrChar)
                stringOrChar = char(stringOrChar);
                stringOrChar = permute(stringOrChar, [3, 2, 1]);
            end
            
            %External file is denoted by the token '@'.
            b = stringOrChar(:, 1) == '@';
                
            %Reshape to original size
            b = reshape(b, sz);
                        
        end
    end
end

