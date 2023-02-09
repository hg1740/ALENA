classdef (ConstructOnLoad) Dimensionable < matlab.mixin.SetGet
    %
    % Provides capability for handling dimensions associated with physical units
    
    properties (SetObservable)
        
        %Physical units associated with this instance ?
        Unit = '';
        
    end
    
    properties (Access = protected)
        
        %Which field of superclass does this physical unit apply ?
        ValueField = 'Value';
        
    end
    
    methods
        
        function [val, unit] = getValue(obj, requiredUnit)
            
            %Start with basic value and unit
            val = obj.(obj.ValueField);
            
            %Polynomial to reduce the value to MKS
            [p, mks_p] = i_tomks(obj.Unit);
            
            %Is that it ?
            if nargin == 1
                
                %Apply conversion
                val = polyval(p, val);
                
                %Send back MKS units
                unit = mks_p;
                
            else
                
                %Get conversion associated with required unit
                [q, mks_q] = i_tomks(requiredUnit);
                
                %Ensure consistency
                assert(strcmp(mks_p, mks_q), 'required unit not consistent with actual');
                
                %Combine conversion polynomials
                r = i_divide(p, q);
                
                %Apply conversion
                val = polyval(r, val);
                
                %Send back the units we now have
                unit = requiredUnit;
                
            end
            
        end
        
        
    end
    
end

function [p, mks] = i_tomks(unit)
%
% Returns the polynomial coefficients P that reduce UNIT to MKS

%Tokenise the unit
tok = strsplit(unit, '/');

%Enlist help from SimScape
[units, conversions, expressions] = pm_getunits;

%Look for match with numerator
b = strcmp(tok{1}, units);
assert(any(b), ['no match with numerator ''', tok{1}, '''']);

%Make a note of conversion polynomial coefficients
p = conversions(b,:);

%Make a note of MKS unit
mks{1} = expressions{b};

%For each denominator
for i = 2:numel(tok)
    
    %Look for match
    b = strcmp(tok{i}, units);
    assert(any(b), ['no match with denominator ''', tok{i}, '''']);
    
    %Combine conversions by polynomial division
    p = i_divide(p, conversions(b,:));
    
    %Make a note of MKS unit
    unit{i} = expressions{b};
    
end

%So the physical unit associated with return value is
mks = strjoin(mks, '/');

end

function r = i_divide(p, q)

%Polynomial division == deconvolution
[r, rem] = deconv(p, q);

%Followed by
r = [r, 0] + rem;

end
