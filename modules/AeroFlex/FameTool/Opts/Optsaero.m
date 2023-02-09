%% OPTSAERO Set aerodynamic options
%

%   Copyright 2016 University of Bristol
%   Private function.
classdef Optsaero
    
    properties
        wing_nSpan       = 32;     % Number of spanwise panels for whole wing
        wing_nChord      = 1;      % Number of chordwise panels for whole wing%12
        wing_addTwist    = true;    % Add twist to the Aero Mesh, 1 = YES, 0 = NO.
        wing_addControls = true;    % Add control surfaces, 1 = YES, 0 = NO. 
        htp_nSpan        = 1;       % Number of spanwise panels for whole wing %10
        htp_nChord       = 1;       % Number of chordwise panels for whole wing%12
        htp_addTwist     = true;     % Add twist to the Aero Mesh, 1 = YES, 0 = NO.
        htp_addControls  = true;     % Add control surfaces, 1 = YES, 0 = NO. 
        vtp_nSpan        = 1;       % Number of spanwise panels for whole wing
        vtp_nChord       = 1;       % Number of chordwise panels for whole wing%12
        vtp_addTwist     = true;     % Add twist to the Aero Mesh, 1 = YES, 0 = NO.
        vtp_addControls  = true;     % Add control surfaces, 1 = YES, 0 = NO. 
        wing_addCamber   = false;
        extendspoiler    = true;    % Extend the spoiler to the trailing edge or not.
    end
    
    methods
        
        function obj = Optsaero(~)
            
        end
        
    end
    
end

