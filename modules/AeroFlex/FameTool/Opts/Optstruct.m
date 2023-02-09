%% OPTSTRUCT Options for structure
%

%   Copyright 2016 University of Bristol
%   Private function.
classdef Optstruct < handle
    
    properties
        addEngine          = true;    % Add engine to the beam model, 1 = YES, 0 = NO.
        rbeOrient          = false;   % RBE0 Orientation 0 streamwise, 1 perpendicular
        remeshStruct       = [];      % Remesh the structure with eta and element numbers defined
        massCorrection     = true;    % Apply mass corrections from Excel file
        massCutPlane       = 'XZ';    % Mass cut plane parallel to XZ-plane ('XZ') or elastic axis ('EA')
         WingStartID        = 110000;
%         HTPStartID         = 120000;
%         VTPStartID         = 130000;
%         FuseStartID        = 100000;
%         StrucOffsetID      = 5000;
%         AeroOffsetID       = 5000;
%         LEOffsetID         = 101000;
%         TEOffsetID         = 102000;
        %WingStartID        = 2000;
        HTPStartID         = 4000;
        VTPStartID         = 3000;
        FuseStartID        = 1000;
        StrucOffsetID      = 500;
        AeroOffsetID       = 5000;
        LEOffsetID         = 100;
        TEOffsetID         = 200;
        WingCAEROID        = 710000;
        HTPCAEROID         = 810000;
        VTPCAEROID         = 910000;
        WingSetId          = 510001;
        HTPSetId          = 520001;
        VTPSetId          = 530001;
        addcarrythrough    = false;
    end
    
    methods
        function obj = Optstruct(~)
            %obj.remeshStruct.eta = [];
             obj.remeshStruct.eta = [0.0000,0.08339,0.2569,0.2868,0.5,0.7946,1.0000];
            %obj.remeshStruct.elementNumbers = [];
             obj.remeshStruct.elementNumbers = [5,17,3,15,10,10];
             
        end
    end
    
end

