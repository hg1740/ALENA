%% FAMEAEROCONVERTER   Extract FAME aerodynamic definitions and create aero 
%                      panels
%
% Extract the FAME aerodynamic definition and create aerodynamic panels for
% NEOCASS and NASTRAN. A FAME model object is passed to this function from
% the FAME2MAT object.
%
% See NEOCASS and NASTRAN manuals for further explanations of the
% aerodynamic decks.
%
% Make sure the following inputs are defined in the FAME2MAT object
% 
%     obj.Dirs.fameFolder       : FAME folder
%     obj.Inp.aeroPlanformFiles : auto generated
%     obj.Opts.Aero.nSpan       : Number of spanwise panels
%     obj.Opts.Aero.nChord      : Number of chordwise panels
% 
%   Numbering conventions:
%     Right hand wing : 710001
%     Left hand wing  : 720001
%
%   See also |getFameAero| method within Fame2mat object.

%   Copyright 2016 University of Bristol
%   Private function.
function [results,Inputfiles] = getFameWingOutline(folderName)

pathName    = ['export',filesep,'fmwpp01',filesep,'flexible'];

fileName    = 'wng_planform.fmwpp01';

Outline     = 0;

Inputfiles.aeroPlanformFiles{1} = fullfile(folderName,pathName,fileName);

keyword = { '%      x.outline       y.outline       z.outline'};

noutline = 0;

if ~exist(Inputfiles.aeroPlanformFiles{1},'file')
    error('Unable to find file %s.\n', Inputfiles.aeroPlanformFiles{1});
else
    fp       = fopen(Inputfiles.aeroPlanformFiles{1}, 'r');
    skipLine = false;
    
    while ~feof(fp)
        if ~skipLine
            tline = fgetl(fp);
        else
            skipLine = false;
        end
        CARD = tline;
        
        switch CARD
            % Outline card
            case keyword{1}
                fgetl(fp);
                fgetl(fp);
                tline    = fgetl(fp);
                noutline = noutline + 1;
                Outline(noutline,1:3) = strread(tline);
                tline    = fgetl(fp);
                
                % continuation detected
                while ~isequal(tline(1),'%')
                    noutline = noutline + 1;
                    Outline(noutline,1:3) = strread(tline);
                    tline    = fgetl(fp);
                end
                
                Outline  = Outline / 1000;
                skipLine = false;
        end
    end
    
    fclose(fp);
end

results = Outline;

end




