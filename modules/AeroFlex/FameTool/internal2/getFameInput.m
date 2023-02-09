%GETFAMEINPUT: A script that reads the *.fm4 file and creates a structure
%that has a similar format the file
%
% ---------------------------------------------------
% Input:     - Fame filename ('*.fm4')
%
% Output:    - Output (structure)
% ---------------------------------------------------
%
% TODO: -Add @include functionality to the script
%       -Fix loadcase definitions (it currently overwrites if the loadtype
%       is the same)
%       -Fix subfields
%
% Date of Creation: 2017
%
% Authors:  Dario Calderon - University of Bristol
%           Etienne Coetzee - Airbus (FPO)
% ---------------------------------------------------
function Output = getFameInput(filename, fcn)

if nargin == 0
    [fName,pName,idx] = uigetfile('Pick a FAME Input file');
    if idx~=0
        filename = [pName,fName];
    else
        error('Unable to proceed without a valid file');
    end
end
if nargin < 2
    fcn = @(s) fprintf('%s\n', s);
end

% ----------------------------------------
% Read the *.fm4 file and convert to a cell of strings
% ----------------------------------------
s = textread(filename,'%s','delimiter','\n','bufsize',10000);
%s = textscan(filename,'%s','delimiter','\n');

NumLines = size(s,1);
Extra_node = 0;
Tree_loc = -1;
checkfields = [];
entry_count = 0;

% ----------------------------------------
% Run through each line of the input file and begin to generate the
% structure
% ----------------------------------------
for i = 1:NumLines
    
    isexc = regexp(s{i,1},'!');
    
    if isexc
        s{i,:} = [s{i,:}(1:isexc-1),''];
    end
    
    if ~isempty(s{i,:}) % Skip line if the line is blank
        
        if ~strcmp(s{i,1}(1),'!') && ~strcmp(s{i,1}(1),'@') % Skip line if the line beings with '!' or '@'
            
            F2 = [];
            F2_end = [];
            
            % Search for field name e.g. <<<< FAME_01_ >>>>
            F0 = strtrim(regexp(s{i,:},'(?<=<<<<).*(?=_\s?>>>>)','match'));
            
            if isempty(F0)
                F1 = strtrim(regexp(s{i,:},'(?<=<<<).*(?=_\s?>>>)','match'));
                F1 = strrep(F1,'-','_');
                if ~isempty(F1)
                    Tree_loc = 1;
                    Output.(F1{1,1}) = [];
                    PrimField = F1;
                else
                    F2 = strtrim(regexp(s{i,:},'(?<=<<)(\s?).*(?=_\s?>>)','match'));
                    F2 = strrep(F2,'-','_');
                    %F2 = strtrim(regexp(s{i,:},'(?<=<<).*(?=_ >>)','match'));
                    if ~isempty(F2)
                        Tree_loc = 2;
                        Output.(PrimField{1,1}).(F2{1,1}) = [];
                        SecField = F2;
                    end
                end
                
                F1_end = strtrim(regexp(s{i,:},'(?<=<<<\s?_).*(?=>>>)','match'));
                if ~isempty(F1_end)
                    Tree_loc = 0;
                else
                    
                    F2_end = strtrim(regexp(s{i,:},'(?<=<<\s?_).*(?=>>)','match'));
                    if ~isempty(F2_end)
                        Tree_loc = 1;
                    end
                end
                
                if isempty(F1) && isempty(F2) && isempty(F1_end) && isempty(F2_end)
                    F_parent    = strtrim(regexp(s{i,:},'(?<=<).*(?=_\s?>)','match'));
                    F_parent_end= strtrim(regexp(s{i,:},'(?<=<)(\s?)_.*(?=>)','match'));
                    F_string    = strtrim(regexp(s{i,:},'(?<=<).*(?=>)','match'));
                    F_value     = strtrim(regexp(s{i,:},'(?<=[).*(?=])','match'));
                    F_stringAmp = strtrim(regexp(s{i,:},'(?<=<(\s*)&).*(?=>)','match'));
                    
                    if ~isempty(F_parent)
                        Extra_node = 1;
                        F_parent_string = F_parent;
                    end
                    
                    if ~isempty(F_parent_end)
                        Extra_node = 0;
                        entry_count = 0;
                        checkfields = [];
                    end
                    
                    if Extra_node == 0 && isempty(F_parent_end)
                        if ~isempty(F_string)
                            if ~isempty(F_stringAmp)
                                F_string = F_stringAmp;
                            end
                            
                            if length(F_string{1,1})>30
                                F_string = regexp(s{i,:},'(?<=<)\w*(?=>)','match');
                            end
                            
                            if Tree_loc == 0
                                Output.(F_string{1,1}) = F_value;
                            elseif Tree_loc == 1
                                Output.(PrimField{1,1}).(F_string{1,1}) = F_value;
                            elseif Tree_loc == 2
                                if strcmp(SecField{1,1}, lcTok)
                                    Output.(PrimField{1,1}).(SecField{1,1})(lcIndex).(F_string{1,1}) = F_value;
                                else
                                    Output.(PrimField{1,1}).(SecField{1,1}).(F_string{1,1}) = F_value;
                                end
                            end
                            checkfields = [];
                        end
                        
                    elseif Extra_node == 1 && isempty(F_parent)
                        if ~isempty(checkfields)
                            entry_count = entry_count + 1;
                            op_br =  strfind(s{i,:},'[');
                            cl_br =  strfind(s{i,:},']');
                            UpdString = s{i,:};
                            for ii = 1:length(op_br)
                                UpdString(op_br(ii):cl_br(ii)) = strrep(s{i,:}(op_br(ii):cl_br(ii)),' ','_');
                            end
                            F_value = regexp(UpdString,'\s*','split');
                            for ii = 1:length(checkfields)
                                F_string = strtrim(regexp(checkfields{ii},'(?<=<).*(?=>)','match'));
                                if length(F_string{1,1})>30
                                    F_string = regexp(checkfields{ii},'(?<=<)\w*(?=>)','match');
                                end
                                %Edit: C.Szczyglowski, 28/10/2018
                                % Treat the load cases differently
                                %   - We can get into this area of the code
                                %   when we want to assign data to the
                                %   embedded structures within the loadcase
                                %   entry. Check for this and action as
                                %   appropriate.
                                if strcmp(SecField{1,1}, 'LOAD_CASES')
                                    Output.(PrimField{1,1}).(SecField{1,1})(lcIndex).(F_parent_string{1,1}).(F_string{1,1}){entry_count} = F_value{ii};
                                else
                                    Output.(PrimField{1,1}).(SecField{1,1}).(F_parent_string{1,1}).(F_string{1,1}){entry_count} = F_value{ii};
                                end
                            end
                        else
                            entry_count = 0;
                            checkfields = regexp(s{i,:}, '<\s*\w*\s*>','match');
                            if length(checkfields) > 1
                                for ii = 1:length(checkfields)
                                    F_string = strtrim(regexp(checkfields{ii},'(?<=<).*(?=>)','match'));
                                    %Edit: C.Szczyglowski, 28/10/2018
                                    % Treat the load cases differently.
                                    %   - We get into this area of the code
                                    %   when there is an embedded structure
                                    %   within the loadcase entry. Check
                                    %   for this and action as appropriate
                                    if strcmp(SecField{1,1}, 'LOAD_CASES')
                                        Output.(PrimField{1,1}).(SecField{1,1})(lcIndex).(F_parent_string{1,1}).(F_string{1,1}) = [];                                        
                                    else
                                        Output.(PrimField{1,1}).(SecField{1,1}).(F_parent_string{1,1}).(F_string{1,1}) = [];
                                    end
                                end
                            else
                                checkfields = [];
                                if Tree_loc == 0
                                    Output.(F_string{1,1}) = F_value;
                                elseif Tree_loc == 1
                                    Output.(PrimField{1,1}).(F_string{1,1}) = F_value;
                                elseif Tree_loc == 2
                                    %Edit: C.Szczyglowski, 28/10/2018
                                    %
                                    % Treat the loadcases differently. 
                                    %   - Loadcases are enclosed in a top
                                    %   level header that defines their
                                    %   type, e.g. MANOEUVRE, 
                                    %   PRATT_GUST_FAR25CH13, etc. We want
                                    %   to store this information as  the
                                    %   load case type and have a structure
                                    %   array of loadcases.
                                    lcTok = 'LOAD_CASES';
                                    if ~strcmp(SecField{1,1}, lcTok)
                                        Output.(PrimField{1,1}).(SecField{1,1}).(F_parent_string{1,1}).(F_string{1,1}) = F_value;                                        
                                    else
                                        
                                        %Expected load case fields
                                        fNames = {'ID'             , 'LABEL', 'TYPE'          , ...
                                            'LOAD_FACTOR'          , 'ZFW'  , 'FUEL_MASS_WING', ...
                                            'WING_FUELING_SEQUENCE', 'MACH' , 'ALTITUDE'      , ...
                                            'ALLOWABLES_SET_ID'    ,'CG_POS', 'COMMENT'       , ....
                                            'FLAP_DEFLECTION'      , ...
                                            'SPOILER_DEFLECTION'   , ...
                                            'AILERON_DEFLECTION'   , 'GAMMA4', 'ROLL_RATE'};
                                        
                                        %Store the type & ID number 
                                        %   - TODO: This relies on ID being
                                        %     the first field in the .fm4
                                        %     loadcase structure.
                                        if ~contains(F_parent_string{1,1}, 'DEFLECTION')
                                            lcType = F_parent_string{1, 1};
                                        end
                                        if strcmp(F_string{1,1}, 'ID')
                                            lcID   = F_value{1,1};
                                        end
                                        
                                        %Initiate the structure & add to
                                        %the parent structure
                                        lc = cell2struct(cell(size(fNames)), fNames, 2);
                                        if isfield(Output.(PrimField{1,1}), lcTok) && isstruct(Output.(PrimField{1,1}).(lcTok))
                                            %Use the load case 'ID' to
                                            %determine if we have moved to
                                            %a new load case.
                                            if lcID == Output.(PrimField{1,1}).(lcTok)(end).ID
                                                lcIndex = numel(Output.(PrimField{1,1}).(lcTok));
                                            else
                                                lcIndex = numel(Output.(PrimField{1,1}).(lcTok)) + 1;
                                                Output.(PrimField{1,1}).(lcTok)(lcIndex) = lc;
                                            end                                            
                                        else     
                                            %First time through?
                                            Output.(PrimField{1,1}).(lcTok) = lc;
                                            lcIndex = 1;
                                        end
                                                                                
                                        %Add data to the structure
                                        Output.(PrimField{1,1}).(lcTok)(lcIndex).TYPE = lcType;
                                        Output.(PrimField{1,1}).(lcTok)(lcIndex).(F_string{1,1}) = F_value{1,1};
                                        
                                    end
                                end
                            end
                        end
                    end
                end
            end
        elseif strcmp(s{i,1}(1),'@')
            fcn([blanks(8) 'LINE IGNORED: ' s{i,:}]);
        end
    end
    
end
end