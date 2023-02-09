function [sn, args] = sheetpicker(fn, args, bPrompt, bPut, str)

%Check file extension
[~,~,e] = fileparts(fn);

%Are we bothered ?
if ~any(strcmpi(e, {'.xls', '.xlsx', '.xlsm'}))
    
    %No
    sn = [];
    return;
    
end

%Look for sheet in args
b = cellfun(@(x)ischar(x) && strcmp(x, 'Sheet'), args);
if any(b)
    
    %Get the corresponding value
    sn = args{find(b) + 1};
    
    %Strip from args
    args(find(b) + [0,1]) = [];
    
    %Caller may be explicitly telling us to prompt
    if strcmp(sn, 'ask')
        bPrompt = true;
        sn = [];
    end
    
else
    
    %We have to prompt
    bPrompt = true;
    sn = [];
    
end
                
%What have we got already ?
if exist(fn, 'file') == 2
    
    %Ask the file
    [~, shts] = xlsfinfo(fn);
    
    %If we're putting
    if bPut
        
        %Add option for new sheet
        shts{end+1} = 'New sheet...';
        
    end
    
else
    shts = {};
end

%Caller may be saying ALL sheets
if strcmp('all', sn)
    
    %We're done
    sn = shts;
    return;
    
end

%Prompt required ?
if bPrompt && ~isempty(shts)
    
    %Ask the user
    [sel, bOK] = listdlg('ListString', shts, ...
        'PromptString', 'Select sheet(s)...', ...
        'Name', str);
    
    %Cancelled ?
    if ~bOK
        
        %Bail out - returning FALSE rather than empty
        % (so caller can tell difference)
        sn = false;
        return;

    elseif numel(sel) == 1
        
        %Unpick
        sn = shts{sel};
        bPrompt = false;

    else
        
        %Send back multiple sheets
        sn = shts(sel);
        bPrompt = false;
        
    end

end

%Still need to ask for sheet name ?
if bPrompt && isempty(shts) || all(strcmp(sn, 'New sheet...'))
    
    %Yes
    val = inputdlg({'New sheet name:'}, str, 1, {''}, 'on');
    
    %Cancelled ?
    if isempty(val)
        
        %Bail - returning FALSE rather than empty
        % (so caller can tell difference)
        sn = [];
        
    else
        
        %Unpick
        sn = val{1};
        
    end
    
end
