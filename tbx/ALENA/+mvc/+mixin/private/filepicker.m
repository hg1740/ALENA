function [fn, idx] = filepicker(fn, spec, bPrompt, bPut, str)

%Prompt required ?
if bPrompt
    
    %Ask the user
    if bPut
        [f, p, idx] = uiputfile(spec, str, fn);
    else
        [f, p, idx] = uigetfile(spec, str, fn);
    end
    
    %Cancelled ?
    if isempty(f) || isequal(f, 0)
        
        %Bail out
        fn = [];
        idx = [];
        return;
        
    end
    
    %Make a note of full filename
    fn = fullfile(p,f);
    
else
    
    %Look at file extension
    [p, ~, e] = fileparts(fn);
    
    %If nothing
    if isempty(e)
        
        %Go with first from spec
        idx = 1;
        
        %Get extension(s)
        e = spec{1,1};

        %Allow for multiple
        e = strsplit(e,';');
        e = e{1};
        
        %Append extension to filename
        fn = [fn, e(2:end)];
        
    else
        
        %Can we work it out mask index ?
        idxx = find(cellfun(@(x)~isempty(strfind(x,e)), spec(:,1)));
        if numel(idxx) == 1
            
            %Yes
            idx = idxx;
            
        else
            error('unable to match file extension against any defined in spec');
        end
        
    end
    
    %If path not explicitly give
    if isempty(p)
        
        %If we are PUTTING a file
        if bPut
        
            %Go with current directory
            fn = fullfile(pwd, fn);
            
        else
            
            %Let 'which' have a go (without losing original name if not found)
            if ~isempty(which(fn))
                fn = which(fn);
            end
        
        end
        
    end
    
end

