%% GETFAMEBOXSTIFFFCN Obtain box stiffness along span from FAME

function [Stiffness,Inputfiles] = getFameStress(folderName)

header_Strings = {'Spannungen in der Oberschale','Spannungen in der Unterschale','Spannungen in den Holmen',...
     'Normalspannung in den Stringern der Oberschale','Normalspannung in den Stringern der Unterschale','Normalkraft in den Wandelementen der Oberschale',...
     'Normalkraft in den Stringern der Oberschale','Normalkraft in den SuperStringern der Oberschale','Normalkraft in den Wandelementen der Unterschale',...
     'Normalkraft in den Stringern der Unterschale','Normalkraft in den SuperStringern der Unterschale'};

component_strings = {'Tau_t','Tau_s','Tau','Tau_s_eff','Tau_eff','Sig_ms','Sig_ms_eff','Sig_allow','Tau_allow','Sig_buck','Tau_buck','RF_strength','RF_buck','Sig_w'};

pathName = ['1_structure',filesep,'1_10_wing',filesep,'l3_fame-w',filesep,'flexible',filesep,'stresses'];

list = dir(fullfile(folderName,pathName,'*mean_stresses.str'));

for i = 1:length(list)
    
    fileName = list(i).name;
    
    FilePath = [folderName,filesep, pathName,filesep, fileName];
    
    Inputfiles.stressFiles{i} = FilePath;
    
    if ~exist(Inputfiles.stressFiles{i},'file')
        error('Unable to find file %s.\n', Inputfiles.stressFiles{i});
    else
        
        s      = textread(Inputfiles.stressFiles{i} ,'%s','delimiter','\n');

        idx = [];
        for j = 1:numel(header_Strings)
            idx(j)    = find(~cellfun(@isempty,regexp(s,header_Strings{j})));
        end
        
        % UpperSkin
        skin_u = s(idx(1):idx(2)-1);
        % Identify the component strings and remove them as I go a long
        
        % Lower Skin
        skin_l_idx = s(idx(2):idx(3)-1);
        
        % Spars
        spar_idx   = s(idx(3):idx(4)-1);
        
        idx    = find(strncmp(s,'%',1));
        s(idx) = [];
        s      = str2num(char(s));
        
        Stiffness.eta  = s(:,1);
        Stiffness.EIxx = s(:,2)/1000000;
        Stiffness.EIxz = s(:,3)/1000000;
        Stiffness.EIzz = s(:,4)/1000000;
        Stiffness.EKxx = s(:,5)/1000000;
        Stiffness.EKzz = s(:,6)/1000000;
        Stiffness.GJ   = s(:,7)/1000000;
        Stiffness.EA   = s(:,8);
        Stiffness.Abox = s(:,9);
        Stiffness.Xec  = s(:,10)/1000;
        Stiffness.Zec  = s(:,11)/1000;
        Stiffness.Xsc  = s(:,12)/1000;
        Stiffness.Zsc  = s(:,13)/1000;
        
    end
end

end

