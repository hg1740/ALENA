%% DIRS Paths to NeoCass and FAME folders, as well as output directory
classdef Dirs
    
    properties
        homeFolder     = neocassroot;
        % --------------------
        % HARTEN
        % --------------------
        
        %         TUXFolder      = 'D:\GITDirectories\AWI_Fame\HARTEN_v2_VLM\01_geometry';
        %         fameFolder     = 'D:\GITDirectories\AWI_Fame\HARTEN_v2_VLM\06_weights\FAME\v3_newCG_LC.fame-w.4.00.m1-1-pv8';
        %         fameFuelFolder = 'D:\GITDirectories\AWI_Fame\HARTEN_v2_VLM\06_weights\FAME\v3_newCG_LC.fame-w.4.00.m1-1-pv8\__mon\wbd2\massdist\flexible';
        %         airfoilFolder  = 'D:\GITDirectories\AWI_Fame\HARTEN_v2_VLM\02_aerodynamics\AWIC5_profiles';
        MainFolder     = 'E:\GITDirectories\AWI_Fame\HARTEN_v2_VLM';
        TUXFolder      = '01_geometry';
        fameFolder     = '06_weights\FAME\*.m1-1-pv8';
        fameFuelFolder = '__mon\wbd2\massdist\flexible';
        airfoilFolder  = '02_aerodynamics\AWIC5_profiles';
        
        
        
        %         %
        %         TUXFolder      = '/Volumes/UOBSTICK/GITDirectories/FAMEInputModel/01_geometry';
        %         fameFolder     = '/Volumes/UOBSTICK/GITDirectories/FAMEInputModel/06_weights/FAME/newCG_LC.fame-w.4.00.m1-1-pv8';
        %         fameFuelFolder = '/Volumes/UOBSTICK/GITDirectories/FAMEInputModel/06_weights/FAME/newCG_LC.fame-w.4.00.m1-1-pv8/__mon/wbd2/massdist/flexible';
        
        % --------------------
        % WOTF
        % --------------------
        %         TUXFolder      = '[]';
        %         fameFolder     = 'D:\Models\WoTF\WoTF_CFRP_v3.3_06.1_Update_Ribs_and_Blade_DC.fame-w.4.00.m1-1-pv8';
        %         fameFuelFolder = 'D:\Models\WoTF\WoTF_CFRP_v3.3_06.1_Update_Ribs_and_Blade_DC.fame-w.4.00.m1-1-pv8\__mon\wbd2\massdist\flexible';
        %         airfoilFolder  = '[]';
        %
        writeFolder    = 'UoBFiles';
        
    end
    
    methods
        
        function val = get.TUXFolder(obj)
            
            val = fullfile(obj.MainFolder,obj.TUXFolder);
            
        end
%         
%         function val = get.fameFolder(obj)
%             
%             d = dir(fullfile(obj.MainFolder,obj.fameFolder));
%             
%             if ~isempty(d)
%                 val = fullfile(d(1).folder,d(1).name);
%             else
%                 val = '';
%             end
%             
%         end
%         
%         function val = get.fameFuelFolder(obj)
%             
%             val = fullfile(obj.fameFolder,obj.fameFuelFolder);
%             
%         end
%         
%         function val = get.airfoilFolder(obj)
%             
%             val = fullfile(obj.MainFolder,obj.airfoilFolder);
%             
%         end
%         
%         function val = get.writeFolder(obj)
%             
%             val = fullfile(obj.MainFolder,obj.writeFolder);
%             
%         end
        
    end
    
end

