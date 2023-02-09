function obj = fame2LoadCase(homeFolder,obj)

fp = fopen(fullfile(homeFolder,'FameModel_LoadCases.m'), 'w');

fprintf(fp,'%%%% Define the required load cases\n');
fprintf(fp,'\n');
fprintf(fp,'%%%% Inputs\n');
fprintf(fp,'%%  LoadCase.M      - Mach Number\n');
fprintf(fp,'%%  LoadCase.Alt    - Altitude (m)\n');
fprintf(fp,'%%  LoadCase.URDD3  - Vertical Acceleration (ms-2)\n');
fprintf(fp,'%%  LoadCase.Fuel   - Fuel Percentage\n');
fprintf(fp,'%%  LoadCase.Payload     - Kgs\n');
fprintf(fp,'%%  LoadCase.PayloadCG   - MAC\n');
fprintf(fp,'%%  LoadCase.CS.label\n');
fprintf(fp,'\n');
fprintf(fp,'\n');

MAXFUEL = max([obj.Fame.LoadCases.AC_mass]);
MINFUEL = min([obj.Fame.LoadCases.AC_mass]);

FUEL = MAXFUEL - MINFUEL;

for idx = 1:numel(obj.Fame.LoadCases)
   
    obj.Fame.LoadCases(idx).payload = obj.Fame.LoadCases(idx).AC_mass - obj.Fame.LoadCases(idx).fuel - obj.Fame.Geometry.Weights.oew; 
    FF = obj.Fame.LoadCases(idx).fuel/FUEL;
    
    fprintf(fp,'LoadCases(%i).Label     = '' %s '';\n',idx,'');
    fprintf(fp,'LoadCases(%i).Type      = '' %s '';\n',idx,'');
    fprintf(fp,'LoadCases(%i).M         = %-8.4f ;\n',idx,obj.Fame.LoadCases(idx).Mach);
    fprintf(fp,'LoadCases(%i).Alt       = %-8.4f ;\n',idx,obj.Fame.LoadCases(idx).Altitude);
    fprintf(fp,'LoadCases(%i).URDD3     = %-8.4f ;\n',idx,obj.Fame.LoadCases(idx).Loadfactor * 9.81);
    fprintf(fp,'LoadCases(%i).Payload   = %-8.4f ;\n',idx,obj.Fame.LoadCases(idx).AC_mass - obj.Fame.LoadCases(idx).fuel - obj.Fame.Geometry.Weights.oew);
    fprintf(fp,'LoadCases(%i).PayloadCG = %-8.4f ;\n',idx,obj.Fame.LoadCases(idx).AC_cg);
    fprintf(fp,'LoadCases(%i).Fuel      = %-8.4f ;\n',idx,FF); %TODO FIX THIS!
    
    % Sort out the control surfaces
    ControlCount = 0;
    
    ControlAngle = [];
    ControlLabel = [];
    ailcount = 0;
    splcount = 0;
    for j = 1:numel(obj.Fame.LoadCases(idx).aileron)
        ailcount = ailcount + 1;
        ControlCount = ControlCount + 1;
        label = ['aileron' num2str(j)'];
        ControlLabel{ControlCount} = ['ail' num2str(ailcount) '1r'];
        
        ControlAngle = [ControlAngle, obj.Fame.LoadCases(idx).aileron(j)];
                
    end
    
    for j = 1:numel(obj.Fame.LoadCases(idx).spoiler)
        ControlCount = ControlCount + 1;
        splcount = splcount + 1; 
        label = ['spoiler' num2str(j)];
        ControlLabel{ControlCount} = ['spl' num2str(splcount) '1r'];
        
        ControlAngle = [ControlAngle, obj.Fame.LoadCases(idx).spoiler(j)];
        
    end
    
    % Write the control surface labels
    fprintf(fp,'LoadCases(%i).CS.Label  = {',idx);
    for j = 1:ControlCount
        if j == 1
            fprintf(fp,'''%s''' ,ControlLabel{j});
        else
            fprintf(fp,',''%s''' ,ControlLabel{j});
        end
    end
    fprintf(fp,'};\n');
    
    % Write the control surface deflections
    fprintf(fp,'LoadCases(%i).CS.Value  = [',idx);
    for j = 1:ControlCount
        if j == 1
            fprintf(fp,'%8.4f' ,ControlAngle(j));
        else
            fprintf(fp,', %8.4f' ,ControlAngle(j));
        end
    end
    fprintf(fp,'];\n');
    
    fprintf(fp,'\n');
end

fclose(fp);
end