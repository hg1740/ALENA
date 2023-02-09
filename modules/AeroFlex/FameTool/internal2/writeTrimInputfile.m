function writeTrimInputfile(homeFolder,obj)

% Run through the trim cases and generate an input file for each of them.

for trimI = 1:numel(obj.Fame.LoadCases)
    
    fp = fopen(fullfile(homeFolder,['\Neo_model\TrimInput_LC' num2str(trimI) '.dat']), 'w');
    fprintf(fp,'SOL 144\n');
    fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
    fprintf(fp,'$------- AERO SURFACE PROPERTIES -----------------------------------------\n');
    fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
    
    AEROS.Span = max(struct2mat(obj.Mdl.Grid,'coord',2)) * 2;
    AEROS.SurfaceArea = trapz((AEROS.Span / 2) * obj.Fame.Geometry.Wing.chord(:,1),obj.Fame.Geometry.Wing.chord(:,2));
    AEROS.Chord = obj.Fame.Geometry.Wing.chord(1,2);
    AEROS.SIMXZ = abs(obj.Opts.Geom.reflect - 1);
    AEROS.SIMYZ = 0;
    
    ch = sprintf('%-8s', num2_8ch(AEROS.Chord));
    if obj.Opts.Geom.reflect == 1
        sp = sprintf('%-8s', num2_8ch(AEROS.Span));
        sa = sprintf('%-8s', num2_8ch(2*AEROS.SurfaceArea));
    else
        sp = sprintf('%-8s', num2_8ch(AEROS.Span / 2));
        sa = sprintf('%s', num2_8ch(AEROS.SurfaceArea));
    end
    xz = sprintf('%-8g', AEROS.SIMXZ);
    yz = sprintf('%-8g', AEROS.SIMYZ);
    fprintf(fp,'AEROS                   %s%s%s%s%s\n', ch, sp, sa, xz, yz);
    
    fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
    fprintf(fp,'$------- SPC PROPERTIES --------------------------------------------------\n');
    fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
    
    spcCases = unique([obj.Mdl.Spc.id]);
    
    for i = 1:numel(spcCases)
        fprintf(fp,'SPC= %1.0f\n',spcCases(i));
    end
    
    for i = 1:numel(obj.Mdl.Spc)
        
        str = genNeocassDeck(obj.Mdl.Spc(i));
        
        for j = 1:numel(str)
            fprintf(fp,'%s\n',str{j});
        end
    end
    
    fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
    fprintf(fp,'$------- TRIM CASES ------------------------------------------------------\n');
    fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
    
    % Populate the body conditions
    TRIM.id = obj.Fame.LoadCases(trimI).LC;
    TRIM.MCH = obj.Fame.LoadCases(trimI).Mach;
    TRIM.ALT = obj.Fame.LoadCases(trimI).Altitude;
    TRIM.URDD3 = obj.Fame.LoadCases(trimI).Loadfactor * 9.81;
    TRIM.MASS = obj.Fame.LoadCases(trimI).AC_mass;
    TRIM.SIDES = 0;
    TRIM.ROLL = 0;
    TRIM.YAW = 0;
    TRIM.URDD2 = 0;
    TRIM.THRUST = 0;
    TRIM.HEAD = 0;
    TRIM.URDD4 = 0;
    TRIM.URDD6 = 0;
    TRIM.CLIMB = 0;
    TRIM.BANK = 0;
    TRIM.URDD5 = 0;
    TRIM.URDD1 = 0;
    TRIM.PITCH = 0;
    
    id      = sprintf('%-8g', TRIM.id);
    ran     = sprintf('%-8g', 1);
    mch     = sprintf('%-8g', TRIM.MCH);
    alt     = sprintf('%-8g', TRIM.ALT);
    sides   = sprintf('%-8g', TRIM.SIDES);
    roll    = sprintf('%-8g', TRIM.ROLL);
    yaw     = sprintf('%-8g', TRIM.YAW);
    urdd2   = sprintf('%-8g', TRIM.URDD2);
    thrust  = sprintf('%-8g', TRIM.THRUST);
    head    = sprintf('%-8g', TRIM.HEAD);
    urdd4   = sprintf('%-8g', TRIM.URDD4);
    urdd6   = sprintf('%-8g', TRIM.URDD6);
    urdd3   = sprintf('%-8g', TRIM.URDD3);
    climb   = sprintf('%-8g', TRIM.CLIMB);
    bank    = sprintf('%-8g', TRIM.BANK);
    urdd5   = sprintf('%-8g', TRIM.URDD5);
    pitch   = sprintf('%-8g', TRIM.PITCH);
    urdd1   = sprintf('%-8g', TRIM.URDD1);
    fprintf(fp,'\n$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
    fprintf(fp,'TRIM= %s\n', id);
    fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
    fprintf(fp,'TRIM    %s%s%s%sSIDES   %sROLL    %s\n', id, ran, mch, alt, sides, roll);
    fprintf(fp,'        YAW     %sURDD2   %sTHRUST  %sHEAD    %s\n', yaw, urdd2, thrust, head);
    fprintf(fp,'        URDD4   %sURDD6   %sURDD3   %sCLIMB   %s\n', urdd4, urdd6, urdd3, climb);
    fprintf(fp,'        BANK    %sURDD5   %sPITCH   %sURDD1   %s\n', bank, urdd5, pitch, urdd1);
    
    % Sort out the control surfaces
    ControlCount = 0;
    idxail = [];
    for j = 1:numel(obj.Fame.LoadCases(trimI).aileron)
        ControlCount = ControlCount + 1;
        label = ['aileron' num2str(j)'];
        TRIM.(label) = obj.Fame.LoadCases(trimI).aileron(j);
        ail = sprintf('%-8.2f', obj.Fame.LoadCases(trimI).aileron(j));
        idxail(j).idx = find(strcmpi({obj.Mdl.Caero.csType},label));
        fprintf(fp,'        %s%s', strjust(sprintf('%8s',obj.Mdl.Caero(idxail(j).idx(1)).csId),'left'),ail);
        
    end
    
    idxspl = [];
    for j = 1:numel(obj.Fame.LoadCases(trimI).spoiler)
        ControlCount = ControlCount + 1;
        label = ['spoiler' num2str(j)];
        TRIM.(label) = obj.Fame.LoadCases(trimI).spoiler(j);
        spl = sprintf('%-8.2f', obj.Fame.LoadCases(trimI).spoiler(j));
        idxspl(j).idx = find(strcmpi({obj.Mdl.Caero.csType},label));
        if mod(ControlCount-1,4) == 0
            fprintf(fp,'\n        %s%s', strjust(sprintf('%8s',obj.Mdl.Caero(idxspl(j).idx(1)).csId),'left'),spl);
        else
            fprintf(fp,'%s%s', strjust(sprintf('%8s',obj.Mdl.Caero(idxspl(j).idx(1)).csId),'left'),spl);
            
        end
    end
    
    if obj.Opts.Geom.addTail
        idxrdr = [];
        for j = 1:numel(obj.Fame.LoadCases(trimI).rudder)
            ControlCount = ControlCount + 1;
            label = ['rudder' num2str(j)];
            TRIM.(label) = obj.Fame.LoadCases(trimI).rudder(j);
            rdr = sprintf('%-8.2f', obj.Fame.LoadCases(trimI).rudder(j));
            idxrdr(j).idx = find(strcmpi({obj.Mdl.Caero.csType},label));
            if mod(ControlCount-1,4) == 0
                fprintf(fp,'\n        %s%s', strjust(sprintf('%8s',obj.Mdl.Caero(idxrdr(j).idx(1)).csId),'left'),rdr);
            else
                fprintf(fp,'%s%s', strjust(sprintf('%8s',obj.Mdl.Caero(idxrdr(j).idx(1)).csId),'left'),rdr);
                
            end
        end
    end
    
    fclose(fp);
end

end
