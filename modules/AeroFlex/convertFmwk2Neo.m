%%CONVERT2Aeroflex: Preliminary function to take the UOB Framework object and
%%begin a sizing using AeroFlex

function [beam_model,Aircraft_param,Optim_Variables,LoadCases] = convertFmwk2Neo(LSObj,BBObj,ac,lc,numelements, logfcn)

%Log function not supplied ?
if nargin < 6
    
    %Go with good old disp
    logfcn = @disp;
    
end

% Just in case!
if nargin < 5
    numelements = 30;
end

if numelements > 60
    numelements = 60;
    logfcn('Number of elements has been reset to 60');
end

% Let's deal with the Target Cl
TargetCl_eta = {lc.TargetCl_eta};
TargetCl     = {lc.TargetCl};
TargetClIdx  = find(~cellfun(@isempty,TargetCl_eta));
if ~isempty(TargetClIdx)
    TargetCl = [TargetCl_eta{TargetClIdx};TargetCl{TargetClIdx}]';
else
    TargetCl = [];
end
%% LOAD LIFTING SURFACES:

PartNames   = {LSObj.Name};
portidx     = cellfun(@isempty,regexpi(PartNames,'port'));

LSObj = LSObj(portidx);

% Store the control surfaces 
ControlSurf = {};
for i = 1:numel(LSObj)
    
    Part = LSObj(i).Name;
    
    % This has been hardcoded in for the time being ensuring that only the
    % 'Wing' is being sized.
    isWing = ~isempty(regexp(Part,'Wing', 'ONCE'));
    
    if isWing
        Part = 'Wing';
        %WingPart = Part; 
    end
      
    %% AIRCRAFT PARTS
    Aircraft_param.(Part).LiftingSurface  = 1;

    Aircraft_param.(Part).Optimisation    = isWing;
    
    % Grab all front spar objects
    Spars = find(LSObj(i).Children,'Type','Spar');
    
    %Grab all control surfaces for this lifting surface
    CS = LSObj(i).ControlSurfaces;
    
    %Find all slats/flaps and remove them
    idx_slat = arrayfun(@(cs) isa(cs, 'awi.model.Slat'), CS);
    idx_flap = arrayfun(@(cs) isa(cs, 'awi.model.Flap'), CS);    
    CS(or(idx_slat, idx_flap)) = [];
    
    % This assumes that the span vector of the VTP will be defined as 'Y'
    SymPlane = ~strcmp(LSObj(i).SpanVector,'Y');
    
    Aircraft_param.(Part).Planform.SymPlane  = SymPlane;
    Aircraft_param.(Part).Parent             = LSObj(1).Parent.Name;
    
    Aircraft_param.(Part).Box.fspar.eta   = Spars(1).Eta';
    Aircraft_param.(Part).Box.fspar.value = Spars(1).XLoc';
    
    Aircraft_param.(Part).Box.aspar.eta   = Spars(2).Eta';
    Aircraft_param.(Part).Box.aspar.value = Spars(2).XLoc';
    
    Aircraft_param.(Part).Planform.c_ref = LSObj(i).RootChord;
%     Aircraft_param.(Part).Planform.c_ref = getBPV(LSObj(i), 'Chord');
    Aircraft_param.(Part).Planform.b_ref = LSObj(i).Span;
    Aircraft_param.(Part).Planform.AR    = [];
    
    %Aircraft_param.(Part).Planform.S_ref = LSObj(i).SurfaceArea; % Dario, do you want this to be populated or not?
    Aircraft_param.(Part).Planform.S_ref = [];
    
    Aircraft_param.(Part).Planform.NumCPanels = 1; % FIX
    
    Aircraft_param.(Part).Planform.NWingSec = numel(LSObj(i).SegLength);
    
    Aircraft_param.(Part).Planform.taper.eta_beg   = LSObj(i).Eta_(1:end-1)';
    Aircraft_param.(Part).Planform.taper.value_beg = (LSObj(i).Chord_(2:end)./LSObj(i).Chord_(1:end-1))';
    
    Aircraft_param.(Part).Box.height_fraction.eta   = [0;1];
    Aircraft_param.(Part).Box.height_fraction.value = [0.15;0.1];

    Aircraft_param.(Part).Planform.dihedral.eta_beg   = LSObj(i).Eta_(1:end-1)';
    Aircraft_param.(Part).Planform.dihedral.value_beg = LSObj(i).Dihedral_(1:end-1)';
    
    Aircraft_param.(Part).Planform.dihedral.eta_end   = LSObj(i).Eta_(2:end)';
    Aircraft_param.(Part).Planform.dihedral.value_end = LSObj(i).Dihedral_(1:end-1)';
    
    if strcmpi(Part,'vtp')
        Aircraft_param.(Part).Planform.dihedral.value_beg(:) = 90;
        Aircraft_param.(Part).Planform.dihedral.value_end(:) = 90;
    end
    
    Aircraft_param.(Part).Planform.QCsweep.eta_beg   = LSObj(i).Eta_(1:end-1)';
    Aircraft_param.(Part).Planform.QCsweep.value_beg = LSObj(i).Sweep_(1:end-1)';
    
    Aircraft_param.(Part).Planform.QCsweep.eta_end   = LSObj(i).Eta_(2:end)';
    Aircraft_param.(Part).Planform.QCsweep.value_end = LSObj(i).Sweep_(1:end-1)';
    
    Aircraft_param.(Part).Planform.QCsweep.loc       = LSObj(i).SweepLoc;
    
    Aircraft_param.(Part).Planform.twist.eta   = [0;1];
    Aircraft_param.(Part).Planform.twist.value = [3;3];
    
    Aircraft_param.(Part).Planform.Offset = LSObj(i).AbsPosition;
    
    Aircraft_param.(Part).Planform.Beamsweep = []; %Dario, do you want this to be populated? We can access this information!
    
    Aircraft_param.(Part).Planform.LEsweep = [];   %Dario, do you want this to be populated? We can access this information!
    
    if ~isWing
        Aircraft_param.(Part).Planform.NBeams = 1; % FIX
    else
        Aircraft_param.(Part).Planform.NBeams = numelements; % FIX
        Aircraft_param.(Part).Planform.CT_y_eta = 0.08; % Chris - The carry-through information needs to exist in the framework somewhere
        if ~isempty(TargetCl)
            Aircraft_param.(Part).Planform.Cl = TargetCl;
        end
    end
    
    Aircraft_param.(Part).Planform.y = LSObj(i).Eta_';
    
    Aircraft_param.(Part).Planform.CS.exist = 1;
    
    % Initialise the thickness
    Aircraft_param.(Part).Box.skin.eta       = [0,1];
    Aircraft_param.(Part).Box.skin.thickness = [0.03,0.03];
    Aircraft_param.(Part).Box.skin.max       = 0.15;
    Aircraft_param.(Part).Box.skin.min       = 0.003;
    
    Aircraft_param.(Part).Box.spar.eta       = [0,1];
    Aircraft_param.(Part).Box.spar.thickness = [0.01,0.01];
    Aircraft_param.(Part).Box.spar.max       = 0.15;
    Aircraft_param.(Part).Box.spar.min       = 0.003;
    
    Aircraft_param.(Part).Box.stringer.eta   = [0,1];
    Aircraft_param.(Part).Box.stringer.area  = [0.0001,0.0001];
    Aircraft_param.(Part).Box.stringer.max   = 0.05;
    Aircraft_param.(Part).Box.stringer.min   = 0.00005;
    
    tempeta = vertcat(CS.Eta);
    tempLE  = vertcat(CS.xLE);
    
    Aircraft_param.(Part).Planform.CS.inboard   = tempeta(:,1:end-1);
    Aircraft_param.(Part).Planform.CS.outboard  = tempeta(:,2:end);
    Aircraft_param.(Part).Planform.CS.chord     = 1-tempLE(:,1:end-1);
    Aircraft_param.(Part).Planform.CS.chord_le  = tempLE(:,1:end-1);
    Aircraft_param.(Part).Planform.CS.name      = {CS.Name};
    
    if isWing
        ControlSurf = {CS.Name};
    end

    Aircraft_param.(Part).Box.smax = 2.80*10^8;
    Aircraft_param.(Part).Box.shmax = 2.80*10^8;
    Aircraft_param.(Part).Box.nu = 0.33;
    Aircraft_param.(Part).Box.emax = 0;
    
    if isWing
        Aircraft_param.(Part).Box.Rho   = 1.6*10^3;
        Aircraft_param.(Part).Box.E     = 7.8437*10^10;
        Aircraft_param.(Part).Box.nu    = 0.4978;
        Aircraft_param.(Part).Box.emax  = 5000/1.5;
        Aircraft_param.(Part).Box.smax  = Aircraft_param.(Part).Box.emax*Aircraft_param.(Part).Box.E/(10^6);
        Aircraft_param.(Part).Box.shmax = 1.8019e+08; % Shear stress 4*10^7
    else
        Aircraft_param.(Part).Box.Rho   = 1700;
        Aircraft_param.(Part).Box.E     = 6.89*10^15;
    end
    
    Aircraft_param.(Part).Box.KAs = 0.36;
    Aircraft_param.(Part).Box.ds = 0.12;
    Aircraft_param.(Part).Box.SP = 0.15;
    Aircraft_param.(Part).Box.RP = 0.6;
    
end

%% BLUFF BODY PROPERTIES
PartNames = arrayfun(@(x) {x.Parent.Name},BBObj);
portidx   = cellfun(@isempty,regexpi(PartNames,'port'));

BBObj = BBObj(portidx);
for i = 1:numel(BBObj)
    
    Part = BBObj(i).Name;
    Parent = BBObj(i).Parent.Name;
    isEngine = ~isempty(regexpi(Part,'engine'));
    
    if isEngine
        Part = 'Engine';
    end
    
    Aircraft_param.(Part).LiftingSurface = 0;
    
    isFuse = strcmpi(Part,'fuselage');
    isWing = ~isempty(regexpi(Parent,'wing','ONCE'));
    
    if isWing
        ParentName = 'Wing';
    else
        ParentName = BBObj(i).Parent.Name;
    end
    
    if isFuse == 1
        Aircraft_param.(Part).Parent         = [];
    else
        Aircraft_param.(Part).Parent         = ParentName;
    end
    
    Aircraft_param.(Part).Optimisation   = 0;
    
    % Fuselage PLANFORM PROPERTIES
    Aircraft_param.(Part).Planform.Length = BBObj(i).Length;
    
    Aircraft_param.(Part).Planform.y_eta  = BBObj(i).Eta_'; 
    
    Aircraft_param.(Part).Planform.radius = BBObj(i).Radius_';
    
    Aircraft_param.(Part).Planform.Offset = [];
    
    if regexpi(Part,'engine','ONCE')
        Aircraft_param.(Part).Properties.Mass = 3000;
        Aircraft_param.(Part).Planform.y_eta_parent = BBObj(i).SOffset;
    else
        Aircraft_param.(Part).Properties.Mass = 26574;
    end
    
    Aircraft_param.(Part).Properties.CG   = BBObj(i).AbsPosition + [BBObj(i).Length/2,0,0];
    
end


box_type = 3;

%% PAYLOAD PROPERTIES
Aircraft_param.Payload.LiftingSurface = 0;
Aircraft_param.Payload.Parent = 'Fuselage';

Aircraft_param.Payload.Planform.y_eta  = [0;1];
Aircraft_param.Payload.Planform.radius = 0.75*max(Aircraft_param.Fuselage.Planform.radius)*[1;1];
Aircraft_param.Payload.Planform.Offset = [];
%Aircraft_param.Payload.Planform.GlobalOffset = [16.4,0,0];
Aircraft_param.Payload.Planform.Length = Aircraft_param.Fuselage.Planform.Length/2;
Aircraft_param.Payload.Planform.NSec    = 1;
%param.Payload.Properties.Mass   = 16059/2;
Aircraft_param.Payload.Properties.Mass   = 2.92e+04 / 2;
Aircraft_param.Payload.Properties.CG     = Aircraft_param.Fuselage.Properties.CG;
Aircraft_param.Payload.Planform.y_eta_parent = 0.0;
Aircraft_param.Payload.Planform.Colour = [0,1,0];

%% WEIGHT PROPERTIES
Aircraft_param.Weights.MTOW       = ac.MTOM;
Aircraft_param.Weights.OEW        = ac.OEM;
Aircraft_param.Weights.FB         = 14800;
Aircraft_param.Weights.Payload    = 2.92e+04;

[beam_model,Aircraft_param,Optim_Variables] = setupbox(Aircraft_param,box_type);

% Extract the loadcases from the lc object
field = {'Label','Type','M','Alt','URDD3','Fuel','Payload','PayloadCG','AC_Mass'};
prop = {'Name','Type','Mach','Altitude','LoadFactor','FuelFraction','PayloadMass','CgMac','AcMass'};

LoadCases = cell2struct(get(lc,prop),field,2);

for i = 1:numel(LoadCases)
    % Convert to acceleration
    LoadCases(i).URDD3 = 9.81*LoadCases(i).URDD3;
    % Find the Control Surface names
    ControlNames     = table2cell(lc(i).CsDeflections(1:numel(lc(i).CsDeflection),1))';
    % Remove the flaps and the slats TEMPORARY FIX
    idx = ~cellfun(@isempty,regexp(ControlNames,('la')));
    ControlNames(idx) = [];
    % Modify the control names so that it adds the string '1r' to the end
    modControlNames  = cellfun(@(x) [x,'1r'],ControlNames,'UniformOutput',false);
    % Downselected CSDelfection
    modCsDeflection = lc(i).CsDeflection(~idx);
    % Specify a logical array so that we can identify the control surfaces
    % that will be fixed
    ControlIdx = logical(zeros(1,numel(modControlNames)));
    % Find the control surfaces that belong to the parts being optimimise
    for j = 1:numel(ControlSurf)
        ControlIdx = or(ControlIdx,~cellfun(@isempty,regexp(modControlNames,['^' ControlSurf{j}])));
    end
    % Assign the values and labels
    LoadCases(i).CS.Label = modControlNames(ControlIdx);
    LoadCases(i).CS.Value = modCsDeflection(ControlIdx);
end
end