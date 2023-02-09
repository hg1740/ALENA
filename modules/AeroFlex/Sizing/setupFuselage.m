% The following function intialises the sizing variables and initialises
% the vortex lattice for a given type of box section

function  [DATA,BAERO] = setupFuselage(DATA,BAERO,PLANFORM)

nbaero = DATA.Info.nbaero + 1;

BAERO.ID(nbaero)        = nbaero;
BAERO.geo.ref_point(nbaero,:) = PLANFORM.ref_point;
BAERO.CP(nbaero)        = 0;
BAERO.geo.L(nbaero)     = PLANFORM.Length;
BAERO.geo.Nelem(nbaero) = numel(PLANFORM.y_eta);

if isfield(PLANFORM,'Colour')
    BAERO.Colour(nbaero,:) = PLANFORM.Colour;
else
    BAERO.Colour(nbaero,:) = [0.5,0.5,0.5];
end

BAERO.SET(nbaero)       = nbaero;
BAERO.geo.fs{nbaero}    = PLANFORM.y_eta;

% Do I need to define these in the global frame?
% BAERO.geo.fsz{nbaero} = PLANFORM.z;

BAERO.geo.Rs{nbaero} = PLANFORM.radius;

BAERO = body_lattice_setupV2(BAERO, []);

DATA.Info.nbaero = nbaero;
end