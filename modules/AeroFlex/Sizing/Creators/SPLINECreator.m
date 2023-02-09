function [CAERO,INFO,IDs] = SPLINECreator(CAERO,CAEROpart,Bar_prop,INFO,startID,setID)

ninterp = INFO.ninterp;
INFO.spline_type = 2;
panelcount = 1;
for i = 1:length(CAEROpart.geo.ny)
    ninterp = ninterp +1;
    CAERO.Interp.Type(ninterp) = 2;
    CAERO.Interp.ID(ninterp) = CAEROpart.ID(i);
    CAERO.Interp.Patch(ninterp) = CAEROpart.ID(i);
    
    subpanelcount = CAEROpart.geo.ny(i)*(CAEROpart.geo.nx(i)+CAEROpart.geo.fnx(i));
    CAERO.Interp.Index(ninterp,1) = 1;
    CAERO.Interp.Index(ninterp,2) = subpanelcount;
    CAERO.Interp.Set(ninterp) = setID;
    CAERO.Interp.Param(ninterp, 1) = 1;
    CAERO.Interp.Param(ninterp, 2) = 2;
    CAERO.Interp.Param(ninterp, 3) = 0;
    CAERO.Interp.Param(ninterp, 4) = 1;
    CAERO.Interp.Param(ninterp, 5) = 1e14;
    panelcount = panelcount + 1;
end
IDs = CAERO.Interp.ID(INFO.ninterp + 1 : INFO.ninterp + length(CAEROpart.geo.ny));
INFO.ninterp = ninterp;
end