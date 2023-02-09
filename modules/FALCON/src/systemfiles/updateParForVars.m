function [x0,Matrices,Aero,Sim] = updateParForVars(x0,Matrices,Aero,Sim,ParForVars)

if isfield(ParForVars,'x0')
    f = fieldnames(ParForVars.x0);
    for ii = 1:length(f)
        x0.(f{ii}) = ParForVars.x0.(f{ii});
    end
end

if isfield(ParForVars,'Matrices')
    f = fieldnames(ParForVars.Matrices);
    for ii = 1:length(f)
        Matrices.(f{ii}) = ParForVars.Matrices.(f{ii});
    end
end

if isfield(ParForVars.Aero,'gust')
    f = fieldnames(ParForVars.Aero.gust);
    for ii = 1:length(f)
        Aero.gust.(f{ii}) = ParForVars.Aero.gust.(f{ii});
    end
end

if isfield(ParForVars,'Sim')
    f = fieldnames(ParForVars.Sim);
    for ii = 1:length(f)
        Sim.(f{ii}) = ParForVars.Sim.(f{ii});
    end
end