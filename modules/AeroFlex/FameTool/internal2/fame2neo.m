%% FAME2NEO Convert FAME model object into NEOCASS format file
%
% This function creates NEOCASS input decks from the FAME model object.
% 
% FAME2NEO(FILENAME,FAME2MATOBJECT)
%
%       FILENAME : Hard coded as |StructureAerodynamicFile.dat| for now. 
%                  TODO : A file name option needs to be created    
% FAME2MATOBJECT : Fame model object
%
%
% See NEOCASS and NASTRAN manuals for further explanations of the decks.
%
%   Restrictions:
%       * Assume a right-handed axis system from FAME where:
%                 X - from nose to tail
%                 Y - right wing
%                 Z - up
%
%
%   See also |writeAdamsFiles| method within Fame2mat object.

%   Copyright 2016 University of Bristol
%   Private function.
function fame2neo(file,obj)

% Check which properties belong to the VTP and ignore

fp = fopen(file, 'w');

fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
fprintf(fp,'$       MATERIAL PROPERTIES \n');
fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');

for i = 1:numel(obj.Mdl.Mat)
    str = genNeocassDeck(obj.Mdl.Mat(i));
    
    for j = 1:numel(str)
        fprintf(fp,'%s\n',str{j});
    end
end

fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
fprintf(fp,'$       GRID NODES \n');
fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');

for i = 1:numel(obj.Mdl.Grid)
    str = genNeocassDeck(obj.Mdl.Grid(i));
    
    for j = 1:numel(str)
        fprintf(fp,'%s\n',str{j});
    end
end
fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
fprintf(fp,'$       BAR PROPERTIES \n');
fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');

for i = 1:numel(obj.Mdl.Pbar)
    str = genNeocassDeck(obj.Mdl.Pbar(i));
    
    for j = 1:numel(str)
        fprintf(fp,'%s\n',str{j});
    end
end

fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
fprintf(fp,'$       BAR ELEMENTS \n');
fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');

for i = 1:numel(obj.Mdl.Cbar)
    str = genNeocassDeck(obj.Mdl.Cbar(i));
    
    for j = 1:numel(str)
        fprintf(fp,'%s\n',str{j});
    end
end

fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
fprintf(fp,'$       AERO BEAMS \n');
fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');

for i = 1:numel(obj.Mdl.Rbe0)
    str = genNeocassDeck(obj.Mdl.Rbe0(i));
    
    for j = 1:numel(str)
        fprintf(fp,'%s\n',str{j});
    end
end

fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
fprintf(fp,'$       AERO CARDS \n');
fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');


for i = 1:numel(obj.Mdl.Caero)
    str = genNeocassDeck(obj.Mdl.Caero(i));
    
    for j = 1:numel(str)
        fprintf(fp,'%s\n',str{j});
    end
end

fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
fprintf(fp,'$       LINK AEROELASTIC VARIABLES \n');
fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');

cslabels= cell(0);
count = 0;
caerocount = numel(obj.Mdl.Caero);
for i = 1:caerocount
    if ~isempty(obj.Mdl.Caero(1,i).csType)
        count = count +1;
        cslabels{count,1} = obj.Mdl.Caero(1,i).csType;
    end
end

uni_cs = unique(cslabels);
aecount = 0;

for j = 1:numel(uni_cs)
    
    idx = find(strcmp({obj.Mdl.Caero.csType},uni_cs(j)));
    count = numel(idx);
    
    for i = 1:count-1
        id   = sprintf('%-8d', i + aecount);
        ae1  = cat(2,obj.Mdl.Caero(idx(1)).csId,blanks(8));
        ae2  = cat(2,obj.Mdl.Caero(idx(1 + i)).csId,blanks(8));
        
        link = sprintf('%-8d', obj.Mdl.Caero(idx(1 + i)).csLink);
        
        fprintf(fp,'AELINK  %s%s%s%s\n', id, ae2(1:8), ae1(1:8), link);
    end
    aecount = aecount + i;
end


fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
fprintf(fp,'$       INTERPOLATION SETS \n');
fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');

for i = 1:numel(obj.Mdl.Sets)
    
    str = genNeocassDeck(obj.Mdl.Sets(i));
    
    for j = 1:numel(str)
        fprintf(fp,'%s\n',str{j});
    end
end

fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
fprintf(fp,'$       SPLINE SETS \n');
fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');

for i = 1:numel(obj.Mdl.Spline)
    
    str = genNeocassDeck(obj.Mdl.Spline(i));
    
    for j = 1:numel(str)
        fprintf(fp,'%s\n',str{j});
    end
    
end

% fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
% fprintf(fp,'$  SUPORT CARD \n');
% fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');

% SUPORT 1  1005 123456

fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
fprintf(fp,'$       CONM CARDS \n');
fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');

for i = 1:numel(obj.Mdl.Conm2)
    
    if ~any(strcmp(obj.Mdl.Conm2(i).type,{'Structural','Engine','Secondary'}))
        continue
    end
    
    str = genNeocassDeck(obj.Mdl.Conm2(i));
    
    for j = 1:numel(str)
        fprintf(fp,'%s\n',str{j});
    end
end

fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
fprintf(fp,'$       RBE2 CARDS \n');
fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');

for i = 1:numel(obj.Mdl.Rbe2)
    str = genNeocassDeck(obj.Mdl.Rbe2(i));
    
    for j = 1:numel(str)
        fprintf(fp,'%s\n\n',str{j});
    end
end

fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
fprintf(fp,'$       PARTID CARDS \n');
fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');

for i = 1:numel(obj.Mdl.PartId)
    
    str = genNeocassDeck(obj.Mdl.PartId(i));
    
    for j = 1:numel(str)
        fprintf(fp,'%s\n',str{j});
    end
end

fclose all;
end
