%% DEV. FWT
%
% TODO - Allow definition of the fold angle
% TODO - Check Nastran normal modes for a hinged wing with a tip mass
%  

clear all
close all

LS = awi.model.LiftingSurface.makeHodgesWing;
%LS = awi.model.LiftingSurface.makeA320Wing;
% draw(LS);
% FEM = convertToFE(LS);
% draw(FEM);

FWT = insertWingFold(LS, 'FlareAngle', 20);

draw(LS);