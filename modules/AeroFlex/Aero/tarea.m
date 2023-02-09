function [panel_area]=tarea(XYZ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tarea: Subsidary function for TORNADO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the area of each panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Tomas Melin, KTH, Department of Aeronautics	%
%				Copyright 2000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Subsidaty function for TORNADO
% Called by:	coeff_create
%
% Calls:			MATLAB 5.2 std fcns
% Loads:	none
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a b c]=size(XYZ);
for i=1:a
    p1=[XYZ(i,1,1) XYZ(i,1,2) XYZ(i,1,3)];	%sets up the vectors
    p2=[XYZ(i,2,1) XYZ(i,2,2) XYZ(i,2,3)];	%to the corners of the
    p3=[XYZ(i,3,1) XYZ(i,3,2) XYZ(i,3,3)];	%panel.
    p4=[XYZ(i,4,1) XYZ(i,4,2) XYZ(i,4,3)];
    
    a=p2-p1;	%sets up the edge vectors
    b=p4-p1;
    c=p2-p3;
    d=p4-p3;
    
    ar1=norm(cross(b,a))/2;	%claculates the ctoss product of
    ar2=norm(cross(c,d))/2;	%two diagonal corners
    
    panel_area(i)=ar1+ar2;	%Sums up the product to make the
end						    %Area
end% function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%