%BEAMMODEL object containing FEM elements for NASTRAN and NeoCass.
% This object is used as a container to store all the elements that are
% needed for an FE model. Data is extracted from FAME and then merged into
% this object.
%
% All the element definitions follow the NASTRAN format, hence this manual
% can be consulted for further information.
classdef BeamModel < matlab.mixin.SetGet
    
    properties
        Caero    = [];      % Defines an aerodynamic macro element (panel)
        Cbar     = [];      % Defines a simple beam element
        Conm2    = [];      % Defines a concentrated mass at a grid point
        Force    = [];      % Defines a static concentrated force at a grid point by specifying a vector
        Grid     = [];      % Defines the location of a geometric grid point, the directions of its displacement, and its permanent single-point constraints
        Info     = [];      % Summary of total amount of elements
        Mat      = [];      % Defines the material properties
        Moment   = [];      % Defines a static concentrated moment at a grid point by specifying a vector
        Pbar     = [];      % Defines the properties of a simple beam element (CBAR entry)
        Rbe0     = [];      % Provides the dependency information for each aeronode to its master
        Rbe2     = [];      % Defines a rigid body with independent degrees-of-freedom that are specified at a single grid point and with dependent degrees-of-freedom that are specified at an arbitrary number of grid points
        Sets     = [];      % Used to defined an interpolation set of nodes which is used by one of the available spatial coupling method to transfer data between structural and aerodynamic meshes
        Spc      = [];      % Defines a set of single-point constraints and enforced motion
        Spline   = [];      % Defines a surface spline for interpolating motion and/or forces for aeroelastic problems on aerodynamic geometries defined by regular arrays of aerodynamic points
        Thrust   = [];      % Thrust follower force used with the non-linear structural solvers
        Loadcase = [];      % Definition of difefrent load cases
        PartId   = [];      % Lists the id numbers for different parts and types
    end
    
    methods
        function obj = BeamModel(~)
            
        end
              
        function h = plotbar(obj,fig)
            
            if nargin == 2
                figure(fig);
                hold on;
            else
                figure;
            end
            
            h = hggroup;
            
            for i = 1:numel(obj.Cbar)
                iStart  = obj.Cbar(i).conn(1) == [obj.Grid.id];
                startPt = obj.Grid(iStart).coord;
                iEnd    = obj.Cbar(i).conn(2) == [obj.Grid.id];
                endPt   = obj.Grid(iEnd).coord;
                coords  = [startPt;endPt];
                
                plot3(coords(:,1), coords(:,2), coords(:,3),'c-','LineWidth',2,'Parent',h);
                hold on;
            end
            axis equal
            
        end
        
        function h = plotrbe0(obj,fig)
            
            if nargin == 2
                figure(fig);
                hold on;
            else
                figure;
            end
            
            h = hggroup;
            
            for i = 1:numel(obj.Rbe0)
                iStart  = obj.Rbe0(i).grid == [obj.Grid.id];
                startPt = obj.Grid(iStart).coord;
                iEnd    = obj.Rbe0(i).data(1) == [obj.Grid.id];
                endPt   = obj.Grid(iEnd).coord;
                coords  = [startPt;endPt];
                
                plot3(coords(:,1), coords(:,2), coords(:,3),'b-','Parent',h);
                hold on;
                iEnd    = obj.Rbe0(i).data(2) == [obj.Grid.id];
                endPt   = obj.Grid(iEnd).coord;
                coords  = [startPt;endPt];
                
                plot3(coords(:,1), coords(:,2), coords(:,3),'b-','Parent',h);
            end
        end
        function h = plotrbe2(obj,fig)
            
            if nargin == 2
                figure(fig);
                hold on;
            else
                figure;
            end
            
            h = hggroup;
            
            for i = 1:numel(obj.Rbe2)
                if ~isempty(obj.Rbe2(i).gn)
                    iStart  = obj.Rbe2(i).gn == [obj.Grid.id];
                    startPt = obj.Grid(iStart).coord;
                    iEnd    = obj.Rbe2(i).gm == [obj.Grid.id];
                    endPt   = obj.Grid(iEnd).coord;
                    coords  = [startPt;endPt];
                    
                    plot3(coords(:,1), coords(:,2), coords(:,3),'k-','Parent',h);
                end
            end
            
        end
        
        function h = plotnodes(obj,fig)
           if nargin == 2
                figure(fig);
                hold on;
            else
                figure;
           end
            
           
             
        end
    end
    
end

