classdef TestOMEGA < matlab.unittest.TestCase
    %TestOmega Contains tests for the Object Orientated Model Generaton for
    %Aircraft (OMEGA) module.
    %
    % TODO - Add test for each of the mixin classes.
    % TODO - Test different build methods and make sure the same geometry
    % is returned each time.  
    % TODO - Write a unit test which invokes each fe object and tries to
    % draw, write to a file, etc.

    properties (TestParameter)
        ObjectsWithDefaultBuilds = { ...
            'awi.model.LiftingSurface'};
        DrawableClasses = getDrawableClasses
    end
    
    properties (SetAccess = private)
        TestFigure
        TestAxes
    end
            
    methods(TestMethodSetup)
        function setup_figure(obj)
            %setup_figure Makes a MATLAB figure with a single axes for
            %plotting data during tests.
            
            obj.TestFigure = figure('Name', 'Test Figure');
            obj.TestAxes   = axes('Parent', obj.TestFigure, 'NextPlot', 'add');            
        end
    end
    
    methods (Test) %Generating geometries
        function geom_make_default_object(obj, ObjectsWithDefaultBuilds)
            %geom_make_default_object Makes each default build object and
            %draws it.
            
            [build_names, ModelObj] = obj.get_default_build_names(ObjectsWithDefaultBuilds);
            
            for ii = 1 : numel(build_names)
                NewModelObj = runDefaultBuild(ModelObj, build_names{ii});
                draw(NewModelObj, obj.TestAxes);
                cla(obj.TestAxes);
            end   
            
        end
    end
    
    methods (Test) %Generating FEM
        function fe_convert_defaults_to_FE(obj, ObjectsWithDefaultBuilds)
            %fe_convert_defaults_to_FE Makes each type of default object
            %and tries to convert it to a FE model which is then plotted.
            
           [build_names, ModelObj] = obj.get_default_build_names(ObjectsWithDefaultBuilds);
            
            for ii = 1 : numel(build_names)
                NewModelObj = runDefaultBuild(ModelObj, build_names{ii});
                if ~NewModelObj.CanConvertToFE
                    continue
                end
                FEM = convertToFE(NewModelObj);
                draw(FEM, obj.TestAxes);
                cla(obj.TestAxes);
            end   
            
        end
        function fe_mesh_lifting_surface_GFEM(obj)
            %fe_mesh_lifting_surface_GFEM Makes a default LiftingSurface
            %object and then converts it to a GFEM FEM.
            
            return
            
            etaRib    = linspace(0, 1, 10);
            lambdaRib = 40 * rand(size(etaRib));
            
            LS = awi.model.LiftingSurface.makeHodgesWing; 
            
            testRibDefinition(LS, etaRib, lambdaRib);            
        end
    end
    
    methods (Test) %Drawing geoemtries
        function drawEmptyClass(obj, DrawableClasses)
            %drawEmptyClass Checks that a drawable class can be
            %instantiated and drawn without any error.
            
            fcn = str2func(DrawableClasses);
            DrawableObj = fcn();
            
            if ~isa(DrawableObj, 'mvc.mixin.Drawable')
                return
            end
            
            drawElement(DrawableObj, obj.TestAxes);
            
        end
    end
    
    methods (TestMethodTeardown)
        function close_figures(obj)
            %close_figures Closes all figures which were opened during the
            %tests.
            
            close(obj.TestFigure)
            
        end
    end
    
    methods (Static)
        function [build_names, ModelObj] = get_default_build_names(cls_name)
            %get_default_build_names
            
            assert(exist(cls_name, 'class') == 8, sprintf([ ...
                'Expected ''cls_name'' to be a valid class name. ', ...
                '%s is not a class in the ALENA framework'], cls_name));
            
            func = str2func(cls_name);
            
            ModelObj    = func();
            build_names = ModelObj.AvailableDefaults;
            
        end
    end
    
end

function modelClassnames = getDrawableClasses
%getDrawableClassnames Returns the classnames of any class in the awi.model
%package folder that inherits from mvc.mixin.Drawable

classFolder = fullfile( ...
    fileparts(fileparts(mfilename('fullpath'))), ...
    'tbx', 'ALENA','+awi', '+model');
Contents = dir(classFolder);

name     = strrep({Contents(~[Contents.isdir]).name}, '.m', '');
idxClass = cellfun(@(x) exist(['awi.model.', x], 'class') == 8, name);

modelClassnames = strcat('awi.model.', name(idxClass));

end