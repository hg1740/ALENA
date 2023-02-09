classdef (ConstructOnLoad) GustResult < awi.model.ResultSet
    %
    % Represents the result of performing a Loac Case analysis on a Aircraft
    
    properties % (SetAccess = immutable)
        
        %Extend metadata with content specific to awi
        LoadCase;
        RotMat;
        ControlDeflection;
        Loads;
        GustLength;
        
    end
    
    properties (Dependent)
        
        Fx
        Fy
        Fz
        Mx
        My
        Mz
        Eta
       
        Legend
    end
    
    methods %get / set
        
        function val = get.Legend(obj)
            
            %Columns of data are stored by loadcase then by linear / non-linear
            val = arrayfun(@(x) [obj.Name,' (',num2str(x),'m)'],obj.GustLength,'UniformOutput',false);
            
        end
        
        function val = get.Eta(obj)
            
            val = linspace(0,1,numel(obj.Loads{1}.x_f{1}(1:6:end,1)))';
            val = [val;flipud(val)];
            
        end
        
        function val = get.Fx(obj)
            %,max(obj.Results(i).Loads{ii}.x_f{1}(1:6:end,:),[],2)
            nloads = numel(obj.Loads);
            valmin = [];
            valmax = [];
            for i = 1:nloads
                valmin = [valmin,min(obj.Loads{i}.x_f{1}(1:6:end,:),[],2)];
                valmax = [valmax,flipud(max(obj.Loads{i}.x_f{1}(1:6:end,:),[],2))];
            end
            
            val = [valmin;valmax];
            
        end
        
        function val = get.Fy(obj)
            
            nloads = numel(obj.Loads);
            valmin = [];
            valmax = [];
            for i = 1:nloads
                valmin = [valmin,min(obj.Loads{i}.x_f{1}(2:6:end,:),[],2)];
                valmax = [valmax,flipud(max(obj.Loads{i}.x_f{1}(2:6:end,:),[],2))];
            end
            
            val = [valmin;valmax];
            
        end
        
        function val = get.Fz(obj)
            
            nloads = numel(obj.Loads);
            valmin = [];
            valmax = [];
            for i = 1:nloads
                valmin = [valmin,min(obj.Loads{i}.x_f{1}(3:6:end,:),[],2)];
                valmax = [valmax,flipud(max(obj.Loads{i}.x_f{1}(3:6:end,:),[],2))];
            end
            
            val = [valmin;valmax];
            
        end
        
        function val = get.Mx(obj)
            
            nloads = numel(obj.Loads);
            valmin = [];
            valmax = [];
            for i = 1:nloads
                valmin = [valmin,min(obj.Loads{i}.x_f{1}(4:6:end,:),[],2)];
                valmax = [valmax,flipud(max(obj.Loads{i}.x_f{1}(4:6:end,:),[],2))];
            end
            
            val = [valmin;valmax];
            
        end
        
        function val = get.My(obj)
            
            nloads = numel(obj.Loads);
            valmin = [];
            valmax = [];
            for i = 1:nloads
                valmin = [valmin,min(obj.Loads{i}.x_f{1}(5:6:end,:),[],2)];
                valmax = [valmax,flipud(max(obj.Loads{i}.x_f{1}(5:6:end,:),[],2))];
            end
            
            val = [valmin;valmax];
            
        end
        
        function val = get.Mz(obj)
            
            nloads = numel(obj.Loads);
            valmin = [];
            valmax = [];
            for i = 1:nloads
                valmin = [valmin,min(obj.Loads{i}.x_f{1}(6:6:end,:),[],2)];
                valmax = [valmax,flipud(max(obj.Loads{i}.x_f{1}(6:6:end,:),[],2))];
            end
            
            val = [valmin;valmax];
            
        end
        
    end
    
    methods % construction / destruction
        
        function obj = GustResult(varargin)
            %
            % Construct a TrimResult object OBJ, based on the Trim Analysis PAR
            
            %Pass it on
            obj@awi.model.ResultSet(varargin{:});
            
            %Extend property groups
            obj.addPropertyGroup('General', ...
                'LoadCase', 'LoadCase');
            
        end
        
    end
    
end
