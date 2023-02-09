classdef (ConstructOnLoad) SizeResult < awi.model.ResultSet
    %
    % Represents the result of performing a Trim analysis on a Aircraft

    properties % (SetAccess = immutable)

        %Extend metadata with content specific to awi
        LoadCase;
        BeamModel;
        OptimVariables;
        AnalysisType;

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
            val = cellfun(@(y)cellfun(@(x)[x, ' (', y, ')'], obj.LoadCase, 'UniformOutput', false), ...
                obj.AnalysisType, 'UniformOutput', false);
            val = horzcat(val{:});

        end

        function val = get.Eta(obj)
            
            val = arrayfun(@(x)obj.OptimVariables{x}.Wing.y_mbox/max(obj.OptimVariables{x}.Wing.y_mbox), 1:numel(obj.OptimVariables), 'UniformOutput', false);
            val = horzcat(val{:});

            %But (for time being at least) there should be one and only one Eta
            if size(val,2) > 1
                assert(all(diff(val, [], 2) == 0), 'inconsistent Eta values');
            end          

            %In which case
            val = val(:,1);

        end

        function val = get.Fx(obj)

            val = arrayfun(@(x)obj.OptimVariables{x}.Wing.Fx, 1:numel(obj.OptimVariables), 'UniformOutput', false);
            val = horzcat(val{:});

        end

        function val = get.Fy(obj)

            val = arrayfun(@(x)obj.OptimVariables{x}.Wing.Fy, 1:numel(obj.OptimVariables), 'UniformOutput', false);
            val = horzcat(val{:});

        end

        function val = get.Fz(obj)
            
            val = arrayfun(@(x)obj.OptimVariables{x}.Wing.Fz, 1:numel(obj.OptimVariables), 'UniformOutput', false);
            val = horzcat(val{:});

        end

        function val = get.Mx(obj)

            val = arrayfun(@(x)obj.OptimVariables{x}.Wing.Mx, 1:numel(obj.OptimVariables), 'UniformOutput', false);
            val = horzcat(val{:});

        end

        function val = get.My(obj)

            val = arrayfun(@(x)obj.OptimVariables{x}.Wing.My, 1:numel(obj.OptimVariables), 'UniformOutput', false);
            val = horzcat(val{:});

        end

        function val = get.Mz(obj)

            val = arrayfun(@(x)obj.OptimVariables{x}.Wing.Mz, 1:numel(obj.OptimVariables), 'UniformOutput', false);
            val = horzcat(val{:});

        end

    end

    methods % construction / destruction

        function obj = SizeResult(varargin)
            %
            % Construct a SizeResult object OBJ, based on the Trim Analysis PAR

            %Pass it on
            obj@awi.model.ResultSet(varargin{:});

            %Extend property groups
            obj.addPropertyGroup('General', ...
                'LoadCase', 'LoadCase');

        end

    end

end