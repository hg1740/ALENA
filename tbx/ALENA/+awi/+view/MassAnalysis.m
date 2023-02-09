classdef MassAnalysis < awi.view.Analysis
    
    properties (Access = protected)
        
        %The table in which we display relevant details
        hTable;
        
    end
    
    methods % constructor
        
        function obj = MassAnalysis( model, varargin )
            
            %Call superclass constructor - do not pass the model in, better to wait
            obj@awi.view.Analysis( model, ...              'Title', 'Mass Analysis', ...
                'OfferLoadCaseSelection', false, ...
                varargin{:} );
            
            %Initialise view
            obj.hTable = uitable('Parent', obj.hPanel, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);            
            
            %Now store the model - which will in turn trigger an update
            obj.Model = model;
            
        end
        
    end
    
    methods ( Access = protected )
        
        function update(obj)
            
            %Start with base class
            update@awi.view.Analysis(obj);
            
            %No Aircraft ?
            if isempty(obj.Aircraft)
                
                %Nothing to show
                dat = {};
                row = {};
                
            else
                
                %Get details
                TotalMass = obj.Aircraft.TotalMass;
                
                %But there is more than one measure of mass
                TakeOffMass = obj.Aircraft.MTOM;
                MaxZeroFuelMass = obj.Aircraft.MZFM;
                MaxLandingWeight = obj.Aircraft.MLM;
                
                %Tabulate
                dat = {TotalMass; TakeOffMass; MaxZeroFuelMass; MaxLandingWeight};
                row = {'Total Mass'; 'Take Off Mass'; 'Max Zero Fuel Mass'; 'Max Landing Weight'};
                
            end
                           
            %Display in table
            set(obj.hTable, 'Data', dat, 'RowName', row, 'ColumnName', {'Value'});
            
        end
        
    end
    
end
