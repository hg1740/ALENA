classdef Cabin < awi.model.Compartment 
    %Cabin Defines a cabin layout uses different parameterisation schemes.
    %
    % TODO - Add build methods and parameter sets.
    
    %Cabin layout properties
    properties
        %Distance between rows
        SeatPitch
        %Width of a single seat
        SeatWidth
        %Height of a single seat
        SeatHeight
        %Width of an aisle
        AisleWidth
        %Key describing the layout for a single row
        RowArrangement
        %Key describing the layout of the cabin
        SeatLayout
    end
    
    properties
        %Total number of seats in the cabin
        NumSeats
        %Space between the edge of the first seat and the edge of the cabin
        EmptyWidth
    end
    
    properties
        %Position of seats w.r.t the forward left edge of the cabin
        SeatPositions
    end    
    
    %Internal properties
    properties (SetAccess = private, Hidden = true)
        %Normalised position of each seat along the length of the cabin
        CabinSeatEta
        %Distance of the start of the seat from the cabin floor [~12in]
        SeatFloatHeight = 0.3; 
        %Coordinates of the seats collected by row
        SeatPositionByRow
    end
    
    %Internal dependent properties
    properties (Dependent, Hidden = true) 
        %Normalised position of each seat along the length of the beam
        SeatEta
        %Coordinates of the patch used to draw the '3D' seat
        SeatPatchCoords
        %Coordinates of the seats collected by row (Beam local frame)
        SeatPosByRowInBeamFrame
        %Cabin datum position - Coordinates of front left corner
        CabinDatum
    end
    
    %Cabin dimensions 
    properties (Dependent)
        %Length of the cabin (in the s-domain)
        CabinLength
        %Width of the cabin at the points where a cross-section has been
        %defined
        CabinWidth
        %Coordinates of the centre of the seats in the parent beam
        %coordinate system. Located on the cabin floor at the moment
        SeatCoords
        %Number of seats per row
        SeatsPerRow
    end   
    
    methods % set / get
        function set.SeatPitch(obj, val)        %set.SeatPitch
            %set.SeatPitch Set method for the property 'SeatPitch'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty', ...
                'positive', 'finite', 'real'}, class(obj), 'SeatPitch');
            obj.SeatPitch = val;
            
        end
        function set.SeatWidth(obj, val)        %set.SeatWidth
            %set.SeatWidth Set method for the property 'SeatWidth'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty', ...
                'positive', 'finite', 'real'}, class(obj), 'SeatWidth');
            obj.SeatWidth = val;
            
        end
        function set.SeatHeight(obj, val)        %set.SeatHeight
            %set.SeatHeight Set method for the property 'SeatHeight'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty', ...
                'positive', 'finite', 'real'}, class(obj), 'SeatHeight');
            obj.SeatHeight = val;
            
        end
        function set.AisleWidth(obj, val)       %set.AisleWidth
            %set.AisleWidth Set method for the property 'AisleWidth'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty', ...
                'positive', 'finite', 'real'}, class(obj), 'AisleWidth');
            obj.AisleWidth = val;
            
        end
        function set.RowArrangement(obj, val)   %set.RowArrangement
            %set.RowArrangement Set method for the property
            %'RowArrangement'.
            
            validateattributes(val, {'char'}, {'row', 'nonempty'}, ...
                class(obj), 'RowArrangement');
            obj.RowArrangement = val;
            
        end
        function set.SeatLayout(obj, val)       %set.SeatLayout
            %set.SeatLayout Set method for the property 'SeatLayout'.
            
            validateattributes(val, {'cell'}, {'column', 'nonempty'}, ...
                class(obj), 'SeatLayout');
            assert(iscellstr(val), ['Expected the ''SeatLayout'' to ', ...
                'be a cell-array of characters.']);
            idxAS = cellfun(@(x) all(or(ismember(x, 'a'), ismember(x, 's'))), val);
            assert(all(idxAS), ['The ''SeatLayout'' must be defined ', ...
                'using only ''a'' (aisle) and ''s'' (seat) tokens.'])
            obj.SeatLayout = val;
            
        end
        function set.NumSeats(obj, val)         %set.NumSeats 
            %set.NumSeats Set method for the property 'SeatPitch'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty', ...
                'positive', 'finite', 'real'}, class(obj), 'NumSeats');
            obj.NumSeats = val;
            
        end
        function set.EmptyWidth(obj, val)       %set.EmptyWidth
            %set.EmptyWidth Set method for the property 'SeatPitch'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty', ...
                'nonnegative', 'finite', 'real'}, class(obj), 'EmptyWidth');
            obj.EmptyWidth = val;
            
        end
        function set.SeatPositions(obj, val)    %set.SeatPositions
            %set.SeatPositions Set method for the property 'SeatPositions'.
            
            validateattributes(val, {'numeric'}, {'2d', 'ncols', 3, ...
                'finite', 'real'}, class(obj), 'SeatPositions');
            obj.SeatPositions = val;
        end
        function val = get.CabinLength(obj)     %get.CabinLength
            %get.CabinLength Get methdod for the dependent property
            %'CabinLength'.
            
            val = [];
            
            if isempty(obj.BeamHandle) %Escape route
                return
            end
            
            %Get the (x,y,z) data of the beam
            [xd, yd, zd] = xyzdata(obj.BeamHandle);
            if isempty(xd) %Escape route
                return
            end
            r    = awi.model.Stick.getLineLength(xd, yd, zd);
            eta  = r ./ r(end);
            
            %Interpolate to find the s-wise positions of the Cabin
            rCabin = interp1(eta, r, obj.EtaLocations);
            
            %The length is the difference between the first and last
            val = rCabin(end) - rCabin(1);
            
        end
        function val = get.CabinWidth(obj)      %get.CabinWidth 
            %get.CabinWidth Get method for the property 'CabinWidth'.
            
            val = [];
            
            if isempty(obj.CompartmentCrossSections)
                return
            end
            
            %Width is just the difference between the max/min x-coords
            csX = vertcat(obj.CompartmentCrossSections.X);
            val = abs(max(csX, [], 2) - min(csX, [], 2));
            
        end
        function val = get.SeatCoords(obj)      %get.SeatCoords 
            %get.SeatCoords Get method for the dependent property
            %'SeatCoords'.
            
            val  = [];
            xyz0 = obj.CabinDatum; 
            
            if isempty(obj.SeatPositions) || isempty(xyz0) %Escape route
                return
            end
                       
%             %Get coordinates of parent beam
%             [xd, yd, zd] = xyzdata(obj.BeamHandle);
%             if isempty(xd) %Escape route
%                 return
%             end
%             r    = awi.model.Stick.getLineLength(xd, yd, zd);
%             eta  = r ./ r(end);
%             
%             %Interpolate to find coordinates at seat eta positions
%             se = obj.SeatEta;
%             if isempty(se) %Escape route
%                 return
%             end
%             xyzS = interp1(eta, [xd ; yd ; zd]', obj.SeatEta);
%                         
%             %Rotate coords through the local rotation matrix....
            
            %Translate into the local beam coordinate system
            val = obj.SeatPositions + xyz0;
            
        end
        function val = get.SeatsPerRow(obj)     %get.SeatsPerRow 
           %get.SeatsPerRow Get method for the dependent property
           %'SeatsPerRow'.
           
           val = [];
           
           if isempty(obj.SeatLayout) %Escape route
               return
           end
           
           val = cellfun(@(x) numel(strfind(x, 's')), obj.SeatLayout);
           
        end
        function val = get.SeatEta(obj)         %get.SeatEta 
            %get.SeatEta Get method for the dependent property 'SeatEta'.
            
            val = [];
            
            if isempty(obj.EtaLocations) || isempty(obj.CabinSeatEta)
                return
            end
            
            dEta = obj.EtaLocations(end) - obj.EtaLocations(1);            
            val  = obj.EtaLocations(1) + (obj.CabinSeatEta .* dEta);
            
        end
        function val = get.SeatPatchCoords(obj) %get.SeatPatchCoords 
           %get.SeatPatchCoords Get method for the seat 'patch' plot.
           
           val = [];
           
           if any(cellfun(@isempty, get(obj, {'SeatWidth', 'SeatHeight', ...
                   'SeatFloatHeight', 'SeatPitch'})))
               return
           end
           
           %Grab data
           sw   = obj.SeatWidth;
           sh   = obj.SeatHeight;
           sl   = obj.SeatPitch ./ 2; 
           dz   = obj.SeatFloatHeight;          
           sw_2 = sw / 2;
           sl_2 = sl / 2;
           maxH = sh + dz;
           
           %The chair will be represented as as 2 rectangles    
           val = [ ...
               -sl_2, -sw_2, dz   ; ...
               sl_2 , -sw_2, dz   ; ...
               sl_2 , -sw_2, maxH ; ...
               sl_2 , sw_2 , maxH ; ...
               sl_2 , sw_2 , dz   ; ...
               -sl_2, sw_2 , dz  ];
           
        end
        function val = get.SeatPosByRowInBeamFrame(obj) %get.SeatPosByRowInBeamFrame
            %get.SeatPosByRowInBeamFrame Get method for the dependent
            %property 'SeatPosByRowInBeamFrame'.
            
            val  = [];
            xyz0 = obj.CabinDatum;
            
            if isempty(obj.SeatPositionByRow) || isempty(xyz0) %Escape route
                return
            end
            
            SeatPosByRow = obj.SeatPositionByRow;
            
            val   = SeatPosByRow;
            val.X = val.X + xyz0(1);
            val.Y = val.Y + xyz0(2);
            val.Z = val.Z + xyz0(3);
            
        end
        function val = get.CabinDatum(obj)      %get.CabinDatum 
            %get.CabinDatum Get method for the dependent property
            %'CabinDatum'.
            
            val = [];
            
            if isempty(obj.CompartmentCrossSections)
                return
            end
            
            %Datum is the forward port position of the start of the cabin
            [xC, yC, zC] = calculateGlobalCoords( ...
                obj.CompartmentCrossSections(1));
            val = [xC(1), min(yC), min(zC)];
            
        end
    end
    
    methods % visualisation
        function hg = drawCabinLayout(obj, ht, varargin)
            %drawElement Draw method for the 'awi.model.Cabin' class.
            %
            % Shows the layout of the seats using a rudimentary 'patch'
            % representation of each seat.
                                   
            hg = [];
            
            %Grab coordinates of the seat
            xyzS = obj.SeatCoords;
            
            if isempty(xyzS) %Escape route
                return
            end
            
            %Draw the origin location of each seat
            hg = plot3(ht, xyzS(:, 1), xyzS(:, 2), xyzS(:, 3), ...
                'Marker', 'o', ...
                'MarkerFaceColor', 'r', ...
                'MarkerEdgeColor', 'k', ...
                'LineStyle'      , 'none', ...
                'Tag', 'Seat Origins');
            
            %Calculate the coordinates of the seat 'patch' objects            
            coords      = obj.SeatPatchCoords; %Relative to each seat origin!        
            if isempty(coords) %Escape route 
                return
            end
            nSeat       = size(xyzS, 1);
            nV          = 6;
            nPoint      = nSeat * nV;
            patchCoords = repmat(coords, [nSeat, 1]);
            vi = [[ ...
                1 : nV : nPoint   ; ...
                2 : nV : nPoint   ; ...
                5 : nV : nPoint   ; ...
                6 : nV : nPoint]' ; ...
                [ ...
                2 : nV : nPoint ; ...
                3 : nV : nPoint ; ...
                4 : nV : nPoint ; ...
                5 : nV : nPoint ]'];
            xyz0 = arrayfun(@(i) repmat(xyzS(i, :), [6, 1]), 1 : nSeat, 'Unif', false);
            xyz0 = vertcat(xyz0{:});
            patchCoords = patchCoords + xyz0;
            
            hg(end + 1) = patch(ht, 'Faces', vi, 'Vertices', patchCoords, ...
                'FaceColor', 'r', ...
                'EdgeColor', 'k', ...
                'Tag'      , 'Cabin Layout');
            
        end
        function hg = drawSeat(obj, ht)
            %drawSeat Draws a single seat object as a simple 'patch'.
            
            hg = [];
            
            patchCoords = obj.SeatPatchCoords;
            if isempty(patchCoords) %Escape route 
                return
            end
            
            if nargin < 2
                hF = figure;
                ht = axes('Parent', hF, 'NextPlot', 'add');
            end
            
            %Draw it!
            hg = patch(ht, ...
                'Faces', [1,2,5,6 ; 2,3,4,5], 'Vertices', patchCoords, ...
                'FaceColor', 'r', ...
                'EdgeColor', 'k');
            
        end
    end
    
    methods % defining the seat locations
        function updateSeatLayout(obj)
            %updateSeatLayout Generates a new seat layout based on the
            %cabin layout parameters.
            %
            %TODO - Update for the case where the Bluff-Body is not aligned
            %with the global X direction!
            
            %Parse - TODO: This will get superceded by Buildable methods
            if any(cellfun(@isempty, get(obj, {'SeatPitch', ...
                    'SeatWidth', 'AisleWidth', 'RowArrangement', 'NumSeats'})))
                return
            end
            
            L      = obj.CabinLength;
            rowKey = obj.RowArrangement;
            
            %Where are the seats and aisles?
            idxS = ismember(rowKey, 's');
            idxA = ismember(rowKey, 'a');
            
            %Calculate the number of rows
            nSeatPerRow = nnz(idxS);
            nRows       = ceil(obj.NumSeats ./ nSeatPerRow);
            
            %Update 'SeatLayout' key
            obj.SeatLayout = repmat({obj.RowArrangement}, [nRows, 1]);
            
            %Calculate coordinates of the centre of each seat
            %   - Position of each row along the cabin length
            rowEta = linspace(0, 1, nRows)';
            seatX   = rowEta .* L;
            %   - Width of each aisle/seat
            w = zeros(size(obj.RowArrangement));
            w(idxS) = obj.SeatWidth;
            w(idxA) = obj.AisleWidth;
            %   - Position of each seat starting at the forward-left edge
            %     of the cabin
            seatY(1)       = obj.EmptyWidth + w(1) ./ 2;
            seatY(2 : numel(w)) = seatY(1) + cumsum(w(1 : end - 1) ./ 2 + w(2 : end) ./ 2);
            %   - Only retain position of the seats!
            seatY = seatY(idxS);
            %   - Matrix of seat positions w.r.t forward left edge of cabin
            seatX = repmat(seatX, [1, nSeatPerRow]);
            seatY = repmat(seatY, [nRows, 1]);
            seatZ = zeros(size(seatX));
            SeatPosByRow.X = seatX;
            SeatPosByRow.Y = seatY;
            SeatPosByRow.Z = seatZ;
            
            %Assign to object
            obj.CabinSeatEta      = repmat(rowEta, [nSeatPerRow, 1]);
            obj.SeatPositions     = [seatX(:), seatY(:), seatZ(:)];
            obj.SeatPositionByRow = SeatPosByRow;
            
        end        
    end
    
end

