classdef (ConstructOnLoad) Spar < awi.model.Component
    %Spar Defines a single spar that belongs to a LiftingSurface.
    
    %Position of spar along the wing
    properties (SetObservable, AbortSet)
        %Normalised position of the spar along the chord of the parent
        %LiftingSurface object.
        XLoc
        %Normalised position of the spar along the spanvector of the parent
        %LiftingSurface object.
        Eta
    end
    
    %Material/thickness properties of the spar
    properties
        %Normalised distribution of the spar thickness along the length of
        %the spar
        EtaThickness
        %Thickness of the spar
        Thickness
        %Handle to the 'awi.model.Material' object that contains the
        %material properties
        Material
    end
            
    methods % set / get
        function set.XLoc(obj, val)         %set.XLoc   
            %set.XLoc Set method for the property 'XLoc'.
            %
            % Rules
            %   - Type : numeric
            %   - Attr : row, nonnegative, <= 1
            
            val = awi.mixin.BeamProperty.validateDistribution(obj, val, 'XLoc');
            obj.XLoc = val;
        end
        function set.Eta(obj, val)          %set.Eta    
            %set.Eta Set method for the property 'Eta'.
            %
            % Rules
            %   - Type : numeric
            %   - Attr : row, nonnegative, <= 1
            
            val = awi.mixin.BeamProperty.validateDistribution(obj, val, 'Eta');
            obj.Eta = val;
        end 
        function set.EtaThickness(obj, val) %set.EtaThickness
            %set.EtaThickness Set method for the property 'EtaThickness'.
            %
            % Rules
            %   - Type : numeric
            %   - Attr : row, nonnegative, <= 1
            
            val = awi.mixin.BeamProperty.validateDistribution(obj, val, 'EtaThickness');
            obj.EtaThickness = val;
        end
        function set.Thickness(obj, val)    %set.Thickness
            %set.Thickness Set method for the property 'Thickness'.
            %
            %Rules
            %   - Type : numeric
            %   - Attr : row, nonnegative, 
            if iscolumn(val)
                val = val';
            end
            validateattributes(val, {'numeric'}, {'row', 'nonnegative', ...
                'finite', 'real', 'nonnan'}, class(obj), 'Thickness');
            obj.Thickness = val;
        end
    end
    
    methods % construction / destruction
        function obj = Spar(varargin)
            
            %Pass it on
            obj@awi.model.Component(varargin{:});
            
            %Configure collectables
            obj.IsLeafNode = true;
            
            %Extend property groups
            obj.addPropertyGroup('Geometry', ...
                'XLoc', ['Normalised position of the spar along the ', ...
                'chord of the parent LiftingSurface object.'], ...
                'Eta' , ['Normalised position of the spar along the ', ...
                'span of the parent LiftingSurface object.']);
                        
        end
    end
    
    methods % visualisation
        function hg = drawElement(obj, ht)
            %drawElement Class specific draw method. 
            %
            % For the 'Spar' class this method will draw a surface defining
            % the plane of the spar midline.
            
            hg = [];
            
            %Get coordinates of the spar in the parent frame
            [x, y, z] = calculateSparCoords(obj);
            if isempty(x) || isempty(y) || isempty(z)
                return
            end
            
            %Specify the (x,y,z) coordinates as a single continuous line 
            x = [x{1}, fliplr(x{2})];
            y = [y{1}, fliplr(y{2})];
            z = [z{1}, fliplr(z{2})];
            
            %Control nodes - Get rid
            plot3(ht, x, y, z, 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineStyle', 'none');
            
            %Define faces data by indexing vertices
            [nSpar, nEta] = size(x);
            ub    = cumsum(repmat(nEta, [1, nSpar]));
            lb    = [1, ub(1 : end - 1) + 1];
            faces = arrayfun(@(i) lb(i) : ub(i), 1 : nSpar, 'Unif', false);
            faces = vertcat(faces{:});
            
            %Store all vertex data            
            x      = num2cell(x', 1);
            y      = num2cell(y', 1);
            z      = num2cell(z', 1);
            vertex = [vertcat(x{:}), vertcat(y{:}), vertcat(z{:})];
            
            %Draw as a single patch object
            hg = patch(ht, 'Faces', faces, 'Vertices', vertex, 'FaceColor', 'r');
            
        end
    end
    
    methods % calculating the spar geometry
        function [x, y, z] = calculateSparCoords(obj)
            %calculateSparCoords Calculates the vertex positions of the
            %Spar object. 
            %
            % Detailed desciption:
            %   - Returns a [1, 2] cell for x, y and z coordinates. 
            %       + First element contains the coordinates of the upper 
            %         vertices of the spar.
            %       + Second element contains the coordinates of the lower
            %         vertices of the spar.
            
            %Start with nothing
            x = [];
            y = [];            
            z = [];

            %Get the parent
            ph = getSparParentLiftingSurface(obj);
            if isempty(ph)
                return
            end
            
            %Check parent LiftingSurface object has profile information
            Profiles = ph.CrossSection;
            if isempty(Profiles)
                return
            end
                        
            %Where are the profiles defined?
            %   - 'eta' here refers to 'r' domain -> convert to span vector
            profileEta   = [Profiles.BeamEta];
            [xd, yd, zd] = ph.xyzdata;
            r = ph.getLineLength(xd, yd, zd);
            coords = interp1(r / r(end), [xd ; yd ; zd]', profileEta);
            switch ph.SpanVector
                case 'Y'
                    profileEta = (coords(:, 2) / coords(end, 2))';
                case 'Z'                    
                    profileEta = (coords(:, 3) / coords(end, 3))';
            end
            
            %Set all spars to be defined at the same spanwise points -
            %Makes plotting marginally easier and allows vectorised ops
            [xLoc, ~, nEta, uEta] = interpolateSparXLoc(obj, profileEta);     
            
            %Spars populated?
            if isempty(xLoc)
                return
            end
            
            %Interpolate the cross-section data at the spar locations
            [xU, xL, zU, zL, ~, xzData] = interpolateNormXZData(Profiles, uEta);
            
            %Interpolate the cross-section data at each spanwise slice to
            %find the spar height at the required eta and x points
            nSp    = numel(obj);
            zSparU = zeros(nEta, nSp);
            zSparL = zeros(nEta, nSp);
            for i = 1 : nEta
                zSparU(i, :) = interp1(xU(i, :), zU(i, :), xLoc(:, i));
                zSparL(i, :) = interp1(xL(i, :), zL(i, :), xLoc(:, i));
            end
                         
            %Scale the coordinates            
            [xLoc, zSparU] = scaleProfileData(xLoc, zSparU', xzData, 1);
            [~   , zSparL] = scaleProfileData(xLoc, zSparL', xzData, 1);            

            %Get local orientation of cross-section at desired eta position
            rotMat  = arrayfun(@(P) P.Orientation.RMatrix(:), Profiles, 'Unif', false);
            rotMat  = horzcat(rotMat{:});
            rotMat_ = interp1([Profiles.BeamEta], rotMat', uEta, 'linear', 'extrap');
            
            %Transform local spar coords into global frame
            sparCoords = zeros(3, 2 * nSp, nEta);
            for iE = 1 : nEta
                r       = reshape(rotMat_(iE, :), [3, 3]);
                uCoords = r * [zeros(nSp, 1) , xLoc(:, iE) , zSparU(:, iE)]';
                lCoords = r * [zeros(nSp, 1) , xLoc(:, iE) , zSparL(:, iE)]';
                sparCoords(:, :, iE) = [uCoords, lCoords];
            end            
            
            %Shift the coordinates along the span
            %   - Assume 'eta' refers to the global 'SpanVector' axis.
            [xd_, yd_, zd_] = xyzdata(ph);
            if isempty(xd_) || isempty(yd_) || isempty(zd_)
                return
            end
            switch ph.SpanVector                
                case 'Y'
                    eta_  = yd_ ./ yd_(end);                    
                case 'Z'
                    eta_ = zd_ ./ zd_(end);
            end
            rData = interp1(eta_, [xd_ ; yd_ ; zd_]', uEta)';
            
            %Translate into local wing frame by adding span offsets
            x = squeeze(sparCoords(1, :, :)) + repmat(rData(1, :), [nSp * 2, 1]);
            y = squeeze(sparCoords(2, :, :)) + repmat(rData(2, :), [nSp * 2, 1]);
            z = squeeze(sparCoords(3, :, :)) + repmat(rData(3, :), [nSp * 2, 1]);
            
            %Return in cell format
            x = {x(1 : nSp, :), x(nSp + 1 : end, :)};
            y = {y(1 : nSp, :), y(nSp + 1 : end, :)};
            z = {z(1 : nSp, :), z(nSp + 1 : end, :)};
                        
        end
        function [xLoc, uEta, nEta, rEta] = interpolateSparXLoc(obj, profileEta, bProfileOnly)
            %interpolateSparXLoc Returns the spar x-offset at the unique
            %distribution of eta points for a collection of Spar objects.   
            %
            % Detailed Description:
            %   - The 'uEta' value specified here is along the parent
            %     LiftingSurface 'SpanVector'.
            %   - The 'rEta' value is the 'uEta' points specified along the
            %     straight-line axis of the beam.
            
            rEta = [];
            
            if nargin < 2
                profileEta = [];
            else
                validateattributes(profileEta, {'numeric'}, {'nonnegative', ...
                    '<=', 1, 'nonnan', 'finite', 'real'}, class(obj), 'profileEta');
            end
            if nargin < 3
                bProfileOnly = false;
            end
            
            eta    = get(obj, {'Eta'});
            xLoc   = get(obj, {'XLoc'});
            bEmpty = any(or(cellfun(@isempty, eta), cellfun(@isempty, xLoc)));
            bNumel = any(diff([cellfun(@numel, eta) , cellfun(@numel, xLoc)], [], 2) ~= 0);
            
            if bEmpty || bNumel
                xLoc = [];
                uEta = [];
                nEta = 0;
                return
            end
            
            if bProfileOnly
                uEta = profileEta;
            else
                uEta = unique([horzcat(eta{:}), profileEta]);
            end
            
            nEta = numel(uEta);
            ind = find(cellfun(@numel, xLoc) ~= nEta);
            for ii = 1 : numel(ind)
                xLoc{ind(ii)} = interp1(eta{ind(ii)}, xLoc{ind(ii)}, uEta, 'linear', 'extrap');
            end
            xLoc = vertcat(xLoc{:});
            
            ph = getSparParentLiftingSurface(obj);
            if isempty(ph)
                return
            end
            span = ph.Span;
            if isempty(ph)
                return
            end
           
            %Convert 'eta' positions to 'r' domain.
            [xd, yd, zd] = xyzdata(ph);
            switch ph.SpanVector
                case 'Y'                    
                    yRib = unique(uEta * span);
                    xRib = interp1(yd, xd, yRib);
                    zRib = interp1(yd, zd, yRib);
                case 'Z'
                    zRib = unique(uEta * span);
                    xRib = interp1(zd, xd, zRib);
                    yRib = interp1(zd, yd, zRib);                    
            end
            r      = ph.getLineLength(xRib, yRib, zRib);
            rEta = r ./ r(end);
            
        end
        function ph = getSparParentLiftingSurface(obj)
            %getSparParentLiftingSurace Returns the handle to the unique
            %'LiftingSurface' object that parents the collection of Spar
            %objects.
            %
            % TODO - Get "unique" working.
            
            ph = [];
            
            par = get(obj, {'Parent'});
            if any(cellfun(@isempty, par))
                return
            end
            par = par{1};
            %ph = unique(horzcat(ph{:}));
            %if numel(ph) > 1
            %    error('Update code to loop through different parents');
            %end
            
            if isa(par, 'awi.model.LiftingSurface') || numel(par) == 1
                ph = par;
            end
            
        end
    end
    
end

