classdef PointLoad < awi.fe.FEBaseClass
    %PointLoad Describes a load applied to a node in a finite element
    %model. The load can either be in a fixed coordinate system (e.g. the
    %global system) or it can follow the local frame as the structure
    %deflects (e.g. a follower load).
    %
    % The definition of the 'PointLoad' object matches that of the FORCE
    % bulk data type from MSC.Nastran.
    
    %Primary properties
    properties
        %Identification number of the load
        SID
        %Coordinate system that the load is applied in {'local', 'global'}
        LoadCoordSys = 'global';
        %Determines whether the load is a follower or non-follower load
        LoadBehaviour = 'non-follower';
        %Type of load {'force', 'moment'}
        LoadType = 'force';
        %Magnitude of the applied load
        Magnitude
        %Orientation of the applied load in the coordinate system defined
        %by 'LoadCoordSys'
        Orientation
        %Identification number of the node the load is associated with
        GID
    end
    
    %Store a reference to the 'awi.fe' objects
    properties (Hidden = true)
        %Handle to the 'awi.fe.Node' object that this load is associated
        %with
        Node
        %Handle to the 'awi.fe.Node' object that defines the orientation of
        %the node - only for 'LoadBehaviour = 'follower'
        OrientationNode
    end
    
    %Visualisation
    properties
        %Maximum vector length for plotting loads
        MaxVectorLength = 1;
    end
    
    methods % set / get
        function set.SID(obj, val)              %set.SID
            %set.GID Set method for the property 'SID'.
            %
            % Passes the value straight to the inherited 'ID' property
            % which will validate and store the ID number.
            obj.ID = val;
        end
        % allows point load to be assigned onto a random nodes e.g. RBE2
        % nodes
        function set.GID(obj, val)              %set.SID
            %set.GID Set method for the property 'SID'.
            %
            % Passes the value straight to the inherited 'ID' property
            % which will validate and store the ID number.
            obj.ID = val;
        end
        function set.LoadCoordSys(obj, val)     %set.LoadCoordSys 
            %set.LoadCoordSys Set method for the property 'LoadCoordSys'.            
            validatestring(val, {'local', 'global'}, class(obj), 'LoadCoordSys');
            obj.LoadCoordSys = lower(val);            
        end
        function set.LoadBehaviour(obj, val)    %set.LoadBehaviour
            %set.LoadBehaviour Set method for the property 'LoadBehaviour'.
            val = lower(val);
            validatestring(val, {'follower', 'non-follower'}, class(obj), 'LoadBehaviour');
            obj.LoadBehaviour = val;
        end
        function set.LoadType(obj, val)         %set.LoadType 
            %set.LoadType Set method for the property 'LoadType'.            
            validatestring(val, {'force', 'moment'}, class(obj), 'LoadType');
            obj.LoadType = lower(val);
        end
        function set.Magnitude(obj, val)        %set.Magnitude
           %set.Magnitude Set method for the property 'Magnitude'.           
           validateattributes(val, {'numeric'}, {'scalar', 'real', 'nonnan'}, ...
               class(obj), 'Magnitude');
           obj.Magnitude = val;
        end
        function set.Orientation(obj, val)      %set.Orientation
            %set.Orientation Set method for the property 'Orientation'.            
            validateattributes(val, {'numeric'}, {'column', 'numel', 3, ...
                'real', 'nonnan', '<=', 1}, class(obj), 'Orientation');
            obj.Orientation = val;
        end
        function set.Node(obj, val)             %set.Node
            %set.Node Set method for the property 'Node'.            
            validateattributes(val, {'awi.fe.Node'}, {'scalar'}, ...
                class(obj), 'Node');
            obj.Node = val;
        end
        function set.OrientationNode(obj, val)  %set.OrientationNode
            %set.OrientationNode Set method for the property 'OrientationNode'.            
            validateattributes(val, {'awi.fe.Node'}, {'scalar'}, ...
                class(obj), 'OrientationNode');
            obj.OrientationNode = val;
        end
        function set.MaxVectorLength(obj, val)  %set.MaxVectorLength
            %set.MaxVectorLength Set method for the property
            %'MaxVectorLength'.            
            validateattributes(val, {'numeric'}, {'scalar', 'real', ...
                'positive'}, class(obj), 'MaxVectorLength');
            obj.MaxVectorLength = val;
        end
        function val = get.SID(obj)             %get.SID
           %get.SID Get method for the property 'SID'.
           val = obj.ID;           
        end
        function val = get.GID(obj)             %get.GID
            %get.GID Get method for the property 'GID'.
            %
            % If the object has been assigned a handle to its 
            % 'awi.fe.Node' object then always use their ID number, else
            % use GID.
            
            if isempty(obj.Node)
                val = obj.PID;
            else
                val = obj.Node.ID;                
            end
            
        end
    end
    
    methods % constructor
        function obj = PointLoad
            
            %Make a note of the property names
            addFEProp(obj, 'SID', 'LoadCoordSys', 'Magnitude', 'Orientation');
            
        end
    end
    
    methods % visualisation
        function hg = drawElement(obj, ha)
            %drawElement Draws the load objects as a line object and
            %returns a single graphics handle for all the loads in the
            %collection. 
            %
            % Accepts a vector of objects.
            
            hg   = [];
            nObj = numel(obj);
            
            %Check the nodes are defined            
            Nodes = [obj.Node];            
            if nObj ~= numel(Nodes)
                warning(['Some of the ''awi.fe.PointLoad'' objects are ', ...
                    'not associated with a ''Node'' object. Update '    , ...
                    'objects before proceeding.'])
                return
            end
            
            %Check orientation is defined
            v = [obj.Orientation];
            if size(v, 2) ~= nObj
                warning(['Some of the ''awi.fe.PointLoad'' objects do ', ...
                    'not have an orientation vector. Update objects '  , ...
                    'before proceeding.'])
                return
            end
            
            %Scale loads
            mag = [obj.Magnitude];
            sf = max([obj.MaxVectorLength]) .* mag ./ max(abs(mag));
            v = v.* sf;
            
            %Rotate load orientation into local system (if applicable)
            idxL  = ismember({obj.LoadCoordSys}, 'local');
            if any(idxL)
                ind = find(idxL);
                %Multiply by rotation matrix for the local coordinate
                %system at the node
                R = getRotationMatrix([Nodes(idxL).OutputCoordSys]);
                v_ = arrayfun(@(i) R(:, :, i) * v(:, ind(i)), 1 : nnz(idxL), 'Unif', false);
                v(:, idxL) = horzcat(v_{:});                
            end
            
            %Generate coordinates for plotting 
            rN = [Nodes.X];
            rL = rN + v;
            xL = obj.padCoordsWithNaN([rN(1, :) ; rL(1, :)]);
            yL = obj.padCoordsWithNaN([rN(2, :) ; rL(2, :)]);
            zL = obj.padCoordsWithNaN([rN(3, :) ; rL(3, :)]);   
            
            %Plot
            hg(1) = line(ha, xL, yL, zL, ...
                'Color'    , 'm', ...
                'LineWidth', 2  , ...
                'Tag', 'Coordinate Systems');
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.Load' object into a
            %text file using the format of the MSC.Nastran 'FORCE',
            %'FORCE1', 'FORCE2', 'MOMENT', 'MOMENT1' or 'MOMENT2' bulk data
            %entries.         
            
            %By default, do not close the file
            bClose = false;
                        
            if nargin < 2 %Ask the user for the file
                fName  = awi.fe.FEBaseClass.getBulkDataFile;
                bClose = true;
                fid    = fopen(fName, 'w');                
            end
            
            if nargin < 3 %Comments by standard
                bComment = true;
            end
            
            %Split into follower loads and global loads
            idxFollower = ismember({obj.LoadBehaviour}, 'follower');
            
            %Follower load (forces & moments)
            FLoad     = obj(idxFollower);
            FLoadType = {FLoad.LoadType};            
            i_writeFollowerLoad(fid, bComment, FLoad(ismember(FLoadType, 'force')));
            i_writeFollowerLoad(fid, bComment, FLoad(ismember(FLoadType, 'moment')));
            
            %Global load (forces & moments)
            GLoad = obj(~idxFollower);
            GLoadType = {GLoad.LoadType}; 
            i_writeGlobalLoad(fid, bComment, GLoad(ismember(GLoadType, 'force')));
            i_writeGlobalLoad(fid, bComment, GLoad(ismember(GLoadType, 'moment')));
            
            if bClose %Close the file?
                fclose(fid);
            end
            
            function i_writeFollowerLoad(fid, bComment, Load)
                %i_writeFollowerLoad Writes a 'FORCE1' or 'MOMENT1' card
                %using the data in 'Load'.
                
                if isempty(Load) %Escape route
                    return
                end
                
                if bComment %Helpful comments?     
                    switch Load(1).LoadType                        
                        case 'force'
                            comment = ['FORCE1 : Defines a concentrated ', ...
                                'force at a grid point by specification ', ...
                                'of a magnitude and two grid points that ', ...
                                'determine the direction.']; 
                            namestr = {'FORCE1'};
                        case 'moment'
                            comment = ['MOMENT1 : Defines a concentrated ', ...
                                'moment at a grid point by specifying a ' , ...
                                'magnitude and two grid points that '     , ...
                                'determine the direction.'];
                            namestr = {'MOMENT1'};
                    end
                    awi.fe.FEBaseClass.writeComment(comment, fid);
                end
                                
                awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
                                               
                %Card name
                nam  = repmat(namestr, [1, numel(Load)]);
                
                %Grab data
                GID_       = {Load.GID};
                NodeOrient = [Load.OrientationNode];
                
                %Set up the format for printing
                data = [nam ; {Load.ID} ; GID_ ; {Load.Magnitude} ; GID_ ; {NodeOrient.ID}];
                
                %Write in 16-character column width as standard
                format = '%-8s%-8i%-8i%#-8.2E%-8i%-8i\r\n';
                
                %Write the data to the file
                fprintf(fid, format, data{:});
                
            end
            
            function i_writeGlobalLoad(fid, bComment, Load)
                %i_writeGlobalLoad Writes a 'FORCE' or 'MOMENT' card using
                %the data in 'Load'.
                
                if isempty(Load)
                    return
                end
                
                if bComment %Helpful comments?     
                    switch Load(1).LoadType                        
                        case 'force'
                            comment = ['FORCE : Defines a static '   , ...
                                'concentrated force at a grid point ', ...
                                'by specifying a vector.']; 
                            namestr = {'FORCE'};
                        case 'moment'
                            comment = ['MOMENT : Defines a static '   , ...
                                'concentrated moment at a grid point ', ...
                                'by specifying a vector.'];
                            namestr = {'MOMENT'};
                    end
                    awi.fe.FEBaseClass.writeComment(comment, fid);
                end
                                
                awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
                                               
                nObj = numel(Load);
                
                %Card name
                nam  = repmat(namestr, [1, nObj]);
                
                %Grab data
                CID        = zeros(1, nObj);
                idxL       = ismember({Load.LoadCoordSys}, 'local');
                LoadNode   = [Load.Node];
                CID(idxL)  = [LoadNode(idxL).CD];
                G          = num2cell([LoadNode.ID]);
                CID        = num2cell(CID);
                F          = num2cell([Load.Magnitude]);
                X          = num2cell([Load.Orientation]);   
                
                %Set up the format for printing
                data = [nam ; {Load.ID} ; G ; CID ; F ; X(1, :) ; X(2, :) ; X(3, :)];
                
                %Write in 16-character column width as standard
                format = '%-8s%-8i%-8i%-8i%#-8.1f%#-8.3f%#-8.3f%#-8.3f\r\n';
                
                %Write the data to the file
                fprintf(fid, format, data{:});

            end
            
        end
    end
    
end

