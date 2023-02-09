function [CAERO,INFO] = sortINTERPArray(CAERO,PARAM,INFO)

INFO.spline_type = unique(CAERO.Interp.Type);

ninterp = INFO.ninterp;
nset = INFO.nset;
ncaero = INFO.ncaero;

if (length(INFO.spline_type)>1)
    error('SPLINE types cannot be mixed up. Use only SPLINE1, SPLINE2 or SPLINE3 cards.');
end

if (INFO.amesh_av_vlm)
    
    dummy_lattice = CAERO.lattice_vlm;
    
    if ninterp > 0
        
        %fprintf(PARAM.FID, '\n\tSetting Aero interpolation database...');
        
        [index, i] = unique(CAERO.Interp.ID);
        
        if (length(i) ~= ninterp)
            
            n = (1 : nset);
            index = setdiff(n, i);
            error('Interpolation set entries have duplicated labels: %d', CAERO.Interp.ID(index(1)));
            
        end
        
        if ~nset
            
            error('No structural set given for interface procedure.');
            
        end
        % allocate data
        for n=1:ncaero
            % number of elements
            ne =  dummy_lattice.DOF(n, CAERO.geo.nelem(n), 2) - dummy_lattice.DOF(n, 1, 1) +1;
            CAERO.IS(n).data = zeros(ne, 1);
            
        end
        
        for n=1:ninterp
            
            i = find(CAERO.ID == CAERO.Interp.Patch(n));
            m = find(CAERO.Set.ID == CAERO.Interp.Set(n));
            
            if isempty(m)
                
                error('Unable to find structural set %d for the interpolation set %d.',...
                    CAERO.Interp.Set(n), CAERO.Interp.ID(n));
                
            else
                
                CAERO.Interp.Set(n)	= m; % substitute with set index for rapid accessing
                
                if isempty(i)
                    
                    error('Unable to find aero patch %d for the interpolation set %d.',...
                        CAERO.Interp.Patch(n), CAERO.Interp.ID(n));
                    
                else
                    
                    % check panel index
                    
                    if CAERO.Interp.Index(n, 2) > (1 + dummy_lattice.DOF(i, CAERO.geo.nelem(i), 2) - dummy_lattice.DOF(i, 1, 1))
                        
                        error('Panel index in interpolation set %d exceeds maximum panel indexing for patch %d.', ...
                            CAERO.Interp.ID(n), CAERO.ID(i));
                        
                    end
                    
                    CAERO.IS(i).data([CAERO.Interp.Index(n, 1) : CAERO.Interp.Index(n, 2)]) = n; % substitute with interp set index
                    
                end
                
            end
            
        end
        
        for n=1:ninterp
            i = find(CAERO.ID == CAERO.Interp.Patch(n));
            CAERO.Interp.Patch(n) = i; % substitute with patch index for rapid accessing
        end
        for n=1:ninterp
            if (isequal(CAERO.Interp.Type(n),1))
                if (CAERO.Interp.Param(n, 1)>0)
                    i1 = find(NODE.ID == CAERO.Interp.Param(n, 1));
                    if (isempty(i1))
                        error(['Node ', num2str(CAERO.Interp.Param(n, 1)),' in SPLINE1 ',num2str(CAERO.Interp.ID(n)),' does not exist.']);
                    end
                    i = find( CAERO.Set.Node(CAERO.Interp.Set(n)).data == i1);
                    if isempty(i)
                        error('Unable to find NODE %d in the interpolation set %d for SPLINE1 %d.',...
                            CAERO.Interp.Param(n, 1), CAERO.Set.ID(CAERO.Interp.Set(n)), CAERO.Interp.ID(n));
                    end
                    CAERO.Interp.Param(n, 1) = i;
                else
                    %fprintf(PARAM.FID, '\n\t\t### Warning: origin node not defined in SPLINE1 %d.\n\t\t    First node in SET1 %d will be used.', ...
                    %                    CAERO.Interp.ID(n), CAERO.Set.ID(CAERO.Interp.Set(n)));
                    CAERO.Interp.Param(n, 1) = 1;
                end
                if (CAERO.Interp.Param(n, 2)>0)
                    i2 = find(NODE.ID == CAERO.Interp.Param(n, 2));
                    if (isempty(i2))
                        error(['Node ', num2str(CAERO.Interp.Param(n, 2)),' in SPLINE1 ',num2str(CAERO.Interp.ID(n)),' does not exist.']);
                    end
                    i = find( CAERO.Set.Node(CAERO.Interp.Set(n)).data == i2);
                    if isempty(i)
                        error('Unable to find NODE %d in the interpolation set %d for SPLINE1 %d.',...
                            CAERO.Interp.Param(n, 1), CAERO.Set.ID(CAERO.Interp.Set(n)), CAERO.Interp.ID(n));
                    end
                    CAERO.Interp.Param(n, 2) = i;
                else
                    %fprintf(PARAM.FID, '\n\t\t### Warning: second node for y-axis not defined in SPLINE1 %d.\n\t\t    Last node in SET1 %d will be used.', ...
                    %                    CAERO.Interp.ID(n), CAERO.Set.ID(CAERO.Interp.Set(n)));
                    CAERO.Interp.Param(n, 2) = length(CAERO.Set.Node(CAERO.Interp.Set(n)).data);
                end
                if (CAERO.Interp.Param(n, 2) == CAERO.Interp.Param(n, 1))
                    error(['Node ', num2str(CAERO.Interp.Param(n, 1)),' is used twice to define SPLINE1 ',num2str(CAERO.Interp.ID(n)),' reference frame.']);
                end
            end
        end
        
        % check if any panel is not interfaced
        
        for n=1:ncaero
            
            ne = length(CAERO.IS(n).data);
            index = find(CAERO.IS(n).data); % get interfaced panels
            index = setdiff((1:ne), index); % get not interfaced panels
            il = length(index);
            %                dof = [dummy_lattice.DOF(CAERO.Interp.Patch(n), 1, 1) : dummy_lattice.DOF(CAERO.Interp.Patch(n), CAERO.geo.nelem(n), 2)];
            if il ~= 0
                fprintf(PARAM.FID, '\n\t\t### Warning: patch %d has %d unsplined panels.', CAERO.ID(n), il);
            end
            
            
        end
        
        %fprintf(PARAM.FID,'done.');
        
        
    end
    
end
end