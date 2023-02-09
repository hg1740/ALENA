
function [lattice, ref] = vlm_setup(fid, geo, state, ref)

%fprintf(fid, '\n\tSetting aerodynamic mesh for VLM\n');

[lattice, ref] = geosetup15(fid, geo, ref);

dim2 = size(lattice.VORTEX,2);

% check if vortex has wake
if dim2==8						
    %discard far wake points
    lattice.VORTEX = lattice.VORTEX(:, 2:7, :); 
    
end

 %appending wake lattice points (farpoints)
if state.AS~=0  
    
    %setting up wake legs.
    lattice = wakesetup2(lattice, state, ref); 
else
    
    terror(13)
    
end

[n, m] = find(geo.flapped');

% Do any flaps have rudder deflections?
if ~isempty(m)        
    
    noof_flaps=sum(sum(geo.flapped));
    
    %Loop all flaps and set them according to setting vector
    for k=1:noof_flaps 
        
        flap_no=k;
        if k == 28
            k;
        end
        deflection=(geo.flap_vector(m(k),n(k)));
        %[lattice]=setrudder3(flap_no,deflection,lattice,geo);
        
    end
    
end
end