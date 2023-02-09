% This function updates aerobeam nodes coordinates
function AERO_POS = update_aerobeam_node_nl(ngrid, NODE, MASTER_POS, MASTER_DR, MASTER_R)

AERO_POS = [];

for n = 1:ngrid

	AERO_POS(n).data = [];
	
	if ~isempty(NODE.Aero) && ~isempty(NODE.Aero.Index(n).data)
	
		AERO_POS(n).data = zeros(size(NODE.Aero.Coord(n).data));
		
		for j=1:length(NODE.Aero.Index(n).data)
		
			AERO_POS(n).data(:, j) = MASTER_POS(n, 1:3)' + MASTER_DR(:,:,n) * MASTER_R(:,:,n) * (NODE.Aero.Coord(n).data(:,j));
		end

	end

end
