function delta = ObjFunction(xSam)

load('C:\Users\Rc15645\Google Drive\AEROGUST\Code\Beam Model\v2\Results Files\AeroGust\UQ\MCS\MCS_sample.mat')
delta = zeros(size(xSam));
alpha = [];

for ii = 1:length(xSam)
    ind = find(rho_MCS == xSam(ii));
    load(['C:\Users\rc15645\Google Drive\AEROGUST\Code\Beam Model\v2\Results Files\AeroGust\UQ\MCS\UAV Trim\trim_results_ff1_rr',num2str(ind),'.mat'])
    delta(ii) = alpha*180/pi;
end

% load('C:\Users\rc15645\Google Drive\AEROGUST\Model DB\UAV v2 wing\UQ\1\Nastran_ff1.mat')
% delta = NastranAoA(1:length(xSam));