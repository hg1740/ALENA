function [hF, PlotLoads] = plotBeamLoads(BeamLoads, FEM, method)
%plotBeamLoads Plots the disrtibution of beam loads across a set of beam
%objects.

%Parse
if nargin < 3
   method = 'Average'; 
end
validatestring(method, {'EndA', 'EndB', 'BothEnds', 'Average'});

nRes = numel(BeamLoads);
assert(numel(FEM) == nRes, 'Expected the number of results to match the number of FEM objects');

if isstruct(BeamLoads)
    %Parse feldnames
    expNam  = {'Beams', 'BM1', 'BM2', 'TS1', 'TS2', 'AF', 'TTRQ'};
    fNames  = fieldnames(BeamLoads(1));
    assert(all(ismember(expNam, fNames)), sprintf(['When defining the ' , ...
        'beam loads results as a MATLAB structure the following fields', ...
        'are required:', repmat('\n\t- %s', [1, numel(expNam)])], expNam{:}));
end


%Make the figures
hF = arrayfun(@(fem) figure('Name', fem.Name), FEM);

fNames = fieldnames(BeamLoads);
fNames(ismember(fNames, 'Beams')) = {'RCoords'};
PlotLoads = repmat(cell2struct(cell(size(fNames)), fNames), size(BeamLoads));

for iR = 1 : nRes
    
    %Grab node coordinates
    BeamNodes = [BeamLoads(iR).Beams.Nodes];
    NodeA     = BeamNodes(1, :);
    NodeB     = BeamNodes(2, :);
    xA = [NodeA.X];
    xB = [NodeB.X];
    
    %Preallocate
    nBeams = numel(BeamLoads(iR).Beams);
    nData  = 2 * nBeams;
%     nSC    = size(BeamLoads(iR).BM1, 1);
    r      = zeros(3, nData);
%     load   = zeros(nSC, nData);
    
    %Assign coordinate data for end-A/end-B
    r(1, 1 : 2 : end) = xA(1, :);
    r(1, 2 : 2 : end) = xB(1, :);
    r(2, 1 : 2 : end) = xA(2, :);
    r(2, 2 : 2 : end) = xB(2, :);
    r(3, 1 : 2 : end) = xA(3, :);
    r(3, 2 : 2 : end) = xB(3, :);
    
    %Plot along the length of the beam
    dXYZ    = diff(r(1, :)).^2 + diff(r(2, :)).^2 + diff(r(3, :)).^2;
    rCoords = cumsum([0, sqrt(dXYZ)]);
    rCoords = rCoords ./ (rCoords(end));
    
    %Make the axes
    qNames = {'BM1', 'BM2', 'TS1', 'TS2', 'AF', 'TTRQ'};
    hAx   = arrayfun(@(i) subplot(3, 2, i, 'Parent', hF(iR), 'NextPlot', 'add'), 1 : numel(qNames));
    
    %Plot each load
    for iQ = 1 : numel(qNames)
        
        %Assign data   
        switch method            
            case 'EndA'
                rPlot = rCoords(1 : end - 1);
                load  = BeamLoads(iR).(qNames{iQ})(:, :, 1);
            case 'EndB'
                rPlot = rCoords(2 : end);
                load  = BeamLoads(iR).(qNames{iQ})(:, :, 2);
            case 'BothEnds'
                rPlot = rCoords;
                load(:, 1 : 2 : end) = BeamLoads(iR).(qNames{iQ})(:, :, 1);
                load(:, 2 : 2 : end) = BeamLoads(iR).(qNames{iQ})(:, :, 2);
            case 'Average'
                rPlot = [rCoords(1 : 2 : end), rCoords(end)];
                endAB = cat(3, ...
                    BeamLoads(iR).(qNames{iQ})(:, 2 : end, 1), ...
                    BeamLoads(iR).(qNames{iQ})(:, 1 : end - 1, 2));
                average = mean(endAB, 3);
                load  = [BeamLoads(iR).(qNames{iQ})(:, 1, 1), average, BeamLoads(iR).(qNames{iQ})(:, end, 2)];
        end

        %Stash
        PlotLoads(iR).(qNames{iQ}) = load;
        PlotLoads(iR).RCoords      = rPlot;
        
        %Plot
        plot(hAx(iQ), rPlot, load);
        title(hAx(iQ), qNames{iQ});
        
    end
    
end

end

