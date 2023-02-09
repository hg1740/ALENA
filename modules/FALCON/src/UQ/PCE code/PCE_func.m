%% PCE
% Expansion Input Parameters
function [y,fy,data] = PCE_func(numVar,order,numSam,mu,stdev,x_in,y_in,dist_flag,data_flag,quantile_vec,xEm)

if ~exist('dist_flag','var')
    dist_flag = 1;
end
if ~exist('data_flag','var')
    data_flag = 0;
end

% numVar = 1; % Number of random variables (Expansion Dimension d)
% order = 4; % Order of Polynomials used (p) - higher order plynomials should give better results
numEm = 10000; % Number of emulations in surrogate model (arbitrarily large number)

% The minimum number of samples required is given by the number of
% expansion coefficients as given by (order + dimension)!/(order!dimension!)
% Typically two or three times this number are required for good results

% numSam = 50; % Number of Samples used for the regression model

distribution = 'normal';
% mu  = [0.151867291     ]; % variation of E and v %CV=0.1
% stdev = [0.151867291/10/3];% 1.03e9 5.9e8 3.66e-6];% 0.03 1.4e9 4e-9] % standard deviations of variables
param = [mu;stdev]; % input required for sampling
% Observations are for input PDF plots which are included for illustrative
% purposes. Delete if not required.
observations = normrnd(repmat(mu,5000,1),repmat(stdev,5000,1));

%--------------------------------------------------------------------------

% EMPIRICAL - may be used for any arbitrary (yet independent) distribution

% Requires param, a numObservations x numVariables matrix containing a
% number of observations of each of the inputs. See below for example use.

% distribution = 'empirical'
% lower = [4800 1.55 67e9 1.9e-7];
% upper = [5200 1.45 73e9 2.1e-7]; 
% observations = betarnd([ones(5000,1)*5 ones(5000,1)*1 ones(5000,1)*2 ones(5000,1)*2],[ones(5000,1)*1 ones(5000,1)*3 ones(5000,1)*5 ones(5000,1)*2]).*repmat(upper - lower,5000,1) + repmat(lower,5000,1);
% % Beta distributions used for illustrative purposes, could be any
% % distribution, lower and upper bounds are not required either
% param = observations;

% UNIFORM (best used when only upper and lower bounds are known)

% Requires input of param, a 2 x numVariables matrix in the form 
% [lowerBound;upperBound] for each of the variables. See below for example use.
% 
% distribution = 'uniform'
% lower = [4800 1.55 67e9 1.9e-7] %lower bounds ([P L E I] in this example)
% upper = [5200 1.45 73e9 2.1e-7] %upper bounds
% param = [lower;upper]; % input required for sampling
% % Observations are for the input PDF plots which are included for
% % illustrative purposes. Delete if not required.
% observations = rand(5000,4).*repmat(upper - lower,5000,1) + repmat(lower,5000,1);

%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate CDF of sample points
% FxSam = lhsdesign(numSam,numVar,'criterion','maximin');
% load('C:\Users\Rc15645\Google Drive\AEROGUST\Code\Beam Model\v2\Results Files\AeroGust\UQ\MCS\MCS_sample.mat')

switch distribution
    case 'normal'
        mu    = repmat(param(1,:),numSam,1);
        stdev = repmat(param(2,:),numSam,1);
        % Take inverse of normal CDF based on input mean and standard
        % deviation values
%         xSam = norminv(FxSam,mu,stdev);
        xSam = x_in;%rho_MCS(1:numSam);
        % Transform sample points onto standard normal variables
        xTrans = (xSam - mu)./stdev;
        % Calculate Hermite Polynomials at Sample Points
        poly = hermite1D(xTrans, order);
        % Generate a number of standard, normally distributed values as
        % input to the surrogate model
        if exist('xEm','var') && ~isempty(xEm)
        else
            xEm = normrnd(zeros(numEm,numVar),ones(numEm,numVar));
        end
        % Calculate the Hermite Polynomials at each emulation input
        polyEm = hermite1D(xEm, order);
    case 'uniform'
        % The CDF for the uniform distribution on [0,1] is Fx = x,
        % therefore the sample points are simply given by Fx
        poly = legendre1D(FxSam, order); % Calculate Legendre Polynomials at sample points
        % Generate a number of uniformly distributed values on [0,1] as
        % input to the surrogate model
        xEm = rand(numEm,numVar);
        % Calculate the Legendre Polynomials at each emulation input
        polyEm = legendre1D(xEm, order);
        % Transform sample points from the interval [0,1] to that defined
        % in the inputs
        xSam = FxSam.*repmat(param(2,:) - param(1,:),numSam,1) + repmat(param(1,:),numSam,1);
    case 'empirical'
        % Transforms the inputs onto a uniform distribution and uses the
        % Legendre polynomials as a basis. (Use of Gram-Schmidt Polynomials 
        % and un-transformed variables may be more efficient).
        xSam = zeros(numSam,numVar);
        % Use the Kernal Density Estimate to calculate the input values
        % which correspond to the required sample CDF
        for i = 1:numVar
            xSam(:,i) = ksdensity(param(:,i),FxSam(:,i),'function','icdf');
        end       
        % Calculated the legendre polynomial each each sample point in
        % transformed, uniform space on [0,1] (note this is simply the sample 
        % CDF since Fx = x for the transformed variable)
        poly = legendre1D(FxSam, order);
        % Generate a number of uniformly distributed values on [0,l] for
        % use in the surrogate model
        xEm = rand(numEm,numVar);
        % Calculate the Legendre Polynomials at each of the emulation input
        % values
        polyEm = legendre1D(xEm, order);
end


basis = constructBasis2(poly);
% Fitted Regression model

% This fits a surrogate model in the form: y = basis * beta + error
% The model is fitted through minimising the error in a least squares
% sense. At this point, the output value of your model at each of the 
%sample points in 'xSam' are required.

%REPLACE THIS WITH YOUR OWN FUNCTION.
% a numSamples x 1 vector of outputs is required
% ySam = ObjFunction(xSam); % Output values for the regression model. 
ySam = y_in;

%Fit regression model through least squares
beta = (basis'*basis)\(basis'*ySam);

% Calculate emulated y values using surrogate model and emulation input
% values
basisEm = constructBasis2(polyEm);
yEm = basisEm*beta;

% Perform Monte Carlo Simulation (MCS) of Model
% 
% % REPLACE THIS WITH YOUR OWN FUNCTION IF REQUIRED
% observations = normrnd(repmat(mu,200,1),repmat(stdev,200,1));
% ss=size(observations)
% yMCS = ObjFunction(observations);

% [fyMCS ybins] = hist(yMCS,25); % Calculate MCS PDF using histogram as above
% %fyMCS = fyMCS/(size(observations,1)*(ybins(1,2)-ybins(1,1)));

if dist_flag
    for ii = 1:size(yEm,2)
        [fy(ii,:),y(ii,:)] = ksdensity(yEm(:,ii));%,min(yEm):(max(yEm)-min(yEm))/(1000-1):max(yEm) Calcuates PCE PDF using kernal density estimate (hist may also be used)
    end
else
    fy = [];
    y  = [];
end
%[fy y] = ksdensity(ySam);

if data_flag
    data.mean = mean(yEm);
    data.std  = std(yEm);
    data.skew = skewness(yEm);
    data.kurt = kurtosis(yEm);
    if exist('quantile_vec','var') && ~isempty(quantile_vec)
        data.quantile = quantile(yEm,quantile_vec);
    end
else
    if exist('quantile_vec','var') && ~isempty(quantile_vec)
        data.quantile = quantile(yEm,quantile_vec);
    else
        data = [];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PDF response plot

% Plot of PCE outputs and comparison with MCS for validation

% figure('Name','Output PDF','NumberTitle','off')
% % bar(ybins,fyMCS,'g')
% hold on
% plot(y,fy,'-r','LineWidth',2)%-trapz(y,y.*fy)
% %line([330 330], [0 0.02],'LineStyle','--');
% %legend ('Deterministic Flutter','FontSize', 12)
% %legend('boxoff')
% %hold off
% title('Polynomial Chaos (PCE)','FontSize', 12)
% xlabel('Deflection (m)','FontSize', 12)
% ylabel('PDF','FontSize', 12)