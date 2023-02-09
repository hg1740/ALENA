function poly = legendre1D(samples, order)

%Outputs 1D legendre polynomials at a number of sample points, for an input
%of any dimension. Inputs are, samples: nSamples x nDimensions matrix, and
%order: maximum order of polynomial to calculate. 

%Output poly, is a nDimensions x 1 cell array, with each entry corresponding
% to a dimension of the input, in the form of a nSamples x order + 1 matrix
%containing each of the legendre polynomials up to the required order,
%evaluated at each sample point, for the dimension in question.

D = size(samples,2);
n = size(samples,1);
poly = cell(D,1);

% Loop through each dimension
for i = 1:D
    currentPoly = zeros(n,order+1);
    % Define the first two polynomials (x0 = 1, x1 = x)
    currentPoly(:,1) = ones(n,1);
    currentPoly(:,2) = samples(:,i);
    for j = 2:order
        % Recursive definition of the Legendre Polynomials
        currentPoly(:,j+1) = (samples(:,i).*currentPoly(:,j)*(2*j-1)-(j-1)*currentPoly(:,j-1))./j;
    end
    poly{i,1} = currentPoly;
end