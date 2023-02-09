function basis = constructBasis2(polys)

% Assemble a d dimensional polynomial basis evaluated a N samples points, 
% through taking the tensor product of d 1D polynomials trunctated up to 
% maximum order p. 

% polys is a numDimensions x 1 cell array, with each entry corresponding to
% a particular dimenision. Each entry contains a nSamples x order + 1
% matrix of the 1D polynomials evaluated at each sample point. The output
% basis is a nSamples x (p + d)!/p!d! matrix of the basis polynomials
% evaluated at each sample point

numDimensions = size(polys,1);
order = size(polys{1,1},2)-1;
N = size(polys{1,1},1);

if numDimensions == 1
    % If there is only one dimension there's no need to run the rest of the
    % code
    basis = polys{1,1};
else  
    % Set it up such that the 'working matrix' for the first iteration is a 
    % vector of the 1st dimension 1D polynomials, to be multiplied by the 
    % vector of 1D polynomials for the 2nd dimension.
    iDimensionProduct = polys{1,1}; % Working matrix of polynomials
    orderMatrix = 0:order; % Counter which keeps track of the polynomial order
    % Loop through remaining dimensions
    for i = 1 : numDimensions - 1
        counter = 1; % Counter which tracks the size of the matrix
        % At each iteration the reset the 'working matrix' to be the
        % polynomials outputed at the previous dimension.
        ithPolynomial = polys{i+1,1}; % Get the next dimension polynomials
        numCoefficients = size(iDimensionProduct,2);
        workingMatrix = iDimensionProduct;
        oldOrderMatrix = orderMatrix; % Keep track of polynomial orders at previous iteration
        % Re-initialise polynomial and order matrices for the new dimension
        iDimensionProduct = zeros(N,factorial(order+i+1)/(factorial(order)*factorial(i+1)));
        orderMatrix = zeros(1,factorial(order+i+1)/(factorial(order)*factorial(i+1))); 
        %Loop over every element in the current dimension 1D polynomials
        for j = 1 : order + 1
            %Loop over every element in the previous dimensions'
            %polynomials
            for k = 1 : numCoefficients
                %Checks the product of two terms is less than or equal to
                %the overall desired order p
                if j - 1 + oldOrderMatrix(1,k) <= order
                    % New polynomials given by product of the two
                    % dimensions
                    iDimensionProduct(:,counter) = ithPolynomial(:,j) .* workingMatrix(:,k);
                    %Keeps track of the order of the new polynomials
                    orderMatrix(1,counter) = oldOrderMatrix(1,k) + j - 1;
                    counter = counter + 1;
                end
            end
        end
    end
    basis = iDimensionProduct;
end