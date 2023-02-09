%% STRUCT2MAT   Convert field value of a strcuture to a vector
%
%  STRUCT2MAT(STRUCTURE,FIELDNAME,INDEX) extracts the vector of values from
%  a field in a structure.
%
%  Inputs:
%      STRUCTURE : A structure
%      FIELDNAME : The field to be extracted - character string
%      INDEX     : The column index if the field contains a matrix
%
%
%   Example:
%      a.b = eye(4);
%      dat = struct2mat(a,'b',2);
%
%      will extract the second column of the field 'b'

%   Copyright 2016 University of Bristol
%   Private function.
function y = struct2mat( Value, fieldName, i )

x = cell2mat({Value.(fieldName)}');

if exist('i','var') == 0
    y = x;
else
    y = x(:,i);
end

end

