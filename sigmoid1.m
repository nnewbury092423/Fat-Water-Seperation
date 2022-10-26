function x = sigmoid(x)
%SIGMOID   Apply sigmoid activation
%
%   Y = SIGMOID(X) computes the sigmoid activation for input data X. X can
%   be a formatted or an unformatted dlarray. If X is a formatted dlarray,
%   output Y is a formatted dlarray with the same dimension labels as X. If
%   X is an unformatted dlarray, output Y is an unformatted dlarray.
%
%   The sigmoid operation is defined by
%       Y = 1 / (1 + exp(-X));
%
%   Example:
%       % Create input data as a formatted dlarray and compute sigmoid
%       % activation
%       x = dlarray(randn(8,8,3,16), 'SSCB'); 
%       y = sigmoid(x);
% 
%   See also BATCHNORM, DLARRAY, LEAKYRELU, RELU

%   Copyright 2019 The MathWorks, Inc.

% Extract the input data
fd = x.FormattedData;
x.FormattedData = [];
[fd,xdata] = extractData(fd);

% The dlarray methods should not accept logical x 
if islogical(xdata)
    error(message('deep:dlarray:LogicalsNotSupported'));
end

% Call the internal API
xdata = internal_sigmoid(xdata);

% Format is guaranteed not to have changed
fd = insertData(fd, xdata);
x.FormattedData = fd;
end

