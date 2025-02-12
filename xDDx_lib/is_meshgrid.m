% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function isMeshgrid = is_meshgrid(X, Y)

if ~isequal(size(X), size(Y))
    isMeshgrid = false;
    return;
end

if ~ismatrix(X)
    isMeshgrid = false;
    return;
end

xVec = unique(squeeze(X(1,:)), 'stable');
yVec = unique(squeeze(Y(:,1)), 'stable');

[Xtest, Ytest] = meshgrid(xVec, yVec);

isMeshgrid = isequal(X, Xtest) && isequal(Y, Ytest);

end