% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function isSphere = check_sphere_fit(xGrid, yGrid, zGrid, tolerance, radiusOfCurvature)

diffMax = max(abs(xGrid(:).^2 + yGrid(:).^2 + (zGrid(:) - radiusOfCurvature).^2 - radiusOfCurvature^2));
isSphere = (diffMax < tolerance);

end
