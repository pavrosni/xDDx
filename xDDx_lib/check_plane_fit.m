% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function isPlane = check_plane_fit(zGrid, tolerance)

isPlane = (max(zGrid(:))-min(zGrid(:)) < tolerance);

end

