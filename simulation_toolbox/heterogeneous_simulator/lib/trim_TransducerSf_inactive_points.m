function T = trim_TransducerSf_inactive_points(TransducerSf)
%TRIM_TRANSDUCERSF_INACTIVE_POINTS Keep only transducer points with nonzero velocity amplitude.
%
%   Returns a copy of TransducerSf with xGrid, yGrid, zGrid and
%   complexVelocityAmplitude restricted to entries where
%   abs(complexVelocityAmplitude) > eps('single'). Other fields are copied
%   unchanged. The caller's struct is not modified.

mask = abs(TransducerSf.complexVelocityAmplitude) > eps('single');
if ~any(mask(:))
    error(['All transducer grid points have negligible complexVelocityAmplitude ', ...
        '(|v| <= eps(''single'')).']);
end
T = TransducerSf;
T.xGrid = TransducerSf.xGrid(mask);
T.yGrid = TransducerSf.yGrid(mask);
T.zGrid = TransducerSf.zGrid(mask);
T.complexVelocityAmplitude = TransducerSf.complexVelocityAmplitude(mask);
end
