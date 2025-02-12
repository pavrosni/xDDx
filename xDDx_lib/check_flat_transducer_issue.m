% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [doAlignmentUpdated, warningFlatTransducerIssue] = check_flat_transducer_issue(isSphericalSource, doAlignment)
doAlignmentUpdated = doAlignment;
warningFlatTransducerIssue = [];
if (~isSphericalSource) && doAlignment
    warningFlatTransducerIssue = 'Quick start auto alignment option is available for spherically shaped transducers only. The script was executed without auto alignment.';
    doAlignmentUpdated = false;
end


