function R = roty(t, deg)
    
    assert((isreal(t) & isscalar(t)) | isa(t, 'sym'), ...
        'RTB:roty:badarg', 'theta must be a real scalar');
    
    if nargin > 1 && strcmp(deg, 'deg')
        t = t *pi/180;
    end
    
    ct = cos(t);
    st = sin(t);
    
    % make almost zero elements exactly zero
    if ~isa(t, 'sym')
        if abs(st) < eps
            st = 0;
        end
        if abs(ct) < eps
            ct = 0;
        end
    end
    
    % create the rotation matrix
    R = [
        ct  0   st
        0   1   0
       -st  0   ct
       ];
end