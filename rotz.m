function R = rotz(t, deg)
    
    assert((isreal(t) & isscalar(t)) | isa(t, 'sym'), ...
        'RTB:rotz:badarg', 'theta must be a real scalar or symbolic');

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
        ct  -st  0
        st   ct  0
        0    0   1
        ];
end