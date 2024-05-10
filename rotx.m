function R = rotx(t, deg)
    
    assert((isreal(t) & isscalar(t)) | isa(t, 'sym'), ...
        'RTB:rotx:badarg', 'theta must be a real scalar');
    
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
        1   0    0
        0   ct  -st
        0   st   ct
        ];
end