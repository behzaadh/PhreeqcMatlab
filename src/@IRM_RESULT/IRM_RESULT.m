classdef IRM_RESULT < uint32
    enumeration
        IRM_OK            (0)   
        IRM_OUTOFMEMORY   (-1) 
        IRM_BADVARTYPE    (-2)  % Failure, Invalid VAR type
        IRM_INVALIDARG    (-3)  % Failure, Invalid argument
        IRM_INVALIDROW    (-4)  % Failure, Invalid row
        IRM_INVALIDCOL    (-5)  % Failure, Invalid column
        IRM_BADINSTANCE   (-6)  % Failure, Invalid rm instance id
        IRM_FAIL          (-7)  % Failure, Unspecified
    end
end

