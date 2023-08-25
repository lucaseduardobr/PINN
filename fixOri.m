function [x_corrected, y_corrected] = fixOri(x, y)
    % check x e y dimensions
    if iscolumn(x)
        x_corrected = x;
    else
        x_corrected = x';
    end

    if iscolumn(y)
        y_corrected = y;
    else
        y_corrected = y';
    end
end