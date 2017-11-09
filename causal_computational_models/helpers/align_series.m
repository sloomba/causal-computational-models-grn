function [new_series] = align_series(series, lag, len)
    count = idivide(length(series), int32(len));
    new_series = [];
    if lag >= 0
        for i=1:count
            new_series = cat(2, new_series, series(len*(i-1) + lag + 1 : len*i));
        end
    else
        for i=1:count
            new_series = cat(2, new_series, series(len*(i-1) + 1 : len*i + lag));
        end
    end
end