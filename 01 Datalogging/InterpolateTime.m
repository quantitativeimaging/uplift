function value_plot = InterpolateTime(time_plot,time_measured,value_measured)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   File name:    InterpolateTime.m
%
%   Purpose  :    Interpolate values from raw data based on a specified
%                 time vector.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

value_plot = zeros(length(time_plot),1);

for row = 1:length(time_plot)
    % find the location of start of interpolation
    index = max(find(time_measured < time_plot(row)));

    if isempty(index)
        % take the first value if out of range of measured time
        value_plot(row) = value_measured(1);
    else
        % linear interpolate to next point
        gradient = (value_measured(index+1) - value_measured(index))/(time_measured(index+1)-time_measured(index));
        value_plot(row) = value_measured(index) + gradient *(time_plot(row)-time_measured(index));
    end
    
end

end