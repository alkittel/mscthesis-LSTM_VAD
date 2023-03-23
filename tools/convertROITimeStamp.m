function t = convertROITimeStamp(timeStamps)
%CONVERTTIMESTAMP Convert MM:SS:FFF timestamps (char) from WSA recording
% corpus labels to 1xN array of values in seconds (double)

t = zeros(height(timeStamps),1);

for ii = 1:size(t)
    t(ii) = (...
        str2double(timeStamps(ii,1:2))*60 + ...
        str2double(timeStamps(ii,4:5)) + ...
        str2double(timeStamps(ii,7:end))*1e-3);
end

end