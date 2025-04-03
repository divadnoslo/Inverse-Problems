function addValuesToBarChart(b)
% Add values to bar chart
%
%

%% Implementation

xtips = b.XEndPoints;
ytips = b.YEndPoints;
values = b.YData;
labels = strings(1, length(values));
for k = 1 : length(values)
    labels(k) = sprintf("%6.3e", values(k));
    if (values(k) >= 0)
        text(...
            xtips(k), ytips(k), labels(k), ...
            "HorizontalAlignment", "center", ...
            "VerticalAlignment", "bottom")
    else
        text(...
            xtips(k), ytips(k), labels(k), ...
            "HorizontalAlignment", "center", ...
            "VerticalAlignment", "top")
    end
end


end

