function [void] = simmaptree_plot(simMapString,params)
    %Plot SimMap tree in "filename" with lineages colored.

tr = phytree_read(simMapString); %phytree_read has been modified to parse lineage annotations in SimMap format
tr = struct(tr);
if(isfield(params,'colorMap'))
    tr.colorMap = params.colorMap;
end
if(isfield(params,'colorBar'))
    tr.colorBar = params.colorBar;
else
    tr.colorBar = false;
end
phytree_plot(tr); %phytree_plot has been modified to plot lineage line segments


end

