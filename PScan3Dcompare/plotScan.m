function plotScan(scan,cmap,dosaturation,withlights)
%% https://github.com/psapirstein/mesh-comparing
% This code is distributed under an Apache License 2.0
% Author: Philip Sapirstein, UNL
%
% The subroutine supports the collection of tools for processing 3D meshes
% and assessing their repeatability accompanying the article:
% "A high-precision photogrammetric recording system for small artifacts"
% Philip Sapirstein, Journal of Cultural Heritage 2017
% https://doi.org/10.1016/j.culher.2017.10.011
%%
    ts = get(gca,'children');
    if dosaturation, cmap = resaturate(cmap,2); end
    if isempty(ts)
        gcf; ts = patch;
        ts.Vertices = scan.Vs; ts.Faces = scan.F3s;
        ts.LineStyle = 'none'; ts.FaceColor = 'interp';
        ts.FaceAlpha = 0.9; ts.FaceLighting = 'gouraud';
        ts.FaceVertexCData = cmap;
        colormap(parula); colorbar; axis equal;
        ts.AmbientStrength = 0.5;
        ts.SpecularStrength = 0.3; ts.SpecularExponent = 20;
        lighting none;
        lightangle(30,45); lightangle(-75,75); lightangle(75,-60);
    else
        lighting none;
        for i=1:length(ts)
            if strcmp(ts(i).Type,'patch')
                ts(i).Vertices = scan.Vs;
                ts(i).Faces = scan.F3s;
                ts(i).FaceVertexCData = cmap;
            end
        end
    end
    caxis([min(min(cmap,[],2)), max(max(cmap,[],2))]);
    if withlights, lighting gouraud; end
end

function x = resaturate(x,tau)
    if nargin<2, tau = 1; end
    x = x-mean(x(:));
    x = x/(2*mean(abs(x)));
    x = max(x,-tau);
    x = min(x,tau);
end