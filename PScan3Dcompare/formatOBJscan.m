function [OBJstruct] = formatOBJscan(filename, path, basename)
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
    %Load a scan from an OBJ file, center it, and return a struct with its attributes
    [Vs, Vcs, Vns, F3s, LbL] = readobj(filename,path);
    %Create the vertex normals if not saved in the file
    if isempty(Vns)
        Vns = vertexNormal(triangulation(F3s,Vs));
    end
    maxDims = [max(Vs(:,1))-min(Vs(:,1)), max(Vs(:,2))-min(Vs(:,2)), ...
        max(Vs(:,3))-min(Vs(:,3))]; %max values in [X, Y, Z]
    
    OBJstruct.filename = filename;
    OBJstruct.basename = basename;
    OBJstruct.shortname = filename((length(basename)+1):end-4);
    OBJstruct.Vs = Vs;
    OBJstruct.Vcs = Vcs; %Vertex color from readobj; -1 if not present
    OBJstruct.Vns = Vns;
    OBJstruct.Qs = ones(size(Vs,1),1); %Vertex quality, default 1
    OBJstruct.F3s = F3s;
    OBJstruct.LineByLine = LbL; %Original OBJ file lines, for saving modifications
    OBJstruct.maxXYZ = maxDims;
end
