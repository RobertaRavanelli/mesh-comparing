function [V,Vc,Vn,F3,linebyline]=readobj(fname,fpath)
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

%The OBJ parsing routine will work with ASCII-encoded files only. 
%It loads vertex coordinates and preserves the rest of the file unchanged.
%Vertex colors, which may be appended to the same line as the vertex coordinates,
%are stored in a separate matrix Vc for writing back into the file.
%2D or 4D vertices (x-y-z-w) are not supported. It should work with
%multiple meshes. Faces are loaded as well, but only those which are
%triangulated. Quads will be noted. If the file does not have // formatting
%for the faces (ie, vertex references only) the program will fail.
%If a parallel pool is open, it will make use of the workers.

fid = fopen(fullfile(fpath,fname),'r');
if fid<0
    error(['Cannot open ' fname '.']);
end

l = textscan(fid,'%s','delimiter','\n');
linebyline = l{1};
fprintf('Read %d lines from %s...', length(linebyline),fname);
fclose(fid);

vertex_lines = linebyline(strncmp('v ',linebyline,2),:);
V = zeros(length(vertex_lines), 3);
Vc = zeros(length(vertex_lines), 3);
%Check if there are 3D vertex coordinates, sometimes followed by vertex color RGB values
[~, lenVs] = sscanf(vertex_lines{1}, '%*c%f');

useparfor = false;
paralleljob = gcp('nocreate');

if ~isempty(paralleljob)
    useparfor = true;
    
    if lenVs == 3
        Vc = [];
        parfor i = 1:length(vertex_lines)
            v = sscanf(char(vertex_lines{i}), 'v %f %f %f');
            V(i,:) = v';
        end
    elseif lenVs == 6
        %Vc = zeros(length(vertex_lines), 3);
        parfor i = 1:length(vertex_lines)
            v = sscanf(char(vertex_lines{i}), 'v %f %f %f %f %f %f');
            V(i, :) = v(1:3)';
            Vc(i, :) = v(4:6)';
        end
    else
        error(' no 3D vertex data encountered');
    end
else
    if lenVs == 3
        Vc = [];
        for i = 1:length(vertex_lines)
            v = sscanf(char(vertex_lines{i}), 'v %f %f %f');
            V(i, :) = v';
        end
    elseif lenVs == 6
        %Vc = zeros(length(vertex_lines), 3);
        for i = 1:length(vertex_lines)
            v = sscanf(char(vertex_lines{i}), 'v %f %f %f %f %f %f');
            V(i, :) = v(1:3)';
            Vc(i, :) = v(4:6)';
        end
    else
        error(' no 3D vertex data encountered');
    end
end
clear vertex_lines;

normal_lines = linebyline(strncmp('vn ',linebyline,2),:);
Vn = zeros(length(normal_lines),3);

if ~isempty(normal_lines)
    if useparfor
        parfor i = 1:length(normal_lines)
            vn = sscanf(char(normal_lines{i}), 'vn %f %f %f');
            Vn(i, :) = vn';
        end
    else
        for i = 1:length(normal_lines)
            vn = sscanf(char(normal_lines{i}), 'vn %f %f %f');
            Vn(i, :) = vn';
        end
    end
end
clear normal_lines;

face_lines = linebyline(strncmp('f ',linebyline,2),:);
F3lines = false(length(face_lines),1);
f4s = 0;

if useparfor
    parfor i = 1:length(face_lines)
        nverts = length(regexp(face_lines{i},'\s\d')); %a new vertex is a space followed by a number
        if nverts == 3, F3lines(i) = true;
        elseif nverts == 4, f4s = f4s + 1; end
    end
else
    for i = 1:length(face_lines)
        nverts = length(regexp(face_lines{i},'\s\d')); %a new vertex is a space followed by a number
        if nverts == 3, F3lines(i) = true;
        elseif nverts == 4, f4s = f4s + 1; end
    end
end

if sum(F3lines) > 0
    F3 = zeros(sum(F3lines),3);
    f3ace_lines = face_lines(F3lines);
    texpattern = 'f %u/%*s%u/%*s%u/%*c';
    nSlashes = regexp(f3ace_lines{1},'\w*/\w*', 'once');
    if isempty(nSlashes), texpattern = 'f %u %u %u*c'; end
    if useparfor
        parfor i = 1:length(f3ace_lines)
            f = sscanf(char(f3ace_lines{i}), texpattern);
            F3(i, :) = f'; %For a face point as v/vt/vn, ignore vt/vn
        end
    else
        for i = 1:length(f3ace_lines)
            f = sscanf(char(f3ace_lines{i}), texpattern);
            F3(i, :) = f'; %For a face point as v/vt/vn, ignore vt/vn
        end
    end
    if f4s > 0
        fprintf(' Mesh contains %d quad faces that have been ignored',f4s);
    end
    clear f3ace_lines;
elseif f4s > 0
    error(' Mesh contains only quad faces');
else
    error(' No recognizable face data encountered');
end
clear face_lines; clear F3lines;
fprintf('with %d vertices and %d triangulated faces\n', size(V,1),size(F3,1));
return