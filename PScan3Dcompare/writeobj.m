function writeobj(modelname,V,Vc,linebyline)
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
    % This will save a copy to an ASCII formatted file with new vertex coordinates
    % (assumes an 8-decimal precision for vertex and 6-decimal color values).
    % Because it is assumed that the model has been repositioned during the
    % scripts calling this function, the vertex normals are discarded.
    fid = fopen(modelname,'w');
    vidx = 1;
    Vs = cell(size(V,1),1);
    
    paralleljob = gcp('nocreate');
    if ~isempty(paralleljob)
        if Vc == -1 %Save vertex coordinates only; -1 flags no color data
            parfor i = 1:size(V,1)
                Vs{i} = sprintf('v %.8f %.8f %.8f', V(i,1), V(i,2), V(i,3)); %#ok<PFBNS>
            end
        elseif size(Vc,2) == 1 %1-D color array indicates vertex quality is being saved
            maxVc = max(Vc); %Save vertex quality as an RGB triplet (Quality 1/good = 0.5 green; Quality 0/bad = 1.0 red);
            parfor i = 1:size(V,1)
                Vs{i} = sprintf('v %.8f %.8f %.8f %.6f %.6f %.6f', V(i,1), V(i,2), V(i,3), (maxVc-Vc(i))/maxVc, 0.5*Vc(i)/maxVc, 0); %#ok<PFBNS>
            end
        else %Save true vertex color data passed as a 3-dimensional RGB array
            parfor i = 1:size(V,1)
                Vs{i} = sprintf('v %.8f %.8f %.8f %.6f %.6f %.6f', V(i,1), V(i,2), V(i,3), Vc(i,1), Vc(i,2), Vc(i,3)); %#ok<PFBNS>
            end
        end
    else
        if Vc == -1 %Save vertex coordinates only; -1 flags no color data
            for i = 1:size(V,1)
                Vs{i} = sprintf('v %.8f %.8f %.8f', V(i,1), V(i,2), V(i,3));
            end
        elseif size(Vc,2) == 1 %1-D color array indicates vertex quality is being saved
            maxVc = max(Vc); %Save vertex quality as an RGB triplet (Quality 1/good = 0.5 green; Quality 0/bad = 1.0 red);
            for i = 1:size(V,1)
                Vs{i} = sprintf('v %.8f %.8f %.8f %.6f %.6f %.6f', V(i,1), V(i,2), V(i,3), (maxVc-Vc(i))/maxVc, 0.5*Vc(i)/maxVc, 0);
            end
        else %Save true vertex color data passed as a 3-dimensional RGB array
            for i = 1:size(V,1)
                Vs{i} = sprintf('v %.8f %.8f %.8f %.6f %.6f %.6f', V(i,1), V(i,2), V(i,3), Vc(i,1), Vc(i,2), Vc(i,3));
            end
        end
    end

    for i = 1:size(linebyline,1)
        if strncmp('v ',linebyline{i},2)
            fprintf(fid,'%s\n',Vs{vidx});
            vidx = vidx+1;
        elseif strncmp('vn ',linebyline{i},2)
            %Discard vertex normals
        else
            fprintf(fid,'%s\n',linebyline{i});
        end
    end

    fclose(fid);
end