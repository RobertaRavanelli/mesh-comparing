function [updatedTriangulation,Vs,F3s,deletedVs] = repairMeshTriangulation(objScan)
%% https://github.com/psapirstein/mesh-comparing
% This code is distributed under an Apache License 2.0
% Author: Philip Sapirstein, UNL
%
% The subroutine supports the collection of tools for processing 3D meshes
% and assessing their repeatability accompanying the article:
% "A high-precision photogrammetric recording system for small artifacts"
% Philip Sapirstein, Journal of Cultural Heritage 2017
% https://doi.org/10.1016/j.culher.2017.10.011
%%  %Merge duplicated or very close vertices and flip edges at sharp triangulations
    %NB: this function only works with watertight meshes (no holes are allowed or it will fail)
    dot_prod = @(a,b) a(:,1).*b(:,1)+a(:,2).*b(:,2)+a(:,3).*b(:,3);
    updatedTriangulation = 0; deletedVs = [];
    Vs = objScan.Vs; F3s = objScan.F3s;
    
    tic();
    fprintf('Unifying normals: '); ncharDel = fprintf('initializing');
    fIDXbyVert = getFaceIndices(size(Vs,1),F3s);
    [fIDXbyEdge,badIDX] = getAdjacentFaces(F3s,fIDXbyVert);
    fprintf(repmat(char(8),1,ncharDel));
    if ~isempty(badIDX)
        fprintf('Stopping due to bad triangulations (non-manifold vertex or holes) at the following faces:\n');
        fprintf('%d ',badIDX);
        updatedTriangulation = -1;
        return;
    end
    
    [F3s,fIDXbyEdge,nFlip,isolatedCs] = unifyNormals(F3s,fIDXbyEdge);
    if nFlip>0, fprintf(' Flipped %d faces. ',nFlip); updatedTriangulation=1;
    else, fprintf(' Consistent triangulation. '); end
    if isolatedCs>0
        fprintf('Found %d separate clusters; cannot unify normals. ',isolatedCs);
        fprintf('\nDelete these clusters before continuing.');
        updatedTriangulation = -1;
        return;
    else %Ensure that the normals are pointing outward
        maxIDX = zeros(3,2); normD(:,:,1) = [-1,0,0;0,-1,0;0,0,-1]; normD(:,:,2) = [1,0,0;0,1,0;0,0,1];
        [~,maxIDX(1,1)] = min(Vs(:,1)); [~,maxIDX(1,2)] = max(Vs(:,1));
        [~,maxIDX(2,1)] = min(Vs(:,2)); [~,maxIDX(2,2)] = max(Vs(:,2));
        [~,maxIDX(3,1)] = min(Vs(:,3)); [~,maxIDX(3,2)] = max(Vs(:,3));
        vNorms = vertexNormal(triangulation(F3s,Vs));
        dotSum = sum(dot_prod(vNorms(maxIDX(:,1),:),normD(:,:,1))) + sum(dot_prod(vNorms(maxIDX(:,2),:),normD(:,:,2)));
        ncharDel = 0;
        if dotSum < 0
            ncharDel = fprintf('Inverting normals');
            F3s = reverseNormals(F3s);
            fIDXbyVert = getFaceIndices(size(Vs,1),F3s);
            [fIDXbyEdge,~] = getAdjacentFaces(F3s,fIDXbyVert);
            updatedTriangulation=1;
        end
        fprintf(repmat(char(8),1,ncharDel));
    end
    fprintf('finished. ');
    
    nVs = size(Vs,1);
    nF3s = size(F3s,1);
    toc();
    
    fprintf('Examining triangulation to repair mesh.\n');
    %Find the closest points in the reference scan, and then find the set of
    %points nearest every vertex of the reference.
    kdOBJ = KDTreeSearcher(Vs);
    nearestPts = knnsearch(kdOBJ,Vs,'K',20);

    edgeThreshShort = 0.025; %Constant: triangles with very short edges will have their vertices merged
    edgeThreshLong = 0.975; %Constant: slivers created by incorrect triangulations
    
    lastChanged = 0;
    iterations = 0; maxIts = 10; minIts = 4;
    while iterations < maxIts
    %%  Several iterations may be needed, although most of the calculation will be spent on the first pass
        iterations = iterations + 1;
        fprintf('Pass %d: ',iterations);
        closeVs = cell(0); mergeVs = cell(0);
        exchangeFIDX = [];
        
        edgeVecs = zeros(size(F3s,1),3,3); edgeLens = zeros(size(F3s));
        for j=1:3
            edgeVecs(:,:,j)=Vs(F3s(:,mod(j,3)+1),:)-Vs(F3s(:,j),:);
            edgeLens(:,j) = sqrt(sum(edgeVecs(:,:,j).^2,2));
        end
        faceNormals = cross(edgeVecs(:,:,1),edgeVecs(:,:,2),2);
        faceNormals = normalizeRows(faceNormals);
        zeroFNs = find(all(faceNormals==0,2));
        if ~isempty(zeroFNs)
            for j=1:length(zeroFNs)
                eVtmp = squeeze(edgeVecs(zeroFNs(j),:,:));
                eVtmp(eVtmp==0) = -1;
                fN2 = cross(eVtmp(:,1),eVtmp(:,2))';
                faceNormals(zeroFNs(j),:) = normalizeRows(fN2);
            end
        end
        FNdots = [dot_prod(faceNormals,faceNormals(fIDXbyEdge(:,1),:)), ...
            dot_prod(faceNormals,faceNormals(fIDXbyEdge(:,2),:)), ...
            dot_prod(faceNormals,faceNormals(fIDXbyEdge(:,3),:))];
        zeroFNs = find(all(faceNormals==0,2));
        
    %%  Isolate short edges or sliver-like faces
        divs = 10;
        notifyStep = floor(nF3s/divs);
        for i=1:nF3s
            [fEdgeLens,edgeIDX] = sort(edgeLens(i,:));
            %In the first pass, close/overlying vertices are removed. There
            %should be no reason to do this iteratively.
            fixCloseVertices = (fEdgeLens(1) < fEdgeLens(2)*edgeThreshShort) & iterations == 1;
            zeroAreaFace = all(faceNormals(i,:)==0) & iterations == 1;
            
            %If any neighboring normals are nearly perpendicular or opposed, try to repair the face triangulation
            fixSliverFace = fEdgeLens(3) > sum(fEdgeLens(1:2))*edgeThreshLong; %Do so even if they are not flipped, but there is a very long edge
            if mod(iterations,minIts)==1, fixFlippedFace = sum(FNdots(i,:) < 0.1)>1;
            else, fixFlippedFace = any(FNdots(i,:) < 0.1); end
            
            if fixCloseVertices %Merge the two points by averaging their position.
                vIDXclose = [mod(edgeIDX(1),3)+1, edgeIDX(1)];
                oldVs = ismember(F3s(i,vIDXclose),reshape(cell2mat(mergeVs),[],1));
                if ~any(oldVs)
                    mergeVs{end+1} = F3s(i,vIDXclose); %#ok<AGROW>
                    meanPt = mean(Vs(F3s(i,vIDXclose),:));
                    Vs(F3s(i,vIDXclose(1)),:) = meanPt;
                    Vs(F3s(i,vIDXclose(2)),:) = meanPt;
                    closeVs{end+1} = i; %#ok<AGROW>
                end
            elseif zeroAreaFace %Caused by all edges of the triangle being parallel
                %Triangles of this sort must always be deleted in the first pass
                vIDXshort = 5-(edgeIDX(3)+mod(edgeIDX(3),3)); %Delete the point on the shorter edges
                mergeVs{end+1} = F3s(i,[mod(vIDXshort,3)+1,vIDXshort]); %#ok<AGROW>
                closeVs{end+1} = i; %#ok<AGROW>
            elseif fixSliverFace || fixFlippedFace
                %Flip edges repair poor triangulations with skinny faces or
                %faces whose normals reverse those in its area
                nextFace = fIDXbyEdge(i,edgeIDX(3));
                vIDXlongEdge = [mod(edgeIDX(3),3)+1, edgeIDX(3)];
                
                nearFs = unique(reshape([fIDXbyVert{unique(reshape(nearestPts(F3s(i,:),:),[],1))}],[],1));
                nearFs(nearFs==i)=[];
                nearFs(ismembc(nearFs,zeroFNs))=[];
                cFs = F3s(nearFs,:);
                connectivityFs = zeros(size(cFs));
                for j=1:3, connectivityFs(cFs==F3s(i,j))=j; end
                oldiSects = checkSurfaceIntersection(Vs(F3s(i,:),:), faceNormals(i,:), ...
                    Vs(cFs(:,1),:),Vs(cFs(:,2),:),Vs(cFs(:,3),:),faceNormals(nearFs,:),connectivityFs);
                %Check whether there were any self-intersecting elements in the nearby triangulation

                
                if iterations < minIts %For the first passes, just flip edges with the neighboring triangle    
                    newFace = true;
                    newFs = sort([i,nextFace]);
                    if sum(ismembc(newFs,unique(reshape(exchangeFIDX,1,[])))), newFace = false; end
                    if newFace
                        newV = setxor(F3s(i,vIDXlongEdge),F3s(nextFace,:));
                        altTriangulation = repmat(F3s(i,:),2,1);
                        altTriangulation(1,vIDXlongEdge(1)) = newV;
                        altTriangulation(2,vIDXlongEdge(2)) = newV;
                        %Normals for these new face triangulations
                        edgeVecsNew = zeros(2,3,3);
                        for j=1:3, edgeVecsNew(:,:,j) = bsxfun(@minus,Vs(altTriangulation(:,mod(j,3)+1),:),Vs(altTriangulation(:,j),:)); end
                        fNnew = cross(edgeVecsNew(:,:,1),edgeVecsNew(:,:,2),2);
                        fNnew = normalizeRows(fNnew);
                        zeroFNs = find(all(fNnew==0,2));
                        if ~isempty(zeroFNs)
                            for j=1:length(zeroFNs)
                                eVtmp = squeeze(edgeVecsNew(zeroFNs(j),:,:));
                                eVtmp(eVtmp==0) = -1;
                                fN2 = cross(eVtmp(:,1),eVtmp(:,2))';
                                fNnew(zeroFNs(j),:) = normalizeRows(fN2);
                            end
                        end
                        if any(all(fNnew==0,2)), continue; end %No improvement possible due to a zero-area face triangulation
                
                        nearFs(nearFs==nextFace)=[];
                        newiSects = cell(0);
                        for iTests=1:2
                            cFs = [altTriangulation(mod(iTests,2)+1,:);F3s(nearFs,:)];
                            connectivityFs = zeros(size(cFs));
                            for j=1:3, connectivityFs(cFs==altTriangulation(iTests,j))=j; end
                            newiSects{iTests} = checkSurfaceIntersection(Vs(altTriangulation(iTests,:),:), ...
                                fNnew(iTests,:), Vs(cFs(:,1),:),Vs(cFs(:,2),:),Vs(cFs(:,3),:), ...
                                [fNnew(mod(iTests,2)+1,:);faceNormals(nearFs,:)],connectivityFs);
                        end
                        newiSects = cellfun(@(x) ~isempty(x),newiSects);
                        
                        %Test consistency of the normals and face areas
                        fRing = unique(reshape(fIDXbyEdge(newFs,:),[],1));
                        fRing = fRing(~ismembc(fRing,newFs));
                        fRingVec = mean(faceNormals(fRing,:),1);
                        newDots = (2*dot_prod(fNnew(1,:),fNnew(2,:)) + dot_prod(mean(fNnew,1),fRingVec))/3;
                        oldDots = (2*dot_prod(faceNormals(i,:),faceNormals(nextFace,:)) + ...
                            dot_prod(mean([faceNormals(i,:);faceNormals(nextFace,:)],1),fRingVec))/3;                        
                        fAreasNew = cellfun(@(x) norm(x),num2cell(cross(edgeVecsNew(:,:,1),edgeVecsNew(:,:,2),2),2))/2;
                        fAreasOld = cellfun(@(x) norm(x),num2cell(cross(edgeVecs(newFs,:,1),edgeVecs(newFs,:,2)),2))/2;
                        aRatio = (max(fAreasOld)-mean(fAreasOld))/(max(fAreasNew)-mean(fAreasNew));
                        
                        improvedTriangulation = (1-newDots) < aRatio*(1-oldDots)/2;
                        originalTriangulation = sum(cellfun(@(x,y) isequal(x,y),num2cell(altTriangulation,2),num2cell(objScan.F3s([i,nextFace],:),2)));
                        if sum(newiSects) < length(oldiSects)
                            improvedTriangulation = true;
                            originalTriangulation = false;
                        end %Always try the new triangulation if intersecting faces are removed
                        %Never apply the new triangulation if it creates new self-intersections
                        if sum(newiSects) > length(oldiSects), continue; end
                        if improvedTriangulation && ~originalTriangulation
                            oldFs = F3s([i,nextFace],:); %Preserve original face ordering in case it must be reversed
                            F3s(i,:) = altTriangulation(1,:);
                            F3s(nextFace,:) = altTriangulation(2,:);
                            [fIDXbyVert,changedTriangulation] = updateFaceIndices(fIDXbyVert,F3s,newFs);
                            if changedTriangulation
                                exchangeFIDX = [exchangeFIDX; newFs]; %#ok<AGROW>
                                fIDXbyEdge = updateEdgeIndices(F3s,fIDXbyVert,fIDXbyEdge,newFs);
                                for j=1:3
                                    edgeVecs(i,:,j)=edgeVecsNew(1,:,j); edgeVecs(nextFace,:,j)=edgeVecsNew(2,:,j);
                                    edgeLens(i,j) = sqrt(sum(edgeVecsNew(1,:,j).^2,2)); edgeLens(nextFace,j) = sqrt(sum(edgeVecsNew(2,:,j).^2,2));
                                end
                                faceNormals(i,:) = fNnew(1,:); faceNormals(nextFace,:) = fNnew(2,:);
                            else %Reverse changes if some problem in the triangulation
                                %prevents edge flipping (such as a non-manifold point)
                                F3s(i,:) = oldFs(1,:);
                                F3s(nextFace,:) = oldFs(2,:);
                            end
                        end
                    end  
                else %For later passes, flip edges in a 3-triangle network
                    %Find a third face attached to the longest edge of the adjacent triangle
                    [~,edgeIDX2] = sort(edgeLens(nextFace,:));
                    thirdFace = fIDXbyEdge(nextFace,edgeIDX2);
                    thirdFace(thirdFace==i) = []; %Delete the original face i from the list
                    if length(intersect(F3s(i,:),F3s(thirdFace(end),:))) > 1
                        %If the third face shares an edge with the first use the other edge
                        thirdFace = thirdFace(1);
                    else
                        thirdFace = thirdFace(end);
                    end
                    
                    newFace = true;
                    newFs = sort([i,nextFace,thirdFace]);
                    if sum(ismembc(newFs,unique(reshape(exchangeFIDX,1,[])))), newFace = false; end
                    if newFace
                        Vring = zeros(5,1);
                        %The resulting points will form a ring p1-p2-p3-p4-p5
                        %where p1 = opposite point face i, p2 / p5 are on its long edge,
                        %p3 is unique to the third face, and p4 unique to the intermediate face
                        Vring(1) = setdiff(F3s(i,:),F3s(nextFace,:));
                        Vring(2) = intersect(F3s(i,:),F3s(thirdFace,:));
                        Vring(3) = setdiff(F3s(thirdFace,:),F3s(nextFace,:));
                        Vring(4) = setdiff(F3s(nextFace,:),F3s(i,:));
                        Vring(5) = setdiff(F3s(i,vIDXlongEdge),Vring(2));
                        
                        altTriangulation = zeros(3,3,2);
                        altTriangulation(:,:,1) = [Vring(1),Vring(3),Vring(2); ...
                            Vring(1),Vring(5),Vring(3); Vring(5),Vring(4),Vring(3)];
                        altTriangulation(:,:,2) = [Vring(4),Vring(1),Vring(5); ...
                            Vring(1),Vring(4),Vring(2); Vring(2),Vring(4),Vring(3)];
                        %These triangulations should maintain the coorect order because p1/p2 set the normal direction,
                        %and they are from the original face i. Likewise with p5/p1 for the next triangle and p3-4-5 on the third.
                        %Try two possible configurations with the new faces radiating from either p3 or p4
                        fRing = unique(reshape(fIDXbyEdge(newFs,:),[],1));
                        fRing = fRing(~ismembc(fRing,newFs));
                        fRingVec = mean(faceNormals(fRing,:),1);
                        edgeVecsNew = zeros(3,3,3,2); fNnew = zeros(3,3,2);
                        scorePass = zeros(2,1); iSectsPass = zeros(2,1);
                        nearFs(nearFs==nextFace)=[];
                        avfNold = mean([dot_prod(faceNormals(i,:),faceNormals(nextFace,:)), ...
                            dot_prod(faceNormals(nextFace,:),faceNormals(thirdFace,:))]);
                        oldDots = (2*avfNold + dot_prod(mean(faceNormals(newFs,:),1),fRingVec))/3;
                        fAreasOld = cellfun(@(x) norm(x),num2cell(cross(edgeVecs(newFs,:,1),edgeVecs(newFs,:,2)),2))/2;
                        for k=1:2
                            for j=1:3, edgeVecsNew(:,:,j,k) = bsxfun(@minus,Vs(altTriangulation(:,mod(j,3)+1,k),:),Vs(altTriangulation(:,j,k),:)); end
                            fNnew(:,:,k) = cross(edgeVecsNew(:,:,1,k),edgeVecsNew(:,:,2,k),2);
                            fNnew(:,:,k) = normalizeRows(fNnew(:,:,k));
                            zeroFNs = find(all(squeeze(fNnew(:,:,k))==0,2));
                            if ~isempty(zeroFNs)
                                for j=1:length(zeroFNs)
                                    eVtmp = squeeze(edgeVecsNew(zeroFNs(j),:,:,k));
                                    eVtmp(eVtmp==0) = -1;
                                    fN2 = cross(eVtmp(:,1),eVtmp(:,2))';
                                    fNnew(zeroFNs(j),:,k) = normalizeRows(fN2);
                                end
                            end
                            
                            avfNnew = mean([dot_prod(fNnew(1,:,k),fNnew(2,:,k)),dot_prod(fNnew(2,:,k),fNnew(3,:,k))]);
                            newDots = (2*avfNnew + dot_prod(mean(fNnew(:,:,k),1),fRingVec))/3;
                            fAreasNew = cellfun(@(x) norm(x),num2cell(cross(edgeVecsNew(:,:,1),edgeVecsNew(:,:,2),2),2))/2;
                            aRatio = (max(fAreasOld)-mean(fAreasOld))/(max(fAreasNew)-mean(fAreasNew));
                            scorePass(k) = (1-newDots)/aRatio;
                            
                            if any(all(squeeze(fNnew(:,:,k))==0,2))
                                iSectsPass(k) = Inf; %Reject the rare cases where a zero-area face is produced
                            else
                                passiSects = cell(0);
                                for iTests=1:3
                                    otherIDX = [mod(iTests,3)+1,mod(iTests+1,3)+1];
                                    cFs = [altTriangulation(otherIDX,:,k);F3s(nearFs,:)];
                                    connectivityFs = zeros(size(cFs));
                                    for j=1:3, connectivityFs(cFs==altTriangulation(iTests,j,k))=j; end
                                    passiSects{iTests} = checkSurfaceIntersection(Vs(altTriangulation(iTests,:,k),:), ...
                                        fNnew(iTests,:,k), Vs(cFs(:,1),:),Vs(cFs(:,2),:),Vs(cFs(:,3),:), ...
                                        [fNnew(otherIDX,:,k);faceNormals(nearFs,:)],connectivityFs);
                                end
                                iSectsPass(k) = sum(cellfun(@(x) ~isempty(x),passiSects));
                            end
                        end
                        
                        [newScore,mintriIDX] = min(scorePass);
                        newiSects = iSectsPass(mintriIDX);
                        if iSectsPass(1)~=iSectsPass(2)
                            [newiSects,mintriIDX] = min(iSectsPass);
                            newScore = scorePass(mintriIDX);
                        end %Triangulations with intersecting faces are always rejected
                        edgeVecsNew = squeeze(edgeVecsNew(:,:,:,mintriIDX));
                        altTriangulation = squeeze(altTriangulation(:,:,mintriIDX));
                        fNnew = squeeze(fNnew(:,:,mintriIDX));
                        if any(all(fNnew==0,2)), continue; end
                        
                        oldFs = F3s([i,nextFace,thirdFace],:); %Preserve original face ordering
                        improvedTriangulation = newScore < (1-oldDots)/2;
                        originalTriangulation = sum(cellfun(@(x,y) isequal(x,y),num2cell(altTriangulation,2),num2cell(objScan.F3s([i,nextFace,thirdFace],:),2)));
                        if newiSects < length(oldiSects)
                            improvedTriangulation = true;
                            originalTriangulation = false;
                        end
                        if newiSects > length(oldiSects), continue; end
                        if improvedTriangulation && ~originalTriangulation
                            F3s(i,:) = altTriangulation(1,:);
                            F3s(nextFace,:) = altTriangulation(2,:);
                            F3s(thirdFace,:) = altTriangulation(3,:);
                            [fIDXbyVert,changedTriangulation] = updateFaceIndices(fIDXbyVert,F3s,newFs);
                            
                            if changedTriangulation
                                exchangeFIDX = [exchangeFIDX; newFs]; %#ok<AGROW>
                                fIDXbyEdge = updateEdgeIndices(F3s,fIDXbyVert,fIDXbyEdge,newFs);
                                for j=1:3
                                    edgeVecs(i,:,j) = edgeVecsNew(1,:,j);
                                    edgeVecs(nextFace,:,j) = edgeVecsNew(2,:,j);
                                    edgeVecs(thirdFace,:,j) = edgeVecsNew(3,:,j);
                                    edgeLens(i,j) = sqrt(sum(edgeVecsNew(1,:,j).^2,2));
                                    edgeLens(nextFace,j) = sqrt(sum(edgeVecsNew(2,:,j).^2,2));
                                    edgeLens(thirdFace,j) = sqrt(sum(edgeVecsNew(3,:,j).^2,2));
                                end
                                faceNormals(i,:) = fNnew(1,:);
                                faceNormals(nextFace,:) = fNnew(2,:);
                                faceNormals(thirdFace,:) = fNnew(3,:);
                            else %Reverse temporary changes to the triangulation
                                F3s(i,:)=oldFs(1,:);
                                F3s(nextFace,:)=oldFs(2,:);
                                F3s(thirdFace,:)=oldFs(3,:);
                            end
                        end
                    end
                end
            end
            if mod(i,notifyStep)==0, fprintf('*'); end
        end        
        fprintf(repmat(char(8),1,divs));
        
        %% Update the triangulation if points will be deleted
        if ~isempty(mergeVs)
            mergeVs = sort(transpose(reshape(cell2mat(mergeVs),2,[])),1);
            for i=1:size(mergeVs,1)
                %Replace the second vertex with the first duplicate vertex in the face matrix
                fIDX = fIDXbyVert{mergeVs(i,2)};
                for j=1:length(fIDX)
                    F3row = F3s(fIDX(j),:);
                    F3row(F3row==mergeVs(i,2)) = mergeVs(i,1);
                    F3s(fIDX(j),:) = F3row;
                end
            end
            deletedVs = sort(mergeVs(:,2));
            
            newVidx = -ones(nVs,1); delIDX=1; shift = 0;
            for i=1:nVs
                if i==deletedVs(delIDX)
                    delIDX = min(delIDX+1,length(deletedVs));
                    shift = shift-1;
                else
                    newVidx(i) = i+shift;
                end
            end
            %Delete the obsolete vertices (the second matching vertex)
            Vs(deletedVs,:) = [];
            
            notifyStep = floor(nF3s/10); ncharDel = 0;
            deletedF3s = cell(0);
            for i=1:nF3s
                F3row = newVidx(F3s(i,:))';
                F3s(i,:) = F3row;
                if length(unique(F3row)) < 3 || any(F3row == -1)
                    deletedF3s{end+1} = i; %#ok<AGROW>
                end %There is probably a faster way to calculate this
                if mod(i,notifyStep)==0
                    fprintf('*');
                    ncharDel = ncharDel+1;
                end
            end
            F3s(cell2mat(deletedF3s),:) = [];
            %Refresh lookup tables whenever any vertices are deleted
            nVs = size(Vs,1); nF3s = size(F3s,1);
            fIDXbyVert = getFaceIndices(nVs,F3s);
            [fIDXbyEdge,~] = getAdjacentFaces(F3s,fIDXbyVert);
            kdOBJ = KDTreeSearcher(Vs);
            nearestPts = knnsearch(kdOBJ,Vs,'K',20);
            fprintf(repmat(char(8),1,ncharDel));
        end
%%      Report the results
        if ~isempty(closeVs) || ~isempty(exchangeFIDX)
            updatedTriangulation = 1;
            if ~isempty(closeVs)
                nVsCut = size(objScan.Vs,1)-nVs;
                nF3sCut = size(objScan.F3s,1)-nF3s;
                vertstr = 'vertices';
                if nVsCut == 1
                    vertstr = 'vertex';
                end
                fprintf('Merged %d nearby %s and %d faces. ',nVsCut,vertstr,nF3sCut);
            end
            if ~isempty(exchangeFIDX)
                justChanged = size(exchangeFIDX,1) + (iterations > minIts)*size(exchangeFIDX,1);
                fprintf('Flipped %d edges of slender or badly oriented faces.',justChanged);
                if (justChanged == lastChanged) && isempty(closeVs) && iterations > minIts
                    fprintf(' No further edge-flipping improvements appear to be possible.\n');
                    return;
                else
                    lastChanged = justChanged;
                end
            end
            fprintf('\n');
        else
            if iterations==1
                fprintf('No changes needed.\n');
            else
                fprintf('No further changes needed.\n');
            end
            return;
        end
    end
end

function [nA] = normalizeRows(A)
    nTmp = sqrt( sum( A.^2, 2 ) );
    nTmp(nTmp == 0) = 1;
    nA = bsxfun(@rdivide, A, nTmp);
end

function [F3sMod,fIDXedgeMod,nFlipped,isolatedClusters] = unifyNormals(F3sIn,fIDXedgeIn)
%%Unify the face normals of a scan with an untrusted vertex ordering
    F3sMod = F3sIn;
    nFlipped = 0;
    matchedFaces = false(size(F3sMod,1),1);
    matchedFaces(1) = true; %Seed process with first face
    unmatchedEdges = true(size(F3sMod));
    fIDXedgeMod = fIDXedgeIn;
    isolatedClusters = 0;
    
    divs = 10; step = ceil(size(F3sIn,1)/divs); notify = size(F3sIn,1):-step:step;
    while ~all(matchedFaces)
        matchFset = find(matchedFaces & any(unmatchedEdges,2));
        if isempty(matchFset)
            unmatchedIDX = find(matchedFaces==false);
            matchedFaces(unmatchedIDX(1)) = true;
            isolatedClusters = isolatedClusters+1;
        end
        for i=1:length(matchFset)
            fComp = matchFset(i);
            for j=find(unmatchedEdges(fComp,:))
                e1 = [F3sMod(fComp,mod(j,3)+1), F3sMod(fComp,j)];
                fComp2 = fIDXedgeMod(fComp,j);
                attachedVs = F3sMod(fComp2,:);
                ept2IDX = find(attachedVs==e1(1));
                ept1IDX = mod(ept2IDX,3)+1;
                if attachedVs(ept1IDX)~=e1(2)
                    %The points must be flipped
                    nFlipped = nFlipped+1;
                    ept1IDX = 1+mod(ept2IDX-2,3);
                    F3sMod(fComp2,ept1IDX)=e1(1);
                    F3sMod(fComp2,ept2IDX)=e1(2);
                    %Reset the matching on other edges
                    unmatchedEdges(fComp2,:) = [true,true,true];
                    flipIDX = find(fIDXedgeMod(fComp2,:)~=fComp);
                    f1 = fIDXedgeMod(fComp2,flipIDX(1));
                    fIDXedgeMod(fComp2,flipIDX(1)) = fIDXedgeMod(fComp2,flipIDX(2));
                    fIDXedgeMod(fComp2,flipIDX(2)) = f1;
                end
                unmatchedEdges(fComp,j) = false;
                unmatchedEdges(fComp2,ept2IDX) = false;
                matchedFaces(fComp2) = true;
            end
        end
        nMatched = sum(matchedFaces);
        if nMatched > notify(end)
            nPassed = false(length(notify));
            for i=1:length(notify)
                if nMatched>notify(i)
                    nPassed(i) = true;
                    fprintf('*');
                end
            end
            notify(nPassed)=[];
        end
    end
    fprintf(repmat(char(8),1,divs-1));
end

function [flipF3s] = reverseNormals(Fs)
    flipF3s = NaN(size(Fs));
    for iflip=1:size(Fs,1)
        flipF3s(iflip,:) = Fs(iflip,3:-1:1);
    end
end

function [connectedFs,badIDX] = getAdjacentFaces(F3s,fIDXbyV)
%%  Service function: compile index of faces by attached edge
    connectedFs = NaN(size(F3s));
    paralleljob = gcp('nocreate');
    badIDX = false(size(F3s,1),1);
    if isempty(paralleljob)
        for i=1:size(F3s,1)
            eFs = getAttachedFace(i,fIDXbyV(F3s(i,:)));
            if length(eFs)~=3, eFs = [0,0,0]; badIDX(i)=true; end
            connectedFs(i,:) = eFs;
        end
    else
        parfor i=1:size(F3s,1)
            eFs = getAttachedFace(i,fIDXbyV(F3s(i,:))); %#ok<PFBNS>
            if length(eFs)~=3, eFs = [0,0,0]; badIDX(i)=true; end
            connectedFs(i,:) = eFs;
        end
    end
    badIDX = find(badIDX);
end

function [connectedFs] = getAttachedFace(face,fIDX)
    fIDXtrim = cellfun(@(x) x(x~=face),fIDX,'UniformOutput',false);
    connectedFs = [fIDXtrim{2}(ismembc(fIDXtrim{2},fIDXtrim{1})), ...
        fIDXtrim{3}(ismembc(fIDXtrim{3},fIDXtrim{2})), ...
        fIDXtrim{1}(ismembc(fIDXtrim{1},fIDXtrim{3}))];
end

function [eIDXmod] = updateEdgeIndices(F3s,fIDX,eIDX,newIDX)
    eIDXmod = eIDX;
    for i=1:length(newIDX), eIDXmod(newIDX(i),:) = getAttachedFace(newIDX(i),fIDX(F3s(newIDX(i),:))); end
    adjacentFaces = unique(reshape(eIDXmod(newIDX,:),[],1));
    adjacentFaces = adjacentFaces(~ismembc(adjacentFaces,newIDX));
    for i=1:length(adjacentFaces)
        eIDXmod(adjacentFaces(i),:) = getAttachedFace(adjacentFaces(i),fIDX(F3s(adjacentFaces(i),:)));
    end
end

function [fIDX] = getFaceIndices(nVs,F3s)
%% Service function: get list of faces by each vertex
    fIDX = cell(nVs,1);
    for i=1:size(F3s,1)
        vIDX = F3s(i,:);
        for j=1:length(vIDX)
            fIDX{vIDX(j)}(end+1) = i;
        end
    end
    %NB: the face indices will always be sorted
end

function [fIDXmod,changed] = updateFaceIndices(fIDX,F3s,newIDX)
    fIDXmod = fIDX;
    changed = true;
    %First, erase the old face reference by vertex
    affectedVs = unique(reshape(F3s(newIDX,:),[],1));
    for i=1:length(affectedVs)
        fID = fIDXmod{affectedVs(i)};
        fIDXmod{affectedVs(i)} = fID(~ismembc(fID,newIDX));
    end
    %Next, add the new face references to each vertex
    for i=1:length(newIDX)
        vIDX = F3s(newIDX(i),:);
        for j=1:length(vIDX)
            fIDXmod{vIDX(j)} = unique([fIDXmod{vIDX(j)},newIDX(i)]);
        end
    end
    %Next, ensure that the triangulation is valid (no vertices connected to
    %fewer than three faces)
    for i=1:length(affectedVs)
        if length(fIDXmod{affectedVs(i)}) < 3
            fIDXmod = fIDX;
            changed = false;
            return;
        end
    end
    %Finally, ensure that the new triangulation has no non-manifold edges
    for i=1:length(newIDX)
        if length(getAttachedFace(newIDX(i),fIDX(F3s(newIDX(i),:)))) ~= 3
            fIDXmod = fIDX;
            changed = false;
            return;
        end
    end
end
