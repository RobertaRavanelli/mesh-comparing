function [vNorms,bNorms,bcMean,bcGauss] = getCurvatureBezierPatch(objScan)
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
    eVecs = @(Vs,F3s,idx) Vs(F3s(:,mod(idx,3)+1),:)-Vs(F3s(:,idx),:);
    Vs = objScan.Vs; F3s = objScan.F3s;
    
    fNorms = cross(eVecs(Vs,F3s,1),eVecs(Vs,F3s,2));
    fNorms = normalizeRows(fNorms); %These will be needed to ensure the estimated normal points the right direction

    %Lookup index for faces by vertex, calculate face normals
    fIDXbyVert = getFacesbyVertex(size(Vs,1),F3s);
    %Higher quality normals estimates by weighted adjacent edge normals
    vNorms = getVertexNormals(Vs,fIDXbyVert,F3s,fNorms);
    [bNorms,bcMean,bcGauss] = getBezierPatchCurvatures(Vs,vNorms,fIDXbyVert,F3s);
end

function [bezierNormals,bezierCurvaturesMean,bezierCurvaturesGauss] = getBezierPatchCurvatures(VsIn,vertexNormals,fIDXbyVs,F3sIn)
%%  Starting with well-oriented normals and an index of neighboring vertices
% Implementing the approach by Anshuman Razdan and MyungSoo Bae,
% "Curvature Estimation Scheme for Triangle Meshes Using Biquadratic Bézier
% Patches" in Computer-Aided Design 37.14 (2005): 1481-91.
% On the bicubic surface fitting here, see Zhong Li, Brian Barsky, and
% Xiaogang Jin, "An effective third-order local fitting patch and its
% application" in IEEE International Conference on Shape Modeling and
% Applications (SMI) 2009, pp 7-14.
%%%%%%%%%%%%%%%%%%%%
    fprintf('Estimating curvature from bicubic bezier patches:  ');
    level = 3;
    nVs = size(VsIn,1);
    bezierNormals = zeros(nVs,3);
    bezierCurvaturesGauss = zeros(nVs,1);
    bezierCurvaturesMean = zeros(nVs,1);
    %Bernstein polynomials and their derivatives for the bicubic surface
    B03 = @(x) (1-x)^3; B13 = @(x) 3*(1-x)^2*x; B23 = @(x) 3*(1-x)*x^2; B33 = @(x) x^3;
    d1B03=@(x) -3*(1-x)^2; d1B13=@(x) 3-12*x+9*x^2; d1B23=@(x) 6*x-9*x^2; d1B33=@(x) 3*x^2;
    d2B03=@(x) 6*(1-x); d2B13=@(x) -12+18*x; d2B23=@(x) 6-18*x; d2B33=@(x) 6*x;
    ShapeM = ... %the shape optimizing matrix
         [1, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0; ...
          0,  1, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0; ...
          1,  0,  0,  0, -2,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0; ...
          0,  0,  0,  0,  1,  0,  0,  0, -2,  0,  0,  0,  1,  0,  0,  0; ...
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -2,  1,  0; ...
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -2,  1; ...
          0,  0,  0,  1,  0,  0,  0, -2,  0,  0,  0,  1,  0,  0,  0,  0; ...
          0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0, -2,  0,  0,  0,  1];
    strLenTimer = 1;
    startTime = clock;
    
    chunk = 1500; %Divide up the data for parallel processing
    for ch=0:floor(nVs/chunk)
        firstIDX = ch*chunk;
        lastIDX = firstIDX+chunk-1;
        if lastIDX>(nVs-1), lastIDX = nVs-1; end
        chLen = 1+lastIDX-firstIDX;
        
        %Gather the faces and the ring of points surrounding the vertex
        CHringPts = cell(chLen,1);
        CHvertNrm = vertexNormals(firstIDX+1:lastIDX+1,:);
        CHbezNorm = zeros(chLen,3);
        CHgaussC = zeros(chLen,1);
        CHmeanC = zeros(chLen,1);
        for i=1:chLen
            connectedVIDX = i+firstIDX;
            for j=1:level
                connectedVIDX = unique([connectedVIDX,reshape(F3sIn([fIDXbyVs{connectedVIDX}],:),1,[])]);
                connectedVIDX(connectedVIDX == i+firstIDX) = [];
            end
            %Translate the neighboring points so p(i) is at the origin
            CHringPts{i} = bsxfun(@minus,VsIn(connectedVIDX,:),VsIn(i+firstIDX,:));
        end
        parfor i=1:chLen
            ringPts = CHringPts{i};
            %Remove some of the farthest points from p(i), down to n (or fewer if there are no more points)
            pointsKept = 24;
            if size(ringPts,1) > pointsKept
                RD = sqrt(sum(ringPts.^2,2));
                ringPts = ringPts(RD <= prctile(RD,pointsKept*100/length(RD)),:);
            end
            [RD, sorted] = sort(sqrt(sum(ringPts.^2,2)));
            ringPts = ringPts(sorted,:);
            %Seed the initial U-V coordinate system for each vertex by crossing
            %the normal with the vector from the first attached vertex
            seedEdge = ringPts(end,:);
            localX = cross(CHvertNrm(i,:),seedEdge); %Put the farthest point at V=0 ...
            localX = localX/norm(localX);
            localY = cross(CHvertNrm(i,:),localX); %and U = negative (max value is negative)
            
            %Project the remaining points in UVW coordinates
            Rot = [localY; localX; CHvertNrm(i,:)];
            projPtsLocal = [0,0,0; transpose(Rot*ringPts')]; %Add the compared point to the front of the list
            
            %Find the UV coordinates for the fitting surface with u,v in [0,1]
            uShift = min(projPtsLocal(:,1));        vShift = min(projPtsLocal(:,2));
            uvprojPts = bsxfun(@minus,projPtsLocal(:,1:2),[uShift,vShift]);
            projUVScale = max(uvprojPts);
            uvprojPts = bsxfun(@rdivide,uvprojPts,projUVScale);
            
            %Solve the linear equation system Ax=B for the Bezier surface knots
            Buv = zeros(size(projPtsLocal,1),16); %The basis functions are next
            alpha = 0.1; %Weight drops by distance from p(i), but p(i) is also lower in weight
            distWeights = [1-alpha; (max(RD)/(1-alpha)-RD)/(max(RD)-min(RD))];
            distWeights = diag(distWeights);
            for j=1:size(Buv,1)
                uJ = uvprojPts(j,1); vJ = uvprojPts(j,2);
                colV = [B03(vJ);B13(vJ);B23(vJ);B33(vJ)]; rowU = [B03(uJ),B13(uJ),B23(uJ),B33(uJ)];
                Buv(j,:) = reshape(colV*rowU,1,[]);
            end
            %Append the shape optimizing matrix
            BuvW = [distWeights*(1-alpha)*Buv; alpha*ShapeM];
            Aobs = [distWeights*(1-alpha)*projPtsLocal; repmat([0,0,0],8,1)];
            coeffsBline = BuvW\Aobs;
            coeffsBij = reshape(coeffsBline,4,4,3);
            
            %Evaluate the tangent and curvature vectors at p(i)
            ptU = uvprojPts(1,1); ptV = uvprojPts(1,2);
            colVs = [B03(ptV);B13(ptV);B23(ptV);B33(ptV)]; rowUs = [B03(ptU),B13(ptU),B23(ptU),B33(ptU)];
            colD1Vs =[d1B03(ptV);d1B13(ptV);d1B23(ptV);d1B33(ptV)]; rowD1Us =[d1B03(ptU),d1B13(ptU),d1B23(ptU),d1B33(ptU)];
            colD2Vs =[d2B03(ptV);d2B13(ptV);d2B23(ptV);d2B33(ptV)]; rowD2Us =[d2B03(ptU),d2B13(ptU),d2B23(ptU),d2B33(ptU)];
            d1Uvec = zeros(3,1); d1Vvec = zeros(3,1); d2UUvec = zeros(3,1); d2VVvec = zeros(3,1); d2UVvec = zeros(3,1);
            for j=1:3
                d1Uvec(j) = sum(sum(squeeze(coeffsBij(:,:,j)) .* (colVs*rowD1Us)));
                d1Vvec(j) = sum(sum(squeeze(coeffsBij(:,:,j)) .* (colD1Vs*rowUs)));
                d2UUvec(j) = sum(sum(squeeze(coeffsBij(:,:,j)) .* (colVs*rowD2Us)));
                d2VVvec(j) = sum(sum(squeeze(coeffsBij(:,:,j)) .* (colD2Vs*rowUs)));
                d2UVvec(j) = sum(sum(squeeze(coeffsBij(:,:,j)) .* (colD1Vs*rowD1Us)));
            end
            d1Wvec = cross(d1Uvec,d1Vvec);
            d1Wvec = d1Wvec/norm(d1Wvec);
            
            %Coefficients for the first and second fundamental forms
            Euu = dot(d1Uvec,d1Uvec); Fuv = dot(d1Uvec,d1Vvec); Gvv = dot(d1Vvec,d1Vvec);
            Luu = dot(d2UUvec,d1Wvec); Muv = dot(d2UVvec,d1Wvec); Nvv = dot(d2VVvec,d1Wvec);
            
            %Rotate the surface normal back into local coordinates
            CHbezNorm(i,:) = transpose(Rot\d1Wvec);
            CHgaussC(i) = (Luu*Nvv-Muv^2)/(Euu*Gvv-Fuv^2);
            CHmeanC(i) = (Nvv*Euu-2*Muv*Fuv+Luu*Gvv)/(2*(Euu*Gvv-Fuv^2));
        end
        bezierNormals(firstIDX+1:lastIDX+1,:) = CHbezNorm;
        bezierCurvaturesMean(firstIDX+1:lastIDX+1) = CHmeanC;
        bezierCurvaturesGauss(firstIDX+1:lastIDX+1) = CHgaussC;
        
        fprintf(repmat(char(8),1,strLenTimer));
        duration = etime(clock,startTime);
        percProgress = (ch+1)*chunk/nVs;
        sRemaining = duration/percProgress-duration;
        strLenTimer = fprintf('%d%% complete (%ds remaining)',floor(100*percProgress),floor(sRemaining));
    end
    fprintf([repmat(char(8),1,strLenTimer),'finished.\n']);
end

function [vertexNormals] = getVertexNormals(VsIn,fIDXbyVs,F3sIn,fNs)
%% An implementation of the normal vector algorithm by Hyoung-Seok Kim & Ho-Sook Kim,
%   "New Computation of Normal Vector and Curvature," in WSEAS Transactions
%   on Computers, 10.8 (October, 2009), pp. 1661-1670. While fairly
%   computationally intensive, the vectors are generally more accurate than
%   those generated by previous methods.
    fprintf('Calculating vertex normals:  ');
    strLenTimer = 1;
    startTime = clock;
    
    nVs = size(VsIn,1);
    vertexNormals = zeros(nVs,3);
    
    chunk = 1500; %Divide up the data for parallel processing
    for ch=0:floor(nVs/chunk)
        firstIDX = ch*chunk;
        lastIDX = firstIDX+chunk-1;
        if lastIDX>(nVs-1), lastIDX = nVs-1; end
        chLen = 1+lastIDX-firstIDX;
        
        %Gather the faces and the ring of points surrounding the vertex
        CHringF3IDX = fIDXbyVs(firstIDX+1:lastIDX+1);
        CHfaceNorms = zeros(chLen,3);
        CHringVIDX = cell(chLen,1); CHfaceVs = cell(chLen,1); CHringEdgeVectors = cell(chLen,1);
        for i=1:chLen
            CHringVIDX{i} = unique(reshape(F3sIn(CHringF3IDX{i},:),1,[]));
            CHringVIDX{i}(CHringVIDX{i}==firstIDX+i) = [];
            CHfaceVs{i} = F3sIn(CHringF3IDX{i},:);
            CHringEdgeVectors{i} = bsxfun(@minus,VsIn(CHringVIDX{i},:),VsIn(firstIDX+i,:));
            CHfaceNorms(i,:) = mean(fNs(CHringF3IDX{i},:),1);
        end
        parfor i=1:chLen
            %Now sort them so the points and faces circle the vertex in order
            lengthRing = length(CHringVIDX{i});
            sortedVs = zeros(lengthRing,1); sortedF3s = zeros(lengthRing,1);
            orderEdges = zeros(lengthRing,1);
            
            sortedF3s(end) = CHringF3IDX{i}(1);
            sortedVs(1) = CHfaceVs{i}(1,1);
            if sortedVs(1) == i+firstIDX
                sortedVs(1) = CHfaceVs{i}(1,2);
            end
            orderEdges(1) = find(CHringVIDX{i}==sortedVs(1));
            for j=2:lengthRing
                %Need to find the two faces that contain the last sorted
                %point and remove the face which has already been sorted
                nextF3IDX = sum(sortedVs(j-1)==CHfaceVs{i},2)';
                nextF3IDX = find(nextF3IDX & ~ismember(CHringF3IDX{i},sortedF3s));
                sortedVs(j) = CHfaceVs{i}(nextF3IDX, CHfaceVs{i}(nextF3IDX,:) ~= i+firstIDX & ...
                    CHfaceVs{i}(nextF3IDX,:) ~= sortedVs(j-1));
                sortedF3s(j-1) = CHringF3IDX{i}(nextF3IDX);
                orderEdges(j) = find(CHringVIDX{i}==sortedVs(j));
            end
            sortedEs = CHringEdgeVectors{i}(orderEdges,:);
            
            %Gather the edge normals and angles of each face at the point
            ringEdgeNormals = zeros(lengthRing,3);
            ringAngles = zeros(lengthRing,1);
            for j=1:length(sortedF3s)
                lastJ = 1+mod(j-2,lengthRing);
                nextJ = 1+mod(j,lengthRing);
                fE0 = sortedEs(lastJ,:); fE1 = sortedEs(j,:); fE2 = sortedEs(nextJ,:);
                fN1 = cross(fE1,fE0); fN1 = fN1/norm(fN1);
                fN2 = cross(fE2,fE1); fN2 = fN2/norm(fN2);
                ringEdgeNormals(j,:) = (fN1+fN2)/norm(fN1+fN2);
                ringAngles(j) = acos(dot(fE1 / norm(fE1), fE2 / norm(fE2)));
            end
            %Calculate omega and lamba for each point
            lamdas = zeros(lengthRing,1);
            for j=1:lengthRing
                lastJ = 1+mod(j-2,lengthRing);
                lamdas(j) = (tan(ringAngles(lastJ)/2)+tan(ringAngles(j)/2)) ...
                    / norm(CHringEdgeVectors{i}(j,:));
            end
            lamdas = lamdas/sum(lamdas);
            lamdaEdgeNormalsum = sum(bsxfun(@times,ringEdgeNormals,lamdas));
            vNorm = lamdaEdgeNormalsum/norm(lamdaEdgeNormalsum);
            if dot(vNorm,CHfaceNorms(i,:)) < 0, vertexNormals(firstIDX+i,:) = -vNorm;
            else, vertexNormals(firstIDX+i,:) = vNorm; end
        end
        fprintf(repmat(char(8),1,strLenTimer));
        duration = etime(clock,startTime);
        percProgress = (ch+1)*chunk/nVs;
        sRemaining = duration/percProgress-duration;
        strLenTimer = fprintf('%d%% complete (%ds remaining)',floor(100*percProgress),floor(sRemaining));
    end
    fprintf([repmat(char(8),1,strLenTimer),'finished.\n']);
end

function [fIDX] = getFacesbyVertex(Vct,F3sIn)
%% Build a lookup table of the face triangulation network
    fIDX = cell(Vct,1);
    for i=1:size(F3sIn,1)
        vIDX = F3sIn(i,:);
        for j=1:length(vIDX), fIDX{vIDX(j)} = [fIDX{vIDX(j)},i]; end
    end
end

function [nA] = normalizeRows(A)
    nTmp = sqrt( sum( A.^2, 2 ) );
    nTmp(nTmp == 0) = 1;
    nA = bsxfun(@rdivide, A, nTmp);
end
