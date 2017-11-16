function [projectedVs, Vquality, err2sigma] = createAveragedModel(floatScan, refScan, averageResults)
%% https://github.com/psapirstein/mesh-comparing
% This code is distributed under an Apache License 2.0
% Author: Philip Sapirstein, UNL
% However, see below the simplied code for triangle-ray intersection
% adapted from a script by Jarek Tuszynski.
%
% The subroutine supports the collection of tools for processing 3D meshes
% and assessing their repeatability accompanying the article:
% "A high-precision photogrammetric recording system for small artifacts"
% Philip Sapirstein, Journal of Cultural Heritage 2017
% https://doi.org/10.1016/j.culher.2017.10.011
%%
    %Create an averaged model from two inputs; the floating scan's
    %points and triangulation are modified for the resulting model.
    if averageResults
        tic();
        Tolerance = max(refScan.maxXYZ)/200000; %Set below the lowest conceivable precision of the model
        PyramidLevels = 6; Visualize = 0; Iterations = 2; %It is assumed the scans are already positioned very close to one another
        %Reorient the "floating" scan to the reference with Tolga's ICP
        FinalPose = Tolga_icp_mod_point_plane_pyr(floatScan.Vs, floatScan.Vns, refScan.Vs, refScan.Vns, Tolerance, Iterations, 3, 1, PyramidLevels, 0, Visualize);
        floatScan.Vs = Tolga_movepoints(FinalPose, floatScan.Vs);
        floatScan.maxXYZ = [max(floatScan.Vs(:,1))-min(floatScan.Vs(:,1)), ...
            max(floatScan.Vs(:,2))-min(floatScan.Vs(:,2)), max(floatScan.Vs(:,3))-min(floatScan.Vs(:,3))];
        floatScan.Vns = vertexNormal(triangulation(floatScan.F3s,floatScan.Vs));
        toc();
    end
    
    %Set up general variables for parallel processing
    Vct = size(floatScan.Vs,1);
    NptsToCompare = 20; %Start by finding faces connected to N vertices near the floating scan
    pass = 0; maxIter = 2; projectionSuccessful = false;
    
    %Generate lookup table with reference scan faces that contain a particular vertex
    refFaceIDXbyV=cell(size(refScan.Vs,1),1);
    for i=1:size(refScan.F3s,1)
        vIDX = refScan.F3s(i,:);
        for j=1:3, refFaceIDXbyV{vIDX(j)} = [refFaceIDXbyV{vIDX(j)}, i]; end
    end
    
    while ~projectionSuccessful && pass<maxIter
        tic();
        initL = fprintf('Initializing...');
        
        %Find the closest points in the reference scan, and then find the set of
        %points nearest every vertex of the reference.
        kdOBJ = KDTreeSearcher(refScan.Vs);
        closeIDXref = knnsearch(kdOBJ,floatScan.Vs,'K',NptsToCompare);
        clear kdOBJ;
        
        projectedPts = refScan.Vs(closeIDXref(:,1),:); %Begin with closest points in the reference scan
        projectedQs = refScan.Qs(closeIDXref(:,1)); %Save quality values for these nearest vertices
        Par4noisect = zeros(Vct,1);
        
        fprintf([repmat(char(8),1,initL),'Projecting float scan points on reference:  ']);
        strLenTimer = 1;
        startTime = clock;
        chunk = 1500; %Divide up the data for parallel processing
        for ch=0:floor(Vct/chunk)
            firstIDX = ch*chunk;
            lastIDX = firstIDX+chunk-1;
            if lastIDX>(Vct-1), lastIDX = Vct-1; end
            chLen = 1+lastIDX-firstIDX;

            CHfltVs=floatScan.Vs(firstIDX+1:lastIDX+1,:);
            CHfltVns=floatScan.Vns(firstIDX+1:lastIDX+1,:);
            
            CHrefFverts=cell(chLen,3); CHfltSmoothVns=zeros(chLen,3);
            for i=1:chLen
                %From the list of closest reference vertices to those of the floating
                %scan, find the attached reference faces & make this list unique.
                refFaces = unique(cell2mat(refFaceIDXbyV(closeIDXref(firstIDX+i,:))'));
                %Compile the reference vertices for each triangular face in this list
                CHrefFverts{i,1} = refScan.Vs(refScan.F3s(refFaces,1),:);
                CHrefFverts{i,2} = refScan.Vs(refScan.F3s(refFaces,2),:);
                CHrefFverts{i,3} = refScan.Vs(refScan.F3s(refFaces,3),:);
                %Create a "smooth" normal combining all the nearby reference vertex normals
                if averageResults
                    avNorm = mean(refScan.Vns(closeIDXref(firstIDX+i,:),:),1);
                    CHfltSmoothVns(i,:) = avNorm/norm(avNorm);
                end
            end
            parfor i=1:chLen
                projPt = TriRayIntersection(CHfltVs(i,:), CHfltVns(i,:), CHrefFverts{i,1}, CHrefFverts{i,2}, CHrefFverts{i,3}); %#ok<PFBNS>
                if averageResults
                    projPtSmooth = TriRayIntersection(CHfltVs(i,:), CHfltSmoothVns(i,:), CHrefFverts{i,1}, CHrefFverts{i,2}, CHrefFverts{i,3}); 
                    dP = dot(CHfltVns(i,:),CHfltSmoothVns(i,:));
                    if dP < 0.8, projectedQs(firstIDX+i) = projectedQs(firstIDX+i) / (1 + 3*(0.8-dP)); end
                    if isnan(projPt)
                        Par4noisect(firstIDX+i) = 0.5;
                        projectedQs(firstIDX+i) = projectedQs(firstIDX+i)/2; %Lower the quality for this point
                        if isnan(projPtSmooth)
                            Par4noisect(firstIDX+i) = 1;
                            projectedQs(firstIDX+i) = projectedQs(firstIDX+i)/4;
                        else
                            projectedPts(firstIDX+i,:) = projPtSmooth;
                        end
                    else
                        if isnan(projPtSmooth)
                            projectedPts(firstIDX+i,:) = projPt;
                            projectedQs(firstIDX+i) = projectedQs(firstIDX+i)/2;
                        else
                            projectedPts(firstIDX+i,:) = (projPt+projPtSmooth)/2;
                        end
                    end
                elseif ~isnan(projPt)
                    projectedPts(firstIDX+i,:) = projPt;
                end
            end
            fprintf(repmat(char(8),1,strLenTimer));
            duration = etime(clock,startTime);
            percProgress = (ch+1)*chunk/Vct;
            sRemaining = duration/percProgress-duration;
            strLenTimer = fprintf('%d%% complete (%ds remaining)',uint8(100*percProgress),uint16(sRemaining));
        end
        fprintf(repmat(char(8),1,strLenTimer));
        if sum(Par4noisect)/Vct > 0.005 && averageResults
            pass = pass+1;
            fprintf('\nWarning (pass %d): %.1f%% of the floating scan vertices were not projected onto the reference scan. ',pass,100*sum(Par4noisect)/Vct);
            NptsToCompare = NptsToCompare*5;
        else, projectionSuccessful = true; pass = maxIter;
        end
        toc();
    end
    clear -regexp ^CH; clear -regexp ^Par4;
    clear closeIDXref;
    
    if averageResults
        tic();
        fprintf('Averaging matched points among the two scans and partially smoothing the results. ');
        
        %Face lookup for the floating scan
        floatFaceIDXbyV=cell(size(floatScan.Vs,1),1);
        for i=1:size(floatScan.F3s,1)
            vIDX = floatScan.F3s(i,:);
            for j=1:3, floatFaceIDXbyV{vIDX(j)} = [floatFaceIDXbyV{vIDX(j)}, i]; end
        end
        
        projectedVs = NaN(Vct,3);
        for i = 1:Vct %Weighted average of the floating scan with its points projected onto the reference surfaces
            floatRingVs = unique(reshape(floatScan.F3s([floatFaceIDXbyV{i}],:),1,[]));
            refQ = (mean(projectedQs(floatRingVs)) + projectedQs(i))/2;
            projW = refQ/(refQ+floatScan.Qs(i));
            floatW = floatScan.Qs(i)/(refQ+floatScan.Qs(i));
            projectedVs(i,:) = projectedPts(i,:)*projW + floatScan.Vs(i,:)*floatW;
        end
        
        %Quality by relative change from original floatScan vertices
        Vquality = zeros(Vct,1);
        distances = NaN(Vct,1);
        for i=1:Vct, distances(i) = sqrt(sum((floatScan.Vs(i,:)-projectedVs(i,:)).^2)); end
        
        maxdist = max(distances);
        fac = maxdist - prctile(distances,5);
        for i=1:Vct
            Vqi = (maxdist - distances(i))/fac;
            if Vqi > 1, Vqi = 1; end
            Vquality(i) = (2*Vqi + projectedQs(i))/3;
        end
        
        %Average the lower quality averaged points with smoothed points
        smoothVs = meannorm_trismooth(projectedVs,floatScan.F3s);
        for i = 1:Vct
           if Vquality(i) < 1, projectedVs(i,:) = ...
               Vquality(i)*projectedVs(i,:) + (1-Vquality(i)) * smoothVs(i,:);
           end
        end
        
        err2sigma = prctile(distances,95.45) / mean(Vquality);
        toc();
    else
        Vquality = []; err2sigma = [];
        projectedVs = projectedPts;
    end
end

function iPt = TriRayIntersection(orig, dir, vertA, vertB, vertC)
%% Adapted and simplified from Jarek Tuszynski (jaroslaw.w.tuszynski@leidos.com) and:
%  *"Fast, minimum storage ray-triangle intersection". Tomas Moeller and
%    Ben Trumbore. Journal of Graphics Tools, 2(1):21--28, 1997.
%    http://www.graphics.cornell.edu/pubs/1997/MT97.pdf
%  * http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/raytri/
%  * http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/raytri/raytri.c
    eps = 1.0e-9;
    iPt = NaN; %dist = Inf;
    orig = repmat(orig,(size(vertA,1)),1);
    dir = repmat(dir,(size(vertA,1)),1);
    
    % Find faces parallel to the ray
    edge1 = vertB-vertA;          % find vectors for two edges sharing vertA
    edge2 = vertC-vertA;
    tvec = orig-vertA;          % vector from vertA to ray origin
    pvec = cross(dir, edge2,2); % begin calculating determinant - also used to calculate U parameter
    det = sum(edge1.*pvec,2);   % determinant of the matrix M = dot(edge1,pvec)
    angleOK = (abs(det)>eps);   % if determinant is near zero then ray lies in the plane of the triangle
    if all(~angleOK), return; end % if all parallel then no intersections
    
    % Calculate variables for useful line/triangle pairs
    det(~angleOK) = NaN;            % change to avoid division by zero
    u = sum(tvec.*pvec,2)./det;     % 1st barycentric coordinate
    v = NaN(size(u)); ds=v;
    ok = (angleOK & u>=0 & u<=1);   % mask
    % if all line/plane intersections are outside the triangle then no intersections
    if ~any(ok), return; end
    qvec = cross(tvec(ok,:), edge1(ok,:),2); % prepare to test V parameter
    v(ok,:) = sum(dir(ok,:).*qvec,2) ./ det(ok,:); % 2nd barycentric coordinate
    ds(ok,:) = sum(edge2(ok,:).*qvec,2)./det(ok,:);
    % test if line/plane intersection is within the triangle
    isect = (ok & v>=0 & u+v<=1);
    [isectIDX,~] = find(isect == true);

    % only return the closest intersection point if there is more than one
    if length(isectIDX) > 1
        [~,closestIDX] = min(abs(ds(isectIDX)));
        isectIDX = isectIDX(closestIDX);
    end
    if sum(isect) > 0
        iPt = vertA(isectIDX,:) + edge1(isectIDX,:).*repmat(u(isectIDX,1),1,3) ...
            + edge2(isectIDX,:).*repmat(v(isectIDX,1),1,3);
    end
end
