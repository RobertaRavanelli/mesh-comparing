function [distances] = estimatePointPlane(refFaces, refVerts, comparedVerts)
%% https://github.com/psapirstein/mesh-comparing
% This code is distributed under an Apache License 2.0
% Author: Philip Sapirstein, UNL
% However, see below the simplied code for triangle-ray distance
% adapted from a script by Shridhar Ravikumar.
%
% The subroutine supports the collection of tools for processing 3D meshes
% and assessing their repeatability accompanying the article:
% "A high-precision photogrammetric recording system for small artifacts"
% Philip Sapirstein, Journal of Cultural Heritage 2017
% https://doi.org/10.1016/j.culher.2017.10.011
%%
    %This returns signed distances for all points measured by point-to-triangle
    initL = fprintf('Initializing...');
    distsFace = cell(size(refFaces,1),1);
    ptsImproved = false(size(comparedVerts,1),1);
    nF3s = size(refFaces,1);
    
    %First, a preliminary measure of distances to establish the maximum
    %distance (just point to point, not using the upsampled reference)
    kdOBJ = KDTreeSearcher(refVerts);
    [minvertIDXref, distances] = knnsearch(kdOBJ,comparedVerts);
    clear kdOBJ;

    %Generate lookup table with index of floating vertices which are closest to each reference face
    %First, the lookup table from reference vertex to reference face
    refFacesByrefVertIDX=cell(size(refVerts,1),1);
    for i=1:size(refFaces,1)
        vIDX = refFaces(i,:);
        for j=1:3
            refFacesByrefVertIDX{vIDX(j)} = [refFacesByrefVertIDX{vIDX(j)}, i];
        end
    end
    %Next, the lookup table from reference face to floating vertex
    floatVertsIDXByrefFace=cell(size(refFaces,1),1);
    for i=1:size(comparedVerts,1)
        fIDX = refFacesByrefVertIDX{minvertIDXref(i)};
        for j=1:length(fIDX)
            floatVertsIDXByrefFace{fIDX(j)} = [floatVertsIDXByrefFace{fIDX(j)}, i];
        end
    end
    
    fprintf([repmat(char(8),1,initL),'Measuring point-to-face distances:  ']);
    strLenTimer = 1;
    startTime = clock;
    chunk = 1500; %Divide up the data for parallel processing
    for ch=0:floor(nF3s/chunk)
        firstIDX = ch*chunk;
        lastIDX = firstIDX+chunk-1;
        if lastIDX>(nF3s-1), lastIDX = nF3s-1; end
        chLen = 1+lastIDX-firstIDX;
        
        CHrefFaceVs=cell(chLen,1);
        CHrefVerts=cell(chLen,1);
        CHdists=cell(chLen,1);
        for i=1:chLen
            CHrefFaceVs{i} = refVerts(refFaces(firstIDX+i,:),:);
            cIDX = floatVertsIDXByrefFace{firstIDX+i};
            if ~isempty(cIDX)
                CHrefVerts{i} = comparedVerts(cIDX,:);
                CHdists{i} = distances(cIDX);
            end
        end
        parfor i=1:chLen
            if ~isnan(CHrefVerts{i})
                distsFace{firstIDX+i} = getPtTriDist(CHrefFaceVs{i},CHrefVerts{i},CHdists{i});
            end
        end
        fprintf(repmat(char(8),1,strLenTimer));
        duration = etime(clock,startTime);
        percProgress = (ch+1)*chunk/nF3s;
        sRemaining = duration/percProgress-duration;
        strLenTimer = fprintf('%d%% complete (%ds remaining)',uint8(100*percProgress),uint16(sRemaining));
    end
    fprintf([repmat(char(8),1,strLenTimer),'Finalizing. ']);
    for i=1:nF3s %Determine minimum distances for each vertex
        cIDX = floatVertsIDXByrefFace{i};
        if ~isempty(cIDX)
            for j=1:length(cIDX)
                if abs(distsFace{i}(j)) < abs(distances(cIDX(j)))
                %Save the closer projected point coordinate and its distance from the floating vertex
                    distances(cIDX(j)) = distsFace{i}(j);
                    ptsImproved(cIDX(j)) = true;
                else %Keep the original distance, but adopt the sign
                    distances(cIDX(j)) = sign(distsFace{i}(j))*distances(cIDX(j));
                end
            end
        end
    end
    if sum(ptsImproved)/size(comparedVerts,1) < 0.99
        fprintf('Warning: %.1f%% of points were measured by point-to-point distances\n',100*(1-sum(ptsImproved)/size(comparedVerts,1)));
    end
end

%%
function [Distances] = getPtTriDist(Triangle,Points,KnownDistances)
%**************************************************************************
% Modified from code by Shridhar Ravikumar, University of Bath based on Mark W. Jones' paper
% '3D Distance from a Point to a Triangle.' URL - http://www-compsci.swan.ac.uk/~csmark/PDFS/dist.pdf
%**************************************************************************
    Distances = Inf(size(Points,1),1);
    
    %Get triangle vertices
    A = Triangle(1,:); B = Triangle(2,:); C = Triangle(3,:);
    %Translate vertex B and C to the origin
    B = B - A;      C = C - A;
    
    % Find theta which rotates vertex B on the X axis to be on XZ plane by setting y coordinate to 0
    % and rotate vertex B and C, so that B is on the XZ plane
    Theta = atan(B(2)/B(3));
    RotX = [1 , 0 , 0; 0, cos(Theta), -sin(Theta); 0, sin(Theta), cos(Theta)];
    B = RotX * B';  C = RotX * C';
    
    % Find the Phi which rotates B on the Y axis to be on Z axis, by setting the x coordinate to 0
    Phi = atan(-B(1)/B(3));
    RotY = [cos(Phi), 0, sin(Phi); 0, 1, 0; -sin(Phi), 0 , cos(Phi)];
    B = RotY * B;   C = RotY * C;
	
    %Find Gamma which rotates vertex C on the Z axis to be on YZ plane, by setting the x coordinate to 0
    Gamma = atan(C(1)/C(2));
    RotZ = [cos(Gamma), -sin(Gamma), 0; sin(Gamma), cos(Gamma), 0; 0, 0 ,1];
    C = RotZ * C;
    
    %Apply the transformation to the points, and create the projected matrix on the plane of the triangle
    rotXYZ = RotZ*RotY*RotX;
    Pts = Points - repmat(A,size(Points,1),1);
    Pts = rotXYZ * Pts';
    
    A = [0;0;0];
    AB = B - A;
    BC = C - B;
    CA = A - C;
    
    E = NaN(3,1);
    for i=1:size(Pts,2)
    %Only examine points likely to be less than the known maximum distance, returning Inf otherwise
        if abs(Pts(1,i)) < KnownDistances(i)
            % Figure out which side of the triangle the points are on
            ProjPt = [0; Pts(2,i); Pts(3,i)];
            E(1) = CheckPointSideOf2DLine(ProjPt,AB,A);
            E(2) = CheckPointSideOf2DLine(ProjPt,BC,B);
            E(3) = CheckPointSideOf2DLine(ProjPt,CA,C);
            if all(E<0) || all(E>0) || any(E==0) %Point is inside triangle or on an edge
                Distances(i) = Pts(1,i);
            else %The point is outside the triangle
                ClosestPt = ClosestTriangleEdgePointToPoint(A,B,C,ProjPt);
                Distances(i) = norm(Pts(:,i) - ClosestPt);
                if Pts(1,i) < 0
                    Distances(i) = -abs(Distances(i));
                end
            end
        end
    end
end

%%
function E = CheckPointSideOf2DLine(CheckPoint,DirectionOfLine,PointOnLine)
    % It's assumed that points and lines being passed are on YZ plane.
    xPt = CheckPoint(2);
    yPt = CheckPoint(3);
    xLn = PointOnLine(2);
    yLn = PointOnLine(3);
    dX = DirectionOfLine(2);
    dY = DirectionOfLine(3);
    E = (xPt - xLn)*dY - (yPt-yLn)*dX;
end

%%
function ClosestPointOnTriangle = ClosestTriangleEdgePointToPoint(A,B,C,OutsidePoint)
% Assuming the projected point lies outside the triangle and the points are all on the YZ plane
    AB = B - A;
    BC = C - B;
    % Check if the point on the left/right side of the first edge (AB)
    pointSide = CheckPointSideOf2DLine(OutsidePoint,AB,B);    
    thirdVertexSide = CheckPointSideOf2DLine(C,AB,B);
    % If the projected point and the third vertex of the triangle lie on
    % the same side of this line, then this line is definitely not the
    % closest line to the point. But if on opposite sides, this is the closest edge.
    if sign(pointSide) ~= sign(thirdVertexSide)
        ClosestPointOnTriangle = ClosestPointOnSegmentToPoint(OutsidePoint,A,B);
    else
        % Try the second edge (BC)
        pointSide = CheckPointSideOf2DLine(OutsidePoint,BC,B);    
        thirdVertexSide = CheckPointSideOf2DLine(A,BC,B);
        if sign(pointSide) ~= sign(thirdVertexSide)
            ClosestPointOnTriangle = ClosestPointOnSegmentToPoint(OutsidePoint,B,C);
        else % It must be closest to edge (CA)
            ClosestPointOnTriangle = ClosestPointOnSegmentToPoint(OutsidePoint,C,A);
        end
    end
end

%%
function ClosestPointOnSegment = ClosestPointOnSegmentToPoint(OutsidePoint, SegmentPoint1, SegmentPoint2)
    % Rotate the side by 90 degrees to get the normal direction
    Theta = degtorad(90); cosTheta = cos(Theta); sinTheta = sin(Theta);
    DirectionOfSegment = SegmentPoint2 - SegmentPoint1;
    x = DirectionOfSegment(2);
    y = DirectionOfSegment(3);
    
    NormalToSide = [x*cosTheta - y*sinTheta; x*sinTheta + y*cosTheta];
    NormalToSide = NormalToSide/norm(NormalToSide);
    
    normalSide1 = CheckPointSideOf2DLine(OutsidePoint,[0, NormalToSide(1), NormalToSide(2)],SegmentPoint1);
    normalSide2 = CheckPointSideOf2DLine(OutsidePoint,[0, NormalToSide(1), NormalToSide(2)],SegmentPoint2);
    
    if(sign(normalSide1) ~= sign(normalSide2) )
        % Point lies in between A and B and outside of the triangle, so it is closest to line AB
        ClosestPointOnSegment = ClosestPointOnLineToPoint(OutsidePoint,SegmentPoint1,SegmentPoint2);
    else
        %Point is close to one of the vertices on the line
        if(norm(OutsidePoint - SegmentPoint1) < norm(OutsidePoint - SegmentPoint2))
            ClosestPointOnSegment = SegmentPoint1;
        else
            ClosestPointOnSegment = SegmentPoint2;
        end
    end
end
   
%%
function ClosestPointOnLine = ClosestPointOnLineToPoint(P,A,B)
    AToB = B - A;
    AToB = AToB/ norm(AToB);
    AToP = P - A;
    
    ProjectionOfAToPOnAToB =  AToB' * AToP;
    ClosestPointOnLine = A + AToB * ProjectionOfAToPOnAToB;    
end