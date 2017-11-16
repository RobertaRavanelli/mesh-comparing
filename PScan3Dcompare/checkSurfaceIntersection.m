function intMatrix = checkSurfaceIntersection(V, vNorm, U1,U2,U3, uNorms, connectivityV)
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
    numU = size(uNorms,1);
    testIDX = 1:numU;
    intMatrix  = false(numU,1);
    zero = eps*10; %Deal with rounding errors
    
    %% === Stage 1 ==========================================================
    % Those triangles which share an edge (two points) need only be tested by their dot products
    connectedTris = cellfun(@(x) nnz(x),num2cell(connectivityV,2));
    edgeTris = find(connectedTris==2);
    vertTris = find(connectedTris==1);
    if ~isempty(edgeTris)
        adjNorms = dot_prod(vNorm,uNorms(edgeTris,:));
        edgeTris(adjNorms > -0.99) = [];
        if ~isempty(edgeTris), intMatrix(edgeTris) = true; end
    end
    if ~isempty(vertTris) %And those sharing one point need only be tested if the dot product of their normals is negative
        adjNorms = dot_prod(vNorm,uNorms(vertTris,:));
        testFace = find(adjNorms < 0);
        if ~isempty(testFace)
            %Check if the opposite edge of both triangles from the shared point
            %intersects with other face
            testTris = vertTris(testFace);
            sharedVIDX = nonzeros(transpose(connectivityV(testTris,:)));
            [sharedUIDX,~] = find(transpose(connectivityV(testTris,:)));
            testEdgesIDXv = otherDim(sharedVIDX);
            testEdgesIDXu = otherDim(sharedUIDX);
            Us = zeros(3,length(testTris),3); Us(1,:,:) = U1(testTris,:); Us(2,:,:) = U2(testTris,:); Us(3,:,:) = U3(testTris,:);
            for i=1:length(testTris)
                testEdgeU = transpose([squeeze(Us(testEdgesIDXu(i,1),i,:)),squeeze(Us(testEdgesIDXu(i,2),i,:))]);
                testEdgeV = [V(testEdgesIDXv(i,1),:);V(testEdgesIDXv(i,2),:)];
                iSect = edgeIsectTri(testEdgeU,V) || edgeIsectTri(testEdgeV,squeeze(Us(:,i,:)));
                if iSect, intMatrix(testTris(i)) = true; end
            end
        end
    end
    testIDX(connectedTris>0) = []; %Do not examine these faces again
    numU = length(testIDX);
    
    %% === Stage 2 ==========================================================
    % If all vertices of one triangle lie to one side of the other's plane, they cannot intersect.
    dV = zeros(numU,3); v1 = repmat(V(1,:),numU,1); v2 = repmat(V(2,:),numU,1); v3 = repmat(V(3,:),numU,1);
    dV(:,1) = dot_prod(v1-U1(testIDX,:),uNorms(testIDX,:));
    dV(:,2) = dot_prod(v2-U1(testIDX,:),uNorms(testIDX,:));
    dV(:,3) = dot_prod(v3-U1(testIDX,:),uNorms(testIDX,:));
    dV(abs(dV)<zero)=0;
    excludedIDX = dV(:,1).*dV(:,2)>0 & dV(:,1).*dV(:,3)>0;
    dV(excludedIDX,:) = [];
    testIDX(excludedIDX) = [];
    if ~all(excludedIDX)
        numU = length(testIDX);

        dU = zeros(numU,3); vNs = repmat(vNorm,numU,1); v1 = repmat(V(1,:),numU,1);
        dU(:,1) = dot_prod(U1(testIDX,:)-v1,vNs);
        dU(:,2) = dot_prod(U2(testIDX,:)-v1,vNs);
        dU(:,3) = dot_prod(U3(testIDX,:)-v1,vNs);
        dU(abs(dU)<zero) = 0;
        excludedIDX = dU(:,1).*dU(:,2)>0 & dU(:,1).*dU(:,3)>0;
        dU(excludedIDX,:) = []; dV(excludedIDX,:) = [];
        testIDX(excludedIDX) = [];
    end
    
    %% === Stage 3 ==========================================================
    % Test coplanar triangles for overlap in two passes.
    coplanarIDX = dV(:,1)==0 & dV(:,2)==0 & dV(:,3)==0;
    if any(coplanarIDX)
        cpTmp = testIDX(coplanarIDX);
        testIDX(coplanarIDX) = [];
        dU(coplanarIDX,:) = []; dV(coplanarIDX,:) = [];
        coplanarIDX = cpTmp;

        for i=1:length(coplanarIDX)
            overlap = true;
            for idim = 1:3
                v = [V(1,idim), V(2,idim), V(3,idim)];
                u = [U1(coplanarIDX(i),idim), U2(coplanarIDX(i),idim), U3(coplanarIDX(i),idim)];
                t1 = min(v,[],2); t2 = max(v,[],2);
                s1 = min(u,[],2); s2 = max(u,[],2);
                overlap = overlap & (s1<=t2 & t1<=s2);
            end
            if overlap %Coplanar triangle tests: edge-to-edge intersection testing
                intMatrix(coplanarIDX(i)) = TriangleIntersection2D(V, vNorm(coplanarIDX(i),:), U1(coplanarIDX(i),:), U2(coplanarIDX(i),:), U3(coplanarIDX(i),:));
            end
        end
    end
    
    %% === Stage 4 ==========================================================
    % Check for 3D intersections for the remaining non-coplanar triangles
    if ~isempty(testIDX)
        intMatrix(testIDX) = TriangleIntersection3D(V, dV, vNorm, U1(testIDX,:), U2(testIDX,:), U3(testIDX,:), dU, uNorms(testIDX,:));
    end
    
    intMatrix = find(intMatrix); %Just return the indices with intersections
end
    
%% ========================================================================
function iSect = edgeIsectTri(eTest,Tri)
    iSect = false;
    zero = eps*10;
    
    dir = eTest(2,:)-eTest(1,:);
    tE1 = Tri(3,:)-Tri(1,:);
    tE2 = Tri(2,:)-Tri(1,:);
    q = cross_prod(dir,tE2);
    a = dot_prod(tE1,q);
    if a>-zero && a<zero, return; end
    
    f = 1/a;
    s = eTest(1,:)-Tri(1,:);
    u = f*dot_prod(s,q);
    if u<0 || u>1, return; end
    
    r = cross_prod(s,tE1);
    v = f*dot_prod(dir,r);
    if v<0 || u+v>1, return; end
    
    iSect = true;
end

function iSect = TriangleIntersection3D(V, dv, vN, U1s, U2s, U3s, du, uN)
    function IDX = getIsolatedPoints(ds)
        Dsigns = num2cell(sign(ds),2);
        Dzeros = cellfun(@(x) 3-nnz(x),Dsigns);
        IDXtmp = cellfun(@(x) -sum(x),Dsigns);
        IDXtmp(Dzeros==2) = -1*IDXtmp(Dzeros==2);
        IDXzero = abs(IDXtmp)==2;
        IDXtmp(Dzeros==1) = 1;
        IDXtmp(IDXzero) = 0;
        IDX = cellfun(@(x,y) find(y==x),Dsigns,num2cell(IDXtmp));
    end
    
    ptProj = @(p1,d1,p2,d2) p1 + (p2-p1)*d1/(d1-d2);
    nEval = size(dv,1); iSect = false(nEval,1);
    vIDX = getIsolatedPoints(dv); uIDX = getIsolatedPoints(du);
    vOrder = [vIDX,otherDim(vIDX)]; uOrder = [uIDX,otherDim(uIDX)];
    
    for i=1:nEval
        Vi = V(vOrder(i,:),:); dvi = dv(i,vOrder(i,:));
        Ui = [U1s(i,:); U2s(i,:); U3s(i,:)];
        Ui = Ui(uOrder(i,:),:); dui = du(i,uOrder(i,:));
        
        Lvec = cross_prod(vN,uN(i,:));
        Lvec = Lvec/norm(Lvec);
        
        P12 = ptProj(Vi(1,:),dvi(1), Vi(2,:),dvi(2));
        P13 = ptProj(Vi(1,:),dvi(1), Vi(3,:),dvi(3));
        Q12 = ptProj(Ui(1,:),dui(1), Ui(2,:),dui(2));
        Q13 = ptProj(Ui(1,:),dui(1), Ui(3,:),dui(3));
        
        tP13 = dot_prod(P13-P12,Lvec);
        tQ12 = dot_prod(Q12-P12,Lvec);
        tQ13 = dot_prod(Q13-P12,Lvec);
        
        Ip = sort([0,tP13]); Iq = sort([tQ12,tQ13]);
        a = max([Ip(1),Iq(1)]); b = min([Ip(2),Iq(2)]);
        if a<b, iSect(i) = true; end
    end
end

%% ========================================================================
function overlap = TriangleIntersection2D(V,N,U1,U2,U3)
    function intersect = EdgesIntersect3D(V1,V2, U1,U2)
        A = V2-V1; B = U1-U2; C = U1-V1;
        % Solve system of equations [A,B,1] * [d;e;0] = C for d and e
        det3 = @(a,b)   a(:,1).*b(:,2)-a(:,3).*b(:,2) + a(:,2).*b(:,3)-a(:,2).*b(:,1) + a(:,3).*b(:,1)-a(:,1).*b(:,3);
        f=det3(A,B); %https://en.wikipedia.org/wiki/Cramer%27s_rule#Explicit_formulas_for_small_systems
        t=det3(C,B)./f; s=det3(A,C)./f;
        intersect = (t>=0 & t<=1 & s>=0 & s<=1);
    end
    function inside = PointInTriangle2D(V1, U) %check if V1 is inside triangle U (U1,U2,U3)
        det2 = @(A,B,C) (A(:,1)-C(:,1))*(B(:,2)-C(:,2)) - (B(:,1)-C(:,1))*(A(:,2)-C(:,2));
        b1 = (det2(U(1,:), U(2,:), V1) > 0); b2 = (det2(U(2,:), U(3,:), V1) > 0); b3 = (det2(U(3,:), U(1,:), V1) > 0);
        inside = ((b1 == b2) & (b2 == b3));
    end
    % Edge-edge intersections
    overlap = false(size(N,1),1);
    i1Idx = [1 1 1 2 2 2 3 3 3]; i2Idx = [3 3 3 1 1 1 2 2 2];
    j1Idx = [1 2 3 1 2 3 1 2 3]; j2Idx = [3 1 2 3 1 2 3 1 2];
    U = zeros(3,size(N,1),3);
    U(:,:,1) = U1; U(:,:,2) = U2; U(:,:,3) = U3;
    for row = 1:size(N,1)
        % When it is necesary to project 3D plane on 2D, dIdx will be the optimal dimensions to use.
        [~, a] = max(abs(N(row,:)));
        [b, c] = otherDim(a); 
        dIdx = [b, c];
        % triangles overlap if any edges of triangle 1 cross those of 2
        edgeMat = EdgesIntersect3D(squeeze(V(row,:,i1Idx))',squeeze(V(row,:,i2Idx))', ...
            squeeze(U(row,:,j1Idx))',squeeze(U(row,:,j2Idx))');
        overlap(row) = any(edgeMat);
        if ~overlap(row)
            % project onto an axis-aligned plane that maximizes the area of the triangles
            % compute indices: dIdx which correspond to 2 smallest N1 components
            V2d = [V(row,dIdx,1); V(row,dIdx,2); V(row,dIdx,3)]; % each row is a 2D vertex
            U2d = [U(row,dIdx,1); U(row,dIdx,2); U(row,dIdx,3)];
            % test if tri1 is totally contained in tri2 or vice varsa
            if PointInTriangle2D(V2d(1,:), U2d) % tri1 is totally contained in tri2
                overlap(row) = true;
            elseif PointInTriangle2D(U2d(1,:), V2d) % tri2 is totally contained in tri1
                overlap(row) = true;
            end
        end
    end
end

% ========================================================================
function CP = cross_prod(a,b), CP = [a(:,2).*b(:,3)-a(:,3).*b(:,2), a(:,3).*b(:,1)-a(:,1).*b(:,3), a(:,1).*b(:,2)-a(:,2).*b(:,1)]; end
function DP = dot_prod(a,b), DP = a(:,1).*b(:,1)+a(:,2).*b(:,2)+a(:,3).*b(:,3); end
function [bc] = otherDim(a)
    bc = zeros(length(a),2);
    for i=1:length(a)
        bc(i,1) = mod(a(i),3)+1;
        bc(i,2) = 6-a(i)-bc(i,1);
    end
end
