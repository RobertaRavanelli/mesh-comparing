function [distance,cGaussRough1,cGaussRough2] = Wang_perceptual_distance(scn1, scn2)
%% This code has been copied from another script
%
% It has been modified for compatibility with the mesh comparison tools
% (https://github.com/psapirstein/mesh-comparing) accompanying:
% "A high-precision photogrammetric recording system for small artifacts"
% Philip Sapirstein, Journal of Cultural Heritage 2017
% https://doi.org/10.1016/j.culher.2017.10.011
%
% The Apache License 2.0 for the rest of the code base does not apply to
% this subroutine, which has been included here for convenience.
%%
% This is the main function for the fast mesh perceptual distance (FMPD)
% computation. FMPD is a roughness-based mesh visual quality metric.
% cf. our paper 'A fast roughness-based approach to the assessment of 3D 
% mesh visual quality', Computers & Graphics 2012
% Authors: Kai Wang, Fakhri Torkhani, and Annick Montanvert
% GIPSA-lab, AGPIG team, CNRS UMR5216, Grenoble, France
% Email: kai.wang@gipsa-lab.grenoble-inp.fr
% Note: We used in our code some functions from the toolboxes of
% Gabriel Payre (cf. http://gpeyre.github.io/numerical-tours/)
% Thanks to Gabriel for sharing his code.
% Usage: distance = perceptual_distance(name_firstmesh,name_secondmesh)
% distance: Fast Mesh Perceptual Distance (FMPD, between 0 and 1)

    %starttime = tic;
    % compute the necessary quantites for the FMPD computation, including the
    % Gaussian curvature, Laplacian matrix and mass matrix
    [cGauss1, Lapl1, Mass1] = compute_perception_quantities(scn1.Vs',scn1.F3s');
    [cGauss2, Lapl2, Mass2] = compute_perception_quantities(scn2.Vs',scn2.F3s');
    SurfA = [sum(spdiags(Mass1)), sum(spdiags(Mass2))];
    
    % compute the roughness as the Laplacian of Gaussian curvature
    [cGaussRough1,cGaussRmean(1)] = compute_roughness(cGauss1,Lapl1,Mass1,SurfA(1));
    [cGaussRough2,cGaussRmean(2)] = compute_roughness(cGauss2,Lapl2,Mass2,SurfA(2));
    % power model for modulating the roughness
    minmaxrough = [0.0005, max(0.20,5.0 * min(cGaussRmean))];
    threshold = min(cGaussRmean); threshold = max(threshold,minmaxrough(1));
    
    cGaussRoughMod1 = modulate_roughness(cGaussRough1,minmaxrough,threshold);
    cGaussRoughMod2 = modulate_roughness(cGaussRough2,minmaxrough,threshold);
    
    % global roughness
    sum_first = cGaussRoughMod1' * spdiags(Mass1,0) / SurfA(1);
    sum_second = cGaussRoughMod2' * spdiags(Mass2,0) / SurfA(2);
    % the scaling factor is just for bringing the distance into [0,1] interval
    c = 8.0; distance = c * abs(sum_first-sum_second);
    if distance>1, distance = 1; end
    
    % end of the processing time counting
    %elapsedtime = toc(starttime);
    %disp(strcat(['The processing time of FMPD is ' num2str(elapsedtime) ' seconds.']));
end

function [cGRout] = modulate_roughness(cGaussRough,MMrough,threshold)
    a = 0.15; b = 0.5; epsilon = MMrough(1);
    cGRout = clamp(cGaussRough,epsilon,MMrough(2));
    % then we modulate the roughness with the power model
    cGRout = (cGRout).^a - (epsilon).^a;
    turnpoint = (threshold).^a - (epsilon).^a;
    cGRout(cGRout>turnpoint) = (cGRout(cGRout>turnpoint)-turnpoint)*b + turnpoint;
end

function [cgauss, laplacian, mass] = compute_perception_quantities(Vs,F3s)
% Usage: [cgauss, laplacian, mass] = compute_perception_quantities(vertices,faces)
% cgauss: Gaussian curvature
% laplacian: Laplacian matrix
% mass: mass matrix
% vertices: the vertex coordinates
% faces: a matrix of the indices of the facet component vertices

    % number of vertices
    n = size(Vs,2);
    % number of facets
    m = size(F3s,2);

    % to store the angles sum on vertices
    sum_angles = zeros(n,1);
    % to store the facet areas
    area = zeros(m,1);
    % for numerical computation stability
    epsilon = 1e-10;

    for i=1:3
       i1 = mod(i-1,3)+1;
       i2 = mod(i  ,3)+1;
       i3 = mod(i+1,3)+1;
       pp = Vs(:,F3s(i2,:)) - Vs(:,F3s(i1,:));
       qq = Vs(:,F3s(i3,:)) - Vs(:,F3s(i1,:));
       % normalize the vectors
       pp_length = sqrt(sum(pp.^2,1));
       qq_length = sqrt(sum(qq.^2,1));
       Ipp_zero = pp_length<epsilon;
       pp_length(Ipp_zero) = 1;
       Iqq_zero = qq_length<epsilon;
       qq_length(Iqq_zero) = 1;
       pp_nor = pp ./ repmat( pp_length, [3 1] );
       qq_nor = qq ./ repmat( qq_length, [3 1] );
       % compute angles and clamped cotans
       cos_ang = sum(pp_nor.*qq_nor,1);
       cos_ang = clamp(cos_ang,-1,1);
       ang = acos( cos_ang );
       eval(['ctan_' num2str(i) ' = cot(ang) / 2.0;']);
       % clamp the cotan value so as to avoid numerical unstability
       eval(['ctan_' num2str(i) ' = clamp(ctan_' num2str(i) ', 0.001,1000);']);
       % store the right indices for the computation of laplacian matrix
       eval(['ii_lap_' num2str(i) ' = F3s(i2,:);']);
       eval(['jj_lap_' num2str(i) ' = F3s(i3,:);']);

       % accumulate the angles on vertices
       for j=1:m
           indextemp = F3s(i1,j);
           sum_angles(indextemp,1) = sum_angles(indextemp,1) + ang(1,j);
       end

       % compute the facet areas
       if i==1
           rr = crossinline(pp',-qq');
           rr = rr';
           area = sqrt( sum(rr.^2,1) ) / 2.0;
       end
    end
    
    % Laplacian matrix (stiffness matrix)
    ii_lap = [ii_lap_1 jj_lap_1 ii_lap_2 jj_lap_2 ii_lap_3 jj_lap_3];
    jj_lap = [jj_lap_1 ii_lap_1 jj_lap_2 ii_lap_2 jj_lap_3 ii_lap_3];
    ss_lap = [ctan_1 ctan_1 ctan_2 ctan_2 ctan_3 ctan_3];
    laplacian = sparse(ii_lap,jj_lap,ss_lap,n,n);
    diag_laplacian = full( sum(laplacian,1) );
    Diag_laplacian = spdiags(diag_laplacian(:),0,n,n);
    laplacian = Diag_laplacian - laplacian;
    
    % lumped mass matrix
    ii_mass = [F3s(1,:) F3s(2,:) F3s(3,:)];
    jj_mass = [F3s(1,:) F3s(2,:) F3s(3,:)];
    area = area / 3.0;
    ss_mass = [area area area];
    mass = sparse(ii_mass,jj_mass,ss_mass,n,n);

    % facet-edge adjacency matrix
    ii_adja = [F3s(1,:) F3s(2,:) F3s(3,:)];
    jj_adja = [F3s(2,:) F3s(3,:) F3s(1,:)];
    ss_adja = [1:m 1:m 1:m];
    adja = sparse(ii_adja,jj_adja,ss_adja,n,n);

    % add missing points
    I_adja = find( adja'~=0 );
    I_adja = I_adja( adja(I_adja)==0 ); 
    adja(I_adja) = -1;

    % find the boundary
    [I,J,V] = find(adja);
    I = I(V==-1);
    J = J(V==-1);
    flag_boundary = false(n,1);
    flag_boundary(I,1) = true;
    flag_boundary(J,1) = true;
    
    % set different constant values (for Gaussian curvature computation) for
    % boundary vertices and non-boundary vertices
    constants = repmat(2*pi,n,1);
    I_boundary = find(flag_boundary==true);
    constants(I_boundary) = repmat(pi,size(I_boundary));
    
    % Gaussian curvature
    cgauss = constants - sum_angles;
end

function [cGaussRough,cGaussRmean] = compute_roughness(cGauss,Lapl,Mass,SurfA)
    cGaussRough = abs(cGauss)' * Lapl;
    cGaussRough = cGaussRough' ./ spdiags(Lapl,0);
    cGaussRough = abs(cGaussRough);
    % average roughness
    cGaussRmean = cGaussRough' * spdiags(Mass,0);
    cGaussRmean = cGaussRmean / SurfA;
end

function z = crossinline(x,y)
    z = x; % x and y are (m,3) dimensional
    z(:,1) = x(:,2).*y(:,3) - x(:,3).*y(:,2);
    z(:,2) = x(:,3).*y(:,1) - x(:,1).*y(:,3);
    z(:,3) = x(:,1).*y(:,2) - x(:,2).*y(:,1);
end

function y = clamp(x,a,b)
    if nargin<2, a = 0; end
    if nargin<3, b = 1; end
    y = max(x,a);
    y = min(y,b);
end
