function [metrics] = compare3Dscans(refScan, floatScan, Iterations, save3Dviz, stageProcessing)
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
    saveOBJ = true; %Always save the results
    getCurvature = false; testScale = true; maxits=1;
    switch stageProcessing
      case 2
        suffix = 'rescale2';
        maxits=2;
      case 3
        suffix = 'final3';
        getCurvature = true;
      otherwise
        suffix = 'stage1';
        testScale = false;
    end
    
    scaleFac = 0.0; pass=0; sF = 0.0; 
    while pass < maxits
        if testScale, [sF,floatScan.Vs] = resizeScan(refScan,floatScan); end
        scaleFac=scaleFac+sF;
        
        tic();
        Tolerance = max(refScan.maxXYZ)/200000; %Set below the lowest conceivable precision of the model
        PyramidLevels = 6; Visualize = 0;
        %Reorient the "floating" scan to the reference with Tolga's ICP
        FinalPose = Tolga_icp_mod_point_plane_pyr(floatScan.Vs, floatScan.Vns, refScan.Vs, refScan.Vns, Tolerance, Iterations, 3, 1, PyramidLevels, 0, Visualize);
        floatScan.Vs = Tolga_movepoints(FinalPose, floatScan.Vs);
        floatScan.maxXYZ = [max(floatScan.Vs(:,1))-min(floatScan.Vs(:,1)), ...
            max(floatScan.Vs(:,2))-min(floatScan.Vs(:,2)), max(floatScan.Vs(:,3))-min(floatScan.Vs(:,3))];
        floatScan.Vns = vertexNormal(triangulation(floatScan.F3s,floatScan.Vs));
        toc();
        
        tic();
        %Estimate the point-plane distances from the reference (discarding the
        %equivalent points from the floating scan to the reference returned by the function)
        floatScan.signedDistances = estimatePointPlane(refScan.F3s, refScan.Vs, floatScan.Vs);
        dists = abs(floatScan.signedDistances);
        toc();
        
        Vquality = zeros(size(floatScan.Vs,1),1);
        maxdist = max(dists);
        fac = maxdist - median(dists); %Any distances below the median will be scored 1 for quality
        for i=1:size(floatScan.Vs,1)
            Vquality(i) = (maxdist - dists(i))/fac;
            if Vquality(i) > 1, Vquality(i) = 1; end
        end
        floatScan.Qs = Vquality;
        
        %Calculate some metrics describing the distribution
        errMedian = 1000*prctile(dists,50);
        err1sigma = 1000*prctile(dists,68.27);
        err2sigma = 1000*prctile(dists,95.45);
        err3sigma = 1000*prctile(dists,99.73);
        errMax = 1000*max(dists);
        metrics = [errMedian, err1sigma, err2sigma, err3sigma, errMax, scaleFac];
        pass = pass+1;
    end
	
    close(gcf);
    fig = figure;
    if save3Dviz %Draw the two models together, coloring the floating scan's faces by the
        hold on; %distances, and plotting the points of the reference scan in blue
        trisurf(floatScan.F3s, floatScan.Vs(:,1), floatScan.Vs(:,2), floatScan.Vs(:,3), ...
            'FaceVertexCData', dists, 'LineStyle','none', 'FaceColor','interp', 'FaceAlpha',0.8);
        caxis([0, max(dists)]);
        colormap(summer);
        scatter3(refScan.Vs(:,1), refScan.Vs(:,2), refScan.Vs(:,3), 4, 'b.');
        colorbar;
        title([floatScan.basename,': surface errors of ',floatScan.shortname,' from reference points of ',refScan.shortname]);
        axis equal;
        view(-45,45);
        drawnow;
        filename = [floatScan.basename,'-',floatScan.shortname,'_fromRef',refScan.shortname,'-3Dviz.fig'];
        savefig(fig, filename);
        hold off;
    end
    %Save a histogram for the distances of the points
    clf(fig); hold on;
    ax = fig.CurrentAxes;
    h = histogram(ax,dists,200);
    h.EdgeColor='none';
    axis(ax,'tight');
    title(ax,[floatScan.basename,': Error of ',floatScan.shortname,' from reference ',refScan.shortname]);
    drawnow;
    filename = [floatScan.basename,'-',floatScan.shortname,'_fromRef',refScan.shortname,'-histogram.fig'];
    savefig(fig, filename);
    hold off;
    
    if getCurvature
        [floatScan.Vns,floatScan.smoothedNormals,floatScan.curvMean,floatScan.curvGauss] = getCurvatureBezierPatch(floatScan);
    end
    
    if saveOBJ, save([floatScan.filename(1:end-4),suffix,'.mat'],'floatScan'); end
end

function [scaleFac,newVs] = resizeScan(refIn,floatIn)
    dot_prod = @(a,b) a(:,1).*b(:,1)+a(:,2).*b(:,2)+a(:,3).*b(:,3);
    
    %The weighted centroid point is averaged for the floating and the reference scan. It is assumed
    %that the two have already be registered to one another by ICP. Weighting is by face areas.
    awPtsRef = repmat(getFaceAreaWeights(refIn.Vs,refIn.F3s),1,3);
    awPtsFloat = repmat(getFaceAreaWeights(floatIn.Vs,floatIn.F3s),1,3);
    weightedCentroid = mean([sum(refIn.Vs.*awPtsRef,1)./sum(awPtsRef(:,1)); ...
        sum(floatIn.Vs.*awPtsFloat,1)./sum(awPtsFloat(:,1))],1);
    %Precompute the two scan's points' distances from the centroid
    centeredFloatPts = floatIn.Vs-repmat(weightedCentroid,size(floatIn.Vs,1),1);
    floatDistctr = sqrt(sum(centeredFloatPts.^2,2));
    relativeDir = bsxfun(@rdivide,centeredFloatPts,floatDistctr); %Normalized vectors
    centeredRefPts = refIn.Vs-repmat(weightedCentroid,size(refIn.Vs,1),1);

    %Project floating points onto the reference along the vector originating in the centroid
    fltCentered.Vs = centeredFloatPts; fltCentered.Vns = relativeDir; fltCentered.F3s = floatIn.F3s;
    refCentered.Vs = centeredRefPts; refCentered.Qs = refIn.Qs; refCentered.F3s = refIn.F3s;
    projectedRefPts = createAveragedModel(fltCentered, refCentered, false);
    relativeDists = sqrt(sum(projectedRefPts.^2,2)) - floatDistctr;

    Vweights = abs(dot_prod(floatIn.Vns,relativeDir));
    kdOBJ = KDTreeSearcher(refIn.Vs);
    [~,closestRefPtDist] = knnsearch(kdOBJ,floatIn.Vs);
    discardedIDX = abs(relativeDists) > 2*closestRefPtDist;
    %Since the points are being projected on the reference scan, we expect the relative distance
    %to be no more than about twice the distance to the closest point on the reference scan.
    %The points rejected as outliers may be projected in the wrong area of the reference scan.
    relativeDists(discardedIDX) = []; Vweights(discardedIDX) = []; floatDistctr(discardedIDX) = [];

    scaleFac = sum(relativeDists.*Vweights)/sum(floatDistctr.*Vweights);
    newVs = centeredFloatPts*(1+scaleFac) + repmat(weightedCentroid,size(floatIn.Vs,1),1);
end

function Vweights = getFaceAreaWeights(Vs,F3s)
    Vweights = zeros(size(Vs,1),1);
    for i=1:size(F3s,1)
        centroid = mean(Vs(F3s(i,:),:),1);
        for j=1:3, Vweights(F3s(i,j)) = Vweights(F3s(i,j)) + sqrt(sum((centroid-Vs(F3s(i,j),:)).^2)); end
    end
    Vweights = Vweights/max(Vweights);
end
