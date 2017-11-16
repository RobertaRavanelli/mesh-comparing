function [curvatureMetrics, curvatureData] = compareCurvature(scn1,scn2)
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
    medianEdgeLength = @(scan) median(sqrt(sum((scan.Vs(scan.F3s(:,2),:)-scan.Vs(scan.F3s(:,1),:)).^2,2)));
    dot_prod = @(a,b) a(:,1).*b(:,1)+a(:,2).*b(:,2)+a(:,3).*b(:,3);
    divideMinByMax = @(x) min(x,[],2)./max(x,[],2);
    
    tic();
    fprintf('Comparing normals, curvature, and roughness: ');
    nVs = size(scn1.Vs,1);
    edgeMedianLen = (medianEdgeLength(scn1)+ medianEdgeLength(scn2)) / 2;
    %linking the first scan to the closest 3 points in the second
    [compIDX,compDs] = knnsearch(scn2.Vs,scn1.Vs,'K',3,'NSMethod','kdtree');
    absM1s = abs(scn1.curvMean); absM2s = abs(scn2.curvMean);
    absG1s = sqrt(abs(scn1.curvGauss)); absG2s = sqrt(abs(scn2.curvGauss));
    [roughD,cGaussRough1,cGaussRough2] = Wang_perceptual_distance(scn1,scn2);
    
    i2s = compIDX(:,1);
    bestNs(:,1) = dot_prod(scn1.Vns, scn2.Vns(i2s,:));
    bestMs(:,1) = divideMinByMax( [absM1s, absM2s(i2s)] );
    bestGs(:,1) = divideMinByMax( [absG1s, absG2s(i2s)] );
    cGaussRough1(cGaussRough1 == 0) = 1; cGaussRough2(cGaussRough2 == 0) = 1;
    bestRs(:,1) = cGaussRough1./cGaussRough2(i2s); %Saved as a ratio first / last, 
    %so if the first scan is rougher than the second, the ratio generally is > 1
    
    pt2triNorms = zeros(nVs,3);
    tmpM2s = zeros(nVs,1); tmpG2s = zeros(nVs,1); tmpR2s = zeros(nVs,1);
    notifyStep = floor(nVs/10);
    for i=1:nVs
        triwts = compDs(i,1)./compDs(i,:);
        triwts = triwts./repmat(sum(triwts),1,3);
        pt2triNorms(i,1:3) = sum(scn2.Vns(compIDX(i,:),:).*repmat(triwts',1,3));
        tmpM2s(i) = sum(absM2s(compIDX(i,:))'.*triwts);
        tmpG2s(i) = sum(absG2s(compIDX(i,:))'.*triwts);
        tmpR2s(i) = sum(cGaussRough2(compIDX(i,:))'.*triwts);
        if mod(i,notifyStep)==0, fprintf('*'); end
    end
    bestNs(:,2) = dot_prod(scn1.Vns, pt2triNorms);
    bestMs(:,2) = divideMinByMax( [absM1s, tmpM2s] );
    bestGs(:,2) = divideMinByMax( [absG1s, tmpG2s] );
    bestRs(:,2) = cGaussRough1./tmpR2s;
    flipDir = bestRs > 1;
    bestRs(flipDir) = 1./(bestRs(flipDir));
    [bestRs,bestRIDX] = max(bestRs,[],2);
    flipBack = (bestRIDX==1 & flipDir(:,1)) | (bestRIDX==2 & flipDir(:,2));
    bestRs(~flipBack) = bestRs(~flipBack)-1;
    bestRs(flipBack) = 1-bestRs(flipBack);
    
    bestComparisons = [max(bestNs,[],2), max(bestMs,[],2), max(bestGs,[],2), bestRs];
    fprintf([repmat(char(8),1,10),'Completed. ']);
    toc();
    
    angles = acosd(bestComparisons(:,1));
    scoreN = (30-angles)/30;
    scoreN(scoreN<0) = 0;
    scores = scoreN.*mean([bestComparisons(:,2),bestComparisons(:,3)],2);
    
    curvatureData.medianEdge = edgeMedianLen;
    curvatureData.scores = scores;
    curvatureData.bestNormals = bestComparisons(:,1);
    curvatureData.bestCurvMean = bestComparisons(:,2);
    curvatureData.bestCurvGauss = bestComparisons(:,3);
    curvatureData.bestRoughness = bestComparisons(:,4);
    
    curvatureMetrics.roughDistance = roughD;
    curvatureMetrics.normals = getMetrics(bestComparisons(:,1));
    curvatureMetrics.curvMean = getMetrics(bestComparisons(:,2));
    curvatureMetrics.curvGauss = getMetrics(bestComparisons(:,3));
    curvatureMetrics.scores = getMetrics(scores);
end

function metrics = getMetrics(bestVals)
    sdDist = 100-[68.27,95.45,99.73];
    metrics.median = median(bestVals);
    metrics.sig1 = prctile(bestVals,sdDist(1));
    metrics.sig2 = prctile(bestVals,sdDist(2));
    metrics.sig3 = prctile(bestVals,sdDist(3));
end
