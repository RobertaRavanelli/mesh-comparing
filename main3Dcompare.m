%%
% https://github.com/psapirstein/mesh-comparing
% Distributed under an Apache License 2.0
%
% The code supplements the article:
% "A high-precision photogrammetric recording system for small artifacts"
% Philip Sapirstein, UNL
% Journal of Cultural Heritage 2017
% https://doi.org/10.1016/j.culher.2017.10.011
%
% The operations are described in the appendices, available as a free
% supplement to the online version of the paper.
%
% This is a collection of Matlab scripts for manipulating 3D mesh data
% saved in OBJ format. Its primary function is to compare multiple scans of
% the same object, and to assess the discrepancies among the meshes. The
% repeatability of the scan data for various parameters is reported, and
% reoriented scans are saved as OBJ and MAT files. The code is designed for
% working with multiple scans, but it might easily be adapted to the
% simpler problem of quantifying the error of a single input mesh relative
% to a reference scan created at a higher accuracy and precision.
%
% The subroutines include several stages of analysis detailed below.
% First, scans are loaded, the geometry repaired, and the scan reoriented
% rigidly relative to a common reference scan selected by the user.
% Second, the scale of each scan is adjusted relative to a reference model,
% which quantifies any small discrepancies in overall scale in the scans.
% Third, an idealized model is created that interpolates points from the
% four scans of highest quality (meaning that at least four models of the
% same object are required for the code to function), including a scoring
% system to reduce the influence of outliers.
% Fourth, the individual scans are compared to this reference, saving final
% estimates of repeatibility. An optional fifth stage assesses curvature.
%%
clear variables;
codeStage = 1; %Manual switch to select among the stages of processing
%Modify the paths for the system running the script
basepath = ''; %Set this to the path where this script resides
addpath([basepath,'PScan3Dcompare'],'-end'); %Algorithms are kept in a subfolder
%Set up the parallel computing
numWorkers = 6; %Suited to an 8-core machine
paralleljob = gcp('nocreate');
if isempty(paralleljob)
    paralleljob = parpool('local8',numWorkers);
elseif paralleljob.NumWorkers ~= numWorkers
    delete(paralleljob);
    paralleljob = parpool('local8',numWorkers);
end

switch codeStage
  case 1
%%  Stage 1: Preliminary inspection of scans, cleaning meshes, approximate positioning near reference
    %The directories where the scans reside, followed by the base name of the scans for comparison
    %The scans must be in OBJ format and start with the "filename_base" followed by individual scan names.
    scan_dirs = {''; ''; ''}; %Set to the appropriate path
    filename_bases = {'Deer'; 'Skull'; 'Pot'}; %For the scans discussed in the article
    %The reference scan will be one chosen from the set, presumably centered near the origin and oriented as desired, to which the
    %others will be approximately fitted. The best scans will be averaged in the following stage, producing a better reference model.
    reference_dir = ''; %Set to the appropriate path
    Iter = 5; %Iterations for the ICP stage; low values are sufficient if the scans are already at approximately the correct orientation
    save3Dviz = true; %If false, only the histogram of the comparison will be saved
    stageProc = 1; %Switch for the different operations invoked at each stage.

    fid = fopen('_firststageComparisonOBJ.csv','at+');
    fprintf(fid, '\nProcess initiated %s\n', datestr(now));
    fprintf(fid, 'Scan, Reference, Float, Median (mm), 1-sigma error, 2-sigma, 3-sigma, Max\n');
    reffiles = dir(fullfile(reference_dir,'*.obj'));
    ref_names = {reffiles.name}';
    if isempty(ref_names), error('Invalid reference directory or no .OBJ files present with the given name.'); end

    for objcounter = 1:length(scan_dirs)
        objfilelist = dir(fullfile(scan_dirs{objcounter}, strcat(filename_bases{objcounter},'*.obj')));
        obj_names = {objfilelist.name}';
        if isempty(obj_names)
            fprintf('%s: Invalid scan directory or no .OBJ files present with the given name.\n', filename_bases{objcounter});
            continue;
        end

        refIDX = find(cellfun(@(x) ~isempty(strfind(x,filename_bases{objcounter})),ref_names));
        refScan = formatOBJscan(ref_names{refIDX},reference_dir,filename_bases{objcounter});
        for i = 1:length(obj_names)
            fprintf('\n');
            tic;
            floatScan = formatOBJscan(obj_names{i},scan_dirs{objcounter},filename_bases{objcounter});
            [repaired,newVs,newF3s,deletedVs] = repairMeshTriangulation(floatScan);
            if repaired>0
                floatScan.Vs = newVs; floatScan.F3s = newF3s;
                floatScan = rmfield(floatScan,'LineByLine');
                if ~isempty(deletedVs)
                    floatScan.Vcs(deletedVs,:) = []; floatScan.Vns(deletedVs,:) = []; floatScan.Qs(deletedVs) = [];
                end
            end

            clear newVs newF3s deletedVs;
            toc;

            if repaired<0
                fprintf('Error: scan must be manually repaired before proceeding.');
                fprintf(fid, '%s must be manually repaired and reprocessed.',floatScan.filename);
            else
                metrics = compare3Dscans(refScan, floatScan, Iter, save3Dviz, stageProc);
                fprintf('Errors from reference (microns):\n  Median: %.1f; 1-sigma: %.1f; 2-sigma: %.1f; 3-sigma: %.1f; max: %.1f\n',metrics(1)*1000,metrics(2)*1000,metrics(3)*1000,metrics(4)*1000,metrics(5)*1000);
                fprintf(fid, '%s,%s,%s,%.6f,%.6f,%.6f,%.6f,%.6f\n', refScan.basename, refScan.shortname, floatScan.shortname, ...
                    metrics(1), metrics(2), metrics(3), metrics(4), metrics(5));
            end
        end
    end
    fclose(fid);
  case 2
%%  Stage 2: examine scale errors and refine the scales of the scans
    scan_dir = 'reorientedScansNew'; %At this stage, all the MAT formatted files to be processed should be kept here
    filename_bases = {'Deer'; 'Skull'; 'Pot'}; %With this prefix
    reference_dir = ''; %Set to the appropriate path
    Iter = 2; save3Dviz = false; stageProc = 2;

    fid = fopen('_secondRescaling.csv','at+');
    fprintf(fid, '\nProcess initiated %s\n', datestr(now));
    fprintf(fid, 'Scan, Reference, Float, Scaling (percentage), Median (mm), 1-sigma error, 2-sigma, 3-sigma, Max\n');
    reffiles = dir(fullfile(reference_dir,'*.obj'));
    ref_names = {reffiles.name}';
    if isempty(ref_names), error('Invalid reference directory or no .OBJ files present with the given name.'); end

    for objcounter=1:length(filename_bases)
        objfilelist = dir(fullfile(pwd,scan_dir, strcat(filename_bases{objcounter},'*.mat')));
        obj_names = {objfilelist.name}';
        if isempty(obj_names)
            fprintf('%s: Invalid scan directory or no .OBJ files present with the given name.\n', filename_bases{objcounter});
            continue;
        end

        refIDX = find(cellfun(@(x) ~isempty(strfind(x,filename_bases{objcounter})),ref_names));
        refScan = formatOBJscan(ref_names{refIDX},reference_dir,filename_bases{objcounter});
        for scans=1:length(obj_names)
            file = load(fullfile(scan_dir,obj_names{scans}));
            floatScan = file.floatScan; clear file;

            metrics = compare3Dscans(refScan, floatScan, Iter, save3Dviz, stageProc);
            fprintf('Errors from reference (microns):\n  Rescaling by: %.3f%%; Median: %.1f; 1-sigma: %.1f; 2-sigma: %.1f; 3-sigma: %.1f; max: %.1f\n',100*metrics(6),metrics(1)*1000,metrics(2)*1000,metrics(3)*1000,metrics(4)*1000,metrics(5)*1000);
            fprintf(fid, '%s,%s,%s,%.10f,%.6f,%.6f,%.6f,%.6f,%.6f\n', refScan.basename, refScan.shortname, floatScan.shortname, ...
                    metrics(6)*100, metrics(1), metrics(2), metrics(3), metrics(4), metrics(5));
        end
    end
    fclose(fid);
  case 3
%%  Stage 3: create averaged scan data
    scan_dir = 'reorientedRescaledScans'; %At this stage, all the MAT files to be processed should be kept here
    filename_bases = {'Deer'; 'Skull'; 'Pot'}; %With this prefix
    best_pairs = [2,3; 6,7]; %The code assumes two pairs selected by the user for each scale
    averagePts = true; %Very few iterations are needed because the scans have already been fitted by ICP numerous times

    for objcounter=1:length(filename_bases)
        objfilelist = dir(fullfile(pwd,scan_dir, strcat(filename_bases{objcounter},'*.mat')));
        obj_names = {objfilelist.name}';
        err2sigma = ones(2,1); %Only set up to work with two pairs
        floatScan = cell(2,1); %Need to save both averaged scans for second averaging
        for pairs = 1:2 %Arbitrarily fitting the second to the first in each pair
            file = load(fullfile(scan_dir,obj_names{best_pairs(pairs,1)}));
            refScan = file.floatScan;
            file = load(fullfile(scan_dir,obj_names{best_pairs(pairs,2)}));
            floatScan{pairs} = file.floatScan; clear file;
            
            fprintf('\nAveraging %s for the %s/%s pairs\n',filename_bases{objcounter}, ...
                refScan.shortname, floatScan{pairs}.shortname);
            [floatScan{pairs}.Vs, floatScan{pairs}.Qs, err2sigma(pairs)] = createAveragedModel(floatScan{pairs}, refScan, averagePts);
            %The floating scan is updated with the averaged & partly
            %smoothed points based on its face triangulation. Vertex color
            %is deleted (it will be replaced by quality) & normals regenerated.
            floatScan{pairs}.Vcs = [];
            floatScan{pairs}.Vns = vertexNormal(triangulation(floatScan{pairs}.F3s,floatScan{pairs}.Vs));
        end
        clear refScan;
        %Adjust the weights based on the pairs' 2-sigma errors
        errRatio = err2sigma(2)/err2sigma(1);
        if errRatio > 1 %The second scan pair is weaker (ie, higher error)
            preferred = 1; reference = 2;
        else
            preferred = 2; reference = 1;
            errRatio = 1/errRatio; %Invert the error ratio so that it is greater than one
        end
        floatScan{reference}.Qs = floatScan{reference}.Qs/errRatio; %Lower the quality of the weaker scan

        fprintf('\nAveraging %s (preferred) with %s\n',floatScan{preferred}.shortname, floatScan{reference}.shortname);
        [floatScan{preferred}.Vs, floatScan{preferred}.Qs, ~] = createAveragedModel(floatScan{preferred}, floatScan{reference}, averagePts);
        fprintf('Saving the idealized reference model.\n');
        refScan = floatScan{preferred};
        
        save([refScan.basename,refScan.shortname,'master.mat'],'refScan');
        [repaired,newVs,newF3s,deletedVs] = repairMeshTriangulation(refScan);
        if repaired>0
            refScan.Vs = newVs; refScan.F3s = newF3s;
            if ~isempty(deletedVs)
                refScan.Vns(deletedVs,:) = []; refScan.Qs(deletedVs) = [];
            end
        end
        refScan = rmfield(refScan,'signedDistances');
        save([refScan.basename,refScan.shortname,'master.mat'],'refScan');
    end
  case 4
%%  Stage 4: Reorient all scans to the new reference scans, save final metrics
    scan_dir = 'reorientedRescaledScans'; %At this stage, all the MAT files to be processed should be kept here
    filename_bases = {'Deer'; 'Skull'; 'Pot'}; %With this prefix
    reference_dir = 'masterScansNew'; %Directory where the reference scans are stored
    Iter = 2; save3Dviz = true; stageProc = 3;
    
    fid = fopen('_finalComparisonOBJ.csv','at+');
    fprintf(fid, 'Scan, Reference, Float, Scaling (percentage), Median (mm), 1-sigma error, 2-sigma, 3-sigma, Max\n');
    
    reffiles = dir(fullfile(reference_dir,'*.mat'));
    ref_names = {reffiles.name}';
    if isempty(ref_names), error('Invalid reference directory or no .mat files present with the given name.'); end
    
    for objcounter=1:length(filename_bases)
        objfilelist = dir(fullfile(pwd,scan_dir, strcat(filename_bases{objcounter},'*.mat')));
        obj_names = {objfilelist.name}';
        if isempty(obj_names)
            fprintf('%s: Invalid scan directory or no .OBJ files present with the given name.\n', filename_bases{objcounter});
            continue;
        end
        
        refIDX = find(cellfun(@(x) ~isempty(strfind(x,filename_bases{objcounter})),ref_names));
        file = load(fullfile(reference_dir,ref_names{refIDX}));
        refScan = file.refScan;
        if ~isfield(refScan,'curvMean')
            [refScan.Vns,refScan.smoothedNormals,refScan.curvMean,refScan.curvGauss] = getCurvatureBezierPatch(refScan);
            save(['master',refScan.basename,refScan.shortname,'curvTest.mat'],'refScan');
        end
        for scans=1:length(obj_names)
            file = load(fullfile(scan_dir,obj_names{scans}));
            floatScan = file.floatScan; clear file;
            
            metrics = compare3Dscans(refScan, floatScan, Iter, save3Dviz, stageProc);
            
            fprintf('Errors from reference (microns):\n  Rescaling by: %.3f%%; Median: %.1f; 1-sigma: %.1f; 2-sigma: %.1f; 3-sigma: %.1f; max: %.1f\n',100*metrics(6),metrics(1)*1000,metrics(2)*1000,metrics(3)*1000,metrics(4)*1000,metrics(5)*1000);
            fprintf(fid, '%s,%s,%s,%.10f,%.6f,%.6f,%.6f,%.6f,%.6f\n', refScan.basename, refScan.shortname, floatScan.shortname, ...
                    metrics(6)*100, metrics(1), metrics(2), metrics(3), metrics(4), metrics(5));
        end
        fprintf('\n');
    end
    fclose(fid);
  case 5
%%  Final switch to compare curvature after it has been estimated
    scan_dir = ''; %Set to the appropriate path
    reference_dir = ''; %Set to the appropriate path
    filename_bases = {'Deer'; 'Skull'; 'Pot'};
    pairIDX(:,1) = 1:14; pairIDX(:,2) = 0;
    pairIDX = [pairIDX; 6,7; 8,5; 6,5; 2,3; 1,4; 1,2; 6,2; 8,4; 9,13; 14,11; 12,10; 9,10; 14,13; 10,11; 12,14; 10,2; 6,10];
    
    reffiles = dir(fullfile(reference_dir,'*.mat'));
    ref_names = {reffiles.name}';
    if isempty(ref_names), error('Invalid reference directory or no .mat files present with the given name.'); end

    fid = fopen('_curveComparisons.csv','at+');
    headerStr = ['\nObject, Scan, Reference, Normal Angle (º) median, Nº-1sig, Nº-2sig, Nº-3sig,', ...
        'curveMean Ratios: median, cMD-2sig, cMD-3sig, curveGauss Ratios: median, cGD-2sig, cGD-3sig,', ...
        'meanScore: median, S-2sig, roughScore: global, mean\n'];
    fprintf(fid, headerStr);
    for objcounter=1:length(filename_bases)
        objfilelist = dir(fullfile(scan_dir,strcat(filename_bases{objcounter},'*.mat')));
        obj_names = {objfilelist.name}';
        if isempty(obj_names)
            fprintf('%s: Invalid scan directory or no .OBJ files present with the given name.\n', filename_bases{objcounter});
            continue;
        end
        
        fprintf('Analyzing model set %s:\n',filename_bases{objcounter});
        refIDX = find(cellfun(@(x) ~isempty(strfind(x,filename_bases{objcounter})),ref_names));
        file = load(fullfile(reference_dir,ref_names{refIDX}));
        refScan = file.refScan;
        for pair=1:size(pairIDX,1)
            scn = cell(0);
            for i=1:2
                if pairIDX(pair,i)==0
                    scn{i} = refScan;
                else
                    file = load(fullfile(scan_dir,obj_names{pairIDX(pair,i)}));
                    scn{i} = file.floatScan;
                end
            end
            clear file;
            compName = scn{2}.shortname;
            if pairIDX(pair,2) == 0, compName = ['ref(',compName,')']; end %#ok<AGROW>
            fprintf('\nComparing %s with %s.\n',scn{1}.shortname,compName);
            [cMtr,curvComparisonData] = compareCurvature(scn{1},scn{2});
            
            fprintf('Median score: %.2f; Normal median/2-s angles: %.1fº/%.1fº; Curv (M/Gav) median/2-s ratios: %.2f/%.2f; Rough global/mean: %.3f/%.3f\n', ...
                cMtr.scores.median, acosd(cMtr.normals.median),acosd(cMtr.normals.sig2), (cMtr.curvMean.median+cMtr.curvGauss.median)/2, ...
                (cMtr.curvMean.sig2+cMtr.curvGauss.sig2)/2, cMtr.roughDistance, mean(curvComparisonData.bestRoughness) );
            fprintf(fid, '%s,%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n', ...
                refScan.basename, scn{1}.shortname, compName, ...
                acosd(cMtr.normals.median),acosd(cMtr.normals.sig1),acosd(cMtr.normals.sig2),acosd(cMtr.normals.sig3), ...
                cMtr.curvMean.median, cMtr.curvMean.sig2, cMtr.curvMean.sig3, ...
                cMtr.curvGauss.median, cMtr.curvGauss.sig2, cMtr.curvGauss.sig3, ...
                cMtr.scores.median, cMtr.scores.sig2, cMtr.roughDistance, mean(curvComparisonData.bestRoughness) );
            
            curvComparisonData.testedScan = scn{1}.shortname;
            curvComparisonData.referenceScan = compName;
            save([refScan.basename,'ccomp',scn{1}.shortname,'-',compName,'data.mat'],'curvComparisonData');
        end
        fprintf('\n');
    end
    fclose(fid);
end

if ~isempty(paralleljob)
    delete(paralleljob);
end
clear paralleljob;
