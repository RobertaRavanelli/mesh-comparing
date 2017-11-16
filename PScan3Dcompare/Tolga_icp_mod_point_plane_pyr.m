function [Pose] = Tolga_icp_mod_point_plane_pyr(SrcPC, SrcN, DstPC, DstN, Tolerence, MaxIterations, RejectionScale, NumNeighborsCorr, NumLevels, SampleType, DrawRegistration)
%% This code has been copied from a script by Tolga Birdal:
% https://www.mathworks.com/matlabcentral/fileexchange/47152-icp-registration-using-efficient-variants-and-multi-resolution-scheme?s_tid=prof_contriblnk
%
% https://github.com/psapirstein/mesh-comparing
% It has been modified for compatibility with the mesh comparison tools
% accompanying the article:
% "A high-precision photogrammetric recording system for small artifacts"
% Philip Sapirstein, Journal of Cultural Heritage 2017
% https://doi.org/10.1016/j.culher.2017.10.011
%
% The Apache License 2.0 for the rest of the code base does not apply to
% this subroutine, whose original form is found at the above link. It has
% been included here for convenience.
%%

    % 1. A smart sampling option is provided (Gelfand et. al.). To use it set
    % SampleType=1. However note that this function is provided seperately at:
    % http://www.mathworks.com/matlabcentral/fileexchange/47138-stable-sampling-of-point-clouds-for-icp-registration
    % 2. Robust outlier filtering is incorporated
    % 3. Linear point-to-plane metric is minimized
    % 4. A correspondence filtering as in picky icp was utilized
    % 5. Multi-Resolution scheme is utilized (Coarse-to-Fine ICP)
    % Morover, a Hartley-Zissermann type of normalization is done.
    % Sample usage: FinalPose = icp_mod_point_plane_pyr(SrcPC, SrcN, DstPC, DstN, 0.05, 100, 3, 1, 8, 0, 1);
    % Parameters:
    %   SrcPC :         Nx3 array of input point cloud model
    %   SrcN :          Nx3 array of the normals of the model
    %   DstPC :         Nx3 array of the target point cloud (scene)
    %   DstN :          Nx3 array of the normals of the target
    %   Tolerence:      If there is insufficient change between iterations, ICP is
    %                   terminated at this level. This insufficiency is determined 
    %                   by Tolerence.
    %   MaxIterations:  Maximum number of iterations. Note that this parameter is
    %                   updated through the pyramid. Less iterations are
    %                   carried out in the coarser levels
    %   RejectionScale: Threshold for outlier rejection. RejectionScale*stddev
    %                   is set as the threshold of rejection. A typically used
    %                   value is 3.
    %   NumNeighborsCorr: Not used for now
    %   NumLevels:      Maximum levels of the pyramid (depth of the search)
    %   SampleType:     0 for uniform sampling, 1 for stable sampling as in 
    %                   Gelfand et. al.
    %   DrawRegistration: 1 for visualization, 0 for no visualizations
    %   Author: Tolga Birdal

    % Assign variables. Leave as default if not provided
    TolP = 0.0001;
    Iterations = 100;
    useRobustRejection = 0;
    outlierScale = 3;
    %numNeighbors = 1; % Allow for N multiple assignments
    visualize = 1;
    nK = 1;
    numPyr = 8;
    samplePCType = 0;
    % Assign parameters
    if (nargin > 4), TolP = Tolerence;                  end
    if (nargin > 5), Iterations = MaxIterations;        end
    if (nargin > 6)
        useRobustRejection = (RejectionScale>0);
        outlierScale = RejectionScale;
    end
    if (nargin > 7)
        numNeighbors = NumNeighborsCorr;

        if (numNeighbors>0)
            nK = numNeighbors;
        end
    end
    if (nargin > 8), numPyr = NumLevels;      end
    if (nargin > 9), samplePCType = SampleType;      end
    if (nargin > 10), visualize = DrawRegistration;      end

    n = length(SrcPC);

    % Hartley-Zissermann Scaling:
    meanSrc = mean(SrcPC);
    meanDst = mean(DstPC);
    meanAvg = (meanSrc+meanDst)*0.5;
    SrcPC = bsxfun(@minus, SrcPC, meanAvg);
    DstPC = bsxfun(@minus, DstPC, meanAvg);

    % compute average dist from origin
    avgDist = (sum(sqrt(sum(SrcPC.^2, 2))) + sum(sqrt(sum(DstPC.^2, 2))))*0.5;

    % scale to unit sphere
    scale = n / avgDist;
    SrcPC = SrcPC.*scale;
    DstPC = DstPC.*scale;

    SrcPCOrig = SrcPC;
    DstPCOrig = DstPC;
    SrcNOrig = SrcN;
    DstNOrig = DstN;

    Pose = [eye(3) zeros(3,1)] ;

    if (visualize)
        figure('units','normalized','outerposition',[0.4 0.3 0.6 0.7]);
    end
    IterationsInit = Iterations;
    TolPInit = TolP;

    % Construct the KD-tree out of destination points (scene)
    DstPC = DstPCOrig;
    DstN = DstNOrig;
    kdtreeobj = KDTreeSearcher(DstPC,'distance','euclidean');

    fprintf('Performing ICP, Level: ');
    
    % Start the registration
    for level=numPyr-1:-1:0
        fprintf('%d ',level);
        
        % Obtain the parameters in this pyramid level
        numSamples = uint32(n./(2^level));
        %TolP = TolPInit.*(2^(2*level));
        TolP = TolPInit.*((level+1)^2);
        Iterations = uint32(IterationsInit./(2^(level)));

        % Obtain the sampled point clouds for this level
        %Pose=PoseInit;
        SrcPC = Tolga_movepoints(Pose, SrcPCOrig);
        rotPose = Pose;
        rotPose(:,4)=[0,0,0]';
        SrcN = Tolga_movepoints(rotPose, SrcNOrig);   

        if (samplePCType==1)
            [SrcPC, ~] = sample_pc_stable(SrcPC, SrcN, numSamples);
        else
            [SrcPC, ~] = sample_pc_uniform(SrcPC, SrcN, numSamples);
        end

        % Distance error in last iteration
        fval_old=inf;

        % Change in distance error between two iterations
        fval_perc=0;

        % Array which contains the transformed points
        Src_Moved=SrcPC;

        %minPose=[]; % we will keep track of the pose with minimum residual
        fval_min=inf;

        i=0;

        if (visualize)
            plot3(Src_Moved(:,1), Src_Moved(:,2), Src_Moved(:,3),'r.');
            hold on, plot3(DstPC(:,1), DstPC(:,2), DstPC(:,3),'b.');

            hold off;
            axis equal;
            %pause(1);
            view(-124,60);
            drawnow;
        end

        PoseX = [1,0,0,0;0,1,0,0;0,0,1,0];

        % Start main loop of ICP
        while( (~(fval_perc<(1+TolP) && fval_perc>(1-TolP))) && i<Iterations)

            % Calculate closest point for all points
            [j,d]=knnsearch(kdtreeobj,Src_Moved,'k',nK);

            % Implement Picky ICP or BC-ICP
            newI = (1:length(j))';
            newJ = j;

            % Step 1 of ICP : Robustly reject outliers
            if (useRobustRejection>0)
                % I don't like this as it is very dependent on the dataset
                dMin = d(:,1);
                med = median(dMin);
                sigma = 1.4826 * mad(dMin,1); % Robust estimate of stddev
                threshold = (outlierScale*sigma+med);
                acceptInd = (dMin<=threshold);
                acceptInd = find(acceptInd(:));
                newJ = j(acceptInd, :);
                newI = newI(acceptInd);
            end

            % Step 2 of Picky ICP:
            % Among the resulting corresponding pairs, if more than one scene point p_i
            % is assigned to the same model point m_j, then select p_i that corresponds
            % to the minimum distance.

            duplicates=get_duplicates(newJ);

            for di=1:length(duplicates)
                dup = duplicates(di);

                if (isempty(dup))
                    continue;
                end

                % Such search could be done much faster i.e. when implemented in C
                % Using, say, a hashtable
                indJ = (find(newJ==dup));

                dists = d(indJ);
                if (isempty(indJ))
                    continue;
                end

                [~, indD] = min(dists);

                tempI = newI( indJ(indD) );
                newJ(indJ) = NaN;
                newI(indJ) = NaN;
                newJ( indJ(1) ) = dup;
                newI( indJ(1) ) = tempI;
            end

            mapSrc = newI;
            mapDst = newJ;

            mapSrc = mapSrc(~isnan(mapSrc));
            mapDst = mapDst(~isnan(mapDst));

            Src_Match=SrcPC(mapSrc,:);
            Dst_Match=DstPC(mapDst,:);
            DstN_Match=DstN(mapDst,:);

            x = minimize_point_to_plane(Src_Match, Dst_Match, DstN_Match);

            % Make the transformation matrix
            PoseX = get_transform_mat(x);

            % Transform the Points
            Src_Moved = Tolga_movepoints(PoseX, SrcPC);

            fval = sum( (Src_Moved - SrcPC).^2, 2);
            fval = sqrt(sum(fval)./length(fval));

            % Calculate change in error between iterations
            fval_perc=fval/fval_old;

            % Store error value
            fval_old=fval;

            if (fval < fval_min)
                fval_min = fval;
                %minPose = Pose;
            end
            
            % release some memory, just in case
            clear Src_Match;
            clear Dst_Match;
            clear DstN_Match;

            i=i+1;
        end

        Pose = [PoseX; 0 0 0 1]*[Pose; 0 0 0 1];
        Pose = Pose(1:3, 1:4);
    end
    
    Pose(1:3, 4) = Pose(1:3, 4)./scale + meanAvg' - Pose(1:3, 1:3)*meanAvg';

    if (visualize)
        plot3(SrcPC(:,1), SrcPC(:,2), SrcPC(:,3),'r.');
        hold on, plot3(DstPC(:,1), DstPC(:,2), DstPC(:,3),'b.');

        hold off;
        axis equal;
        %pause(1);
        %view(-124,60);
        drawnow;
    end
end

%%
function [SrcSample, SrcSampleNormals] = sample_pc_uniform(Src, NormalsSrc, numPoints)
    % by Tolga Birdal
    % Uniform sampling of point clouds

    n = length(Src);

    % Really sample them.
    sampledIndices=1:fix(n/numPoints):n;

    SrcSample = Src(sampledIndices, :);
    SrcSampleNormals = NormalsSrc(sampledIndices, :);
end

%%
function [x] = minimize_point_to_plane(Src, Dst, Normals)
    % Minimize the point to plane metric according to Kok Lim Low : Linear
    % Least Squares Optimization for Point-to-Plane ICP Surface Registration
    % Also see Efficient Variants of the ICP Algorithm by Szymon Rusinkiewicz
    % Author: Tolga Birdal

    b = dot(Dst-Src, Normals, 2);
    A1 = cross(Src, Normals);
    A2 = Normals;
    A=[A1 A2];
    x = (A\b)';
end

%%
function duplicates = get_duplicates(X)
    % Returns the duplicates of a vector
    % 
    % Author: Tolga Birdal
    uniqueX = unique(X);
    countOfX = hist(X,uniqueX);
    index = (countOfX~=1 & countOfX~=0);
    duplicates = uniqueX(index);
end

%%
function M = get_transform_mat(par)
    % Given the parameters as [\theta_x, \theta_y, \theta_z, t_x, t_y, t_z]
    % Returns the pose matrix [R | t]
    % Author: Tolga Birdal

    r=par(1:3);
    t=par(4:6);
    Rx=[1 0 0 ;
        0 cos(r(1)) -sin(r(1)) ;
        0 sin(r(1)) cos(r(1)) ];

    Ry=[cos(r(2)) 0 sin(r(2)) ;
        0 1 0;
        -sin(r(2)) 0 cos(r(2))];

    Rz=[cos(r(3)) -sin(r(3)) 0;
        sin(r(3)) cos(r(3)) 0;
        0 0 1];

    M=[Rx*Ry*Rz t'];
end

