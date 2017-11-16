function Po = Tolga_movepoints(M,P)
    % This function MOVEPOINTS will transform a N x 3 array of 3D points
    % with a 4x4  affine transformation matrix.
    %
    % inputs,
    %    P : M x 3 array with XYZ points
    %    M : Affine transformation matrix 4 x 4
    %
    % outputs,
    %    PO : the transformed points

    Po=P;
    Po(:,1)=P(:,1)*M(1,1)+P(:,2)*M(1,2)+P(:,3)*M(1,3)+M(1,4);
    Po(:,2)=P(:,1)*M(2,1)+P(:,2)*M(2,2)+P(:,3)*M(2,3)+M(2,4);
    Po(:,3)=P(:,1)*M(3,1)+P(:,2)*M(3,2)+P(:,3)*M(3,3)+M(3,4);

    if (size(P,2)>5)
        Po(:,4)=P(:,4)*M(1,1)+P(:,5)*M(1,2)+P(:,6)*M(1,3);
        Po(:,5)=P(:,4)*M(2,1)+P(:,5)*M(2,2)+P(:,6)*M(2,3);
        Po(:,6)=P(:,4)*M(3,1)+P(:,5)*M(3,2)+P(:,6)*M(3,3);
    end
end

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