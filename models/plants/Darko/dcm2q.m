function q = dcm2q( dcm )
%  DCM2QUAT Convert direction cosine matrix to quaternion.
%   Q = DCM2QUAT( N ) calculates the quaternion, Q, for given
%   direction cosine matrix, N.   Input N is a 3-by-3-by-M matrix of
%   orthogonal direction cosine matrices.  The direction cosine matrix performs the
%   coordinate transformation of a vector in inertial axes to a vector in
%   body axes.  Q returns an M-by-4 matrix containing M quaternions. Q has
%   its scalar number as the first column.
%
%   Examples:
%
%   Determine the quaternion from direction cosine matrix:
%      dcm = [0 1 0; 1 0 0; 0 0 1];
%      q = dcm2quat(dcm)
%
%   Determine the quaternions from multiple direction cosine matrices:
%      dcm        = [ 0 1 0; 1 0 0; 0 0 1];
%      dcm(:,:,2) = [ 0.4330    0.2500   -0.8660; ...
%                     0.1768    0.9186    0.3536; ...
%                     0.8839   -0.3062    0.3536];
%      q = dcm2quat(dcm)
%
%   See also ANGLE2DCM, DCM2ANGLE, QUAT2DCM, QUAT2ANGLE, ANGLE2QUAT.

%   Copyright 2000-2010 The MathWorks, Inc.

q = zeros(4,1);

tr = trace(dcm(:,:));

if (tr > 0)
    sqtrp1 = sqrt( tr + 1.0 );
    
    q(1) = 0.5*sqtrp1;
    q(2) = (dcm(2, 3) - dcm(3, 2))/(2.0*sqtrp1);
    q(3) = (dcm(3, 1) - dcm(1, 3))/(2.0*sqtrp1);
    q(4) = (dcm(1, 2) - dcm(2, 1))/(2.0*sqtrp1);
else
    d = diag(dcm(:,:));
    if ((d(2) > d(1)) && (d(2) > d(3)))
        % max value at dcm(2,2,i)
        sqdip1 = sqrt(d(2) - d(1) - d(3) + 1.0 );
        
        q(3) = 0.5*sqdip1;
        
        if ( sqdip1 ~= 0 )
            sqdip1 = 0.5/sqdip1;
        end
        
        q(1) = (dcm(3, 1) - dcm(1, 3))*sqdip1;
        q(2) = (dcm(1, 2) + dcm(2, 1))*sqdip1;
        q(4) = (dcm(2, 3) + dcm(3, 2))*sqdip1;
    elseif (d(3) > d(1))
        % max value at dcm(3,3,i)
        sqdip1 = sqrt(d(3) - d(1) - d(2) + 1.0 );
        
        q(4) = 0.5*sqdip1;
        
        if ( sqdip1 ~= 0 )
            sqdip1 = 0.5/sqdip1;
        end
        
        q(1) = (dcm(1, 2) - dcm(2, 1))*sqdip1;
        q(2) = (dcm(3, 1) + dcm(1, 3))*sqdip1;
        q(3) = (dcm(2, 3) + dcm(3, 2))*sqdip1;
    else
        % max value at dcm(1,1,i)
        sqdip1 = sqrt(d(1) - d(2) - d(3) + 1.0 );
        
        q(2) = 0.5*sqdip1;
        
        if ( sqdip1 ~= 0 )
            sqdip1 = 0.5/sqdip1;
        end
        
        q(1) = (dcm(2, 3) - dcm(3, 2))*sqdip1;
        q(3) = (dcm(1, 2) + dcm(2, 1))*sqdip1;
        q(4) = (dcm(3, 1) + dcm(1, 3))*sqdip1;
    end
    
end