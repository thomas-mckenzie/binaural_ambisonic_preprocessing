function R_M = getRSH(M, dirs_deg)
%GETRSH Get vector of real orthonormal spherical harmonic values up to order N
%
% Inputs:
%   N:      maximum order of harmonics
%   dirs:   [azimuth_1 elevation_1; ...; azimuth_K elevation_K] angles
%           in degs for each evaluation point, where elevation is the
%           polar angle from the horizontal plane
%
% Outpus:
%   R_N:    [(N+1)^2 x K] matrix of SH values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 10/10/2013
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mdirs = size(dirs_deg, 1);
Mharm = (M+1)^2;

% convert to rads
dirs = dirs_deg*pi/180;
% initialize SH matrix
R_M = zeros(Mharm, Mdirs);
% zero order
R_M(1,:) = 1/sqrt(4*pi);
% higher orders
if M>0 
    idx_R = 1;
    for m=1:M
        
        n = (0:m)';
        % vector of unnormalised associated Legendre functions of current order
        Pmn = legendre(m, sin(dirs(:,2)'));
        % cancel the Condon-Shortley phase from the definition of
        % the Legendre functions to result in signless real SH
        uncondon = (-1).^[n(end:-1:2);n] * ones(1,Mdirs);
        Pmn = uncondon .* [Pmn(end:-1:2, :); Pmn];
        
        % normalisations
        norm_real = sqrt( (2*m+1)*factorial(m-n) ./ (4*pi*factorial(m+n)) );
        
        % convert to matrix, for direct matrix multiplication with the rest
        Nmn = norm_real * ones(1,Mdirs);
        Nmn = [Nmn(end:-1:2, :); Nmn];
        
        CosSin = zeros(2*m+1,Mdirs);
        % zero degree
        CosSin(m+1,:) = ones(1,size(dirs,1));
        % positive and negative degrees
        CosSin(n(2:end)+m+1,:) = sqrt(2)*cos(n(2:end)*dirs(:,1)');
        CosSin(-n(end:-1:2)+m+1,:) = sqrt(2)*sin(n(end:-1:2)*dirs(:,1)');
        Rmn = Nmn .* Pmn .* CosSin;
        
        R_M(idx_R + (1:2*m+1), :) = Rmn;
        idx_R = idx_R + 2*m+1;
    end
end

end
