% Generate a 2D sylindrically/axially symmetric distributoin and apply Abel projection to 1D
% Can be used to generate test imagesc for Inv_Abel_Linear.m
% Each hit assumes a Gaussian distribution.
% Note that rand() excludes 0 and 1

function proj_1D = projAbelGen1D(L,C,R,N)

% L = dimension
% C = distribution centre
% R = distribution radius < L/2
% N = # of hits, should be large.

W = 5; % spread in # of pixels in 1D
sig = 2; % = sqrt(2*variance), size of hit in 2D

projGen = zeros(L,1);
for j=1:N
    D = 0 + R*cos(rand(1)*2*pi);
    projGen(C+floor(D)) = erf((D-floor(D))/sig) + erf((ceil(D)-D)/sig) + projGen(C+floor(D));
    for k = 1:W
        projGen(C+floor(D)+k,1) = erf((ceil(D)-D+k)/sig) - erf((ceil(D)-D+k-1)/sig) + projGen(C+floor(D)+k,1);
        projGen(C+floor(D)-k,1) = erf((D-floor(D)+k)/sig) - erf((D-floor(D)+k-1)/sig) + projGen(C+floor(D)-k,1);
    end
end
proj_1D = projGen/4;
    
end

% To make a 2D distribution with hits at specific coloums Ci and radii Ri:
% proj_2D = zeros(L);
% proj_2D(:,Ci) = projAbelGen1D(L,C,Ri,N);

