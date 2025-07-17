% gad_u3_adv_x routine from MITgcm

clear all
close all

% #include "GAD_OPTIONS.h"
%
% CBOP
% C !ROUTINE: GAD_U3_ADV_X
%
% C !INTERFACE: ==========================================================
%       SUBROUTINE GAD_U3_ADV_X(
%      I           bi,bj,k,
%      I           uTrans, maskLocW,
%      I           tracer,
%      O           uT,
%      I           myThid )
%
% C !DESCRIPTION:
% C Calculates the area integrated zonal flux due to advection of a tracer
% C using upwind biased third-order interpolation (or the $\kappa=1/3$ scheme):
% C \begin{equation*}
% C F^x_{adv} = U \overline{ \theta  - \frac{1}{6} \delta_{ii} \theta }^i
% C                 + \frac{1}{12} |U| \delta_{iii} \theta
% C \end{equation*}
% C Near boundaries, mask all the gradients ==> still 3rd O.
%
%% C !USES: ===============================================================
%       IMPLICIT NONE
% #include "SIZE.h"
% c#include "GRID.h"
% #include "GAD.h"

% From GAD.h
% C loop range for computing vertical advection tendency
% C  iMinAdvR,iMaxAdvR  :: 1rst index (X-dir) loop range for vertical advection
% C  jMinAdvR,jMaxAdvR  :: 2nd  index (Y-dir) loop range for vertical advection
%       INTEGER iMinAdvR, iMaxAdvR, jMinAdvR, jMaxAdvR
% c     PARAMETER ( iMinAdvR = 1-OLx , iMaxAdvR = sNx+OLx )
% c     PARAMETER ( jMinAdvR = 1-OLy , jMaxAdvR = sNy+OLy )
% C- note: we use to compute vertical advection tracer tendency everywhere
% C        (overlap included) as above, but really needs valid tracer tendency
% C        in interior only (as below):
%       PARAMETER ( iMinAdvR = 1 , iMaxAdvR = sNx )
%       PARAMETER ( jMinAdvR = 1 , jMaxAdvR = sNy )

%% Define grid size
% From SIZE.h
% C     Voodoo numbers controlling data layout:
% C     sNx :: Number of X points in tile.
% C     sNy :: Number of Y points in tile.
% C     OLx :: Tile overlap extent in X.
% C     OLy :: Tile overlap extent in Y.
% C     nSx :: Number of tiles per process in X.
% C     nSy :: Number of tiles per process in Y.
% C     nPx :: Number of processes to use in X.
% C     nPy :: Number of processes to use in Y.
% C     Nx  :: Number of points in X for the full domain.
% C     Ny  :: Number of points in Y for the full domain.
% C     Nr  :: Number of points in vertical direction.
sNx =  30;
sNy =  15;
OLx =   2;
OLy =   2;
nSx =   2;
nSy =   4;
nPx =   1;
nPy =   1;
Nx  = sNx*nSx*nPx;
Ny  = sNy*nSy*nPy;
Nr  =   4;

%MITgcm cell ids
[cidx,cidy] = meshgrid(1-OLx:sNx+OLx,sNy+OLy:-1:1-OLy);

%% C !INPUT PARAMETERS: ===================================================
% C  bi,bj        :: tile indices
% C  k            :: vertical level
% C  uTrans       :: zonal volume transport
% C  maskLocW     :: mask (either 0 or 1) at grid-cell western edge
% C  tracer       :: tracer field
% C  myThid       :: my thread Id number
%       INTEGER bi,bj,k
%       _RL uTrans  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RS maskLocW(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RL tracer  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       INTEGER myThid

% Velocity and volume transport
%velocity
%U = 0.1*rand(1-OLx:sNx+OLx,1-OLy:sNy+OLy);
%            (  -1 : 32    ,  -1 : 17    )
U = 0.1*rand(sNx+2*OLx,sNy+2*OLy);
%           (   34    ,   19    )
%Need grid cell width in y-direction, dy
%in MOM5 this is dyte
%median(GRD.dxtn(:)) = 8.9959e+04
%median(GRD.dyte(:)) = 1.0804e+05
%dy = 1.08e+05;
dy=1;
uTrans = U*dy;

% Grid mask
maskLocW = ones(sNx+2*OLx,sNy+2*OLy);

% Tracer
tracer = zeros(sNx+2*OLx,sNy+2*OLy);
tracer(5:10,7:9) = rand;

%% C !LOCAL VARIABLES: ====================================================
% C  i,j                  :: loop indices
% C  Rjm,Rj,Rjp           :: differences at i-1,i,i+1
% C  Rjjm,Rjjp            :: second differences at i-1,i
%       INTEGER i,j
%       _RL Rjm,Rj,Rjp,Rjjm,Rjjp
% CEOP


%% x advection
uT = NaN(sNx+2*OLx,sNy+2*OLy);
%for j=1-OLy:sNy+OLy     %for j=-1:17 (the tile plus the 2-cell buffer in both y-directions)
for j=1:sNy+2*OLy
    uT(1,j)=0.0;            %uT(-1,j)
    uT(2,j)=0.0;            %uT(0,j)
    uT(sNx+2*OLx,j)=0.0;    %uT(32,j)
end
%%
%for j=1-OLy:sNy+OLy         %for j=-1:17 (the tile plus the 2-cell buffer in both y-directions)
for j=1:sNy+2*OLy            %for j=1:19
    %for i=1-OLx+2:sNx+OLx-1 %for i=1:31 (the tile plus 1 cell of buffer in W dir)
    for i=1+OLx:sNx+2*OLx-1  %for i=3:33
        Rjp = (tracer(i+1,j)-tracer( i ,j))*maskLocW(i+1,j);
        Rj  = (tracer( i ,j)-tracer(i-1,j))*maskLocW( i ,j);
        Rjm = (tracer(i-1,j)-tracer(i-2,j))*maskLocW(i-1,j);
        Rjjp = Rjp-Rj;
        Rjjm = Rj-Rjm;
        % F^x_{adv} = U * bar^1{ theta  - 1/6 * delta_{ii}(theta) }
        %             + 1/12 |U| * delta_{iii}(theta)
        uT(i,j) = 0.5 * uTrans(i,j) * ...
            ( tracer(i,j) + tracer(i-1,j) - (1/6)*(Rjjp+Rjjm) ) + ...
            abs(uTrans(i,j)) * 0.5 * (1/6) * (Rjjp-Rjjm);
    end
end

%% C !OUTPUT PARAMETERS: ==================================================
% C  uT           :: zonal advective flux
%       _RL uT      (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%

sum(tracer(:))
sum(uT(:))
% GET NEGATIVE VALUES AND LOSE MASS
% Maybe something wrong with the random velocities?


