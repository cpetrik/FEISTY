% #include "GAD_OPTIONS.h"
%
% CBOP
% C !ROUTINE: GAD_CALC_RHS
%
% C !INTERFACE: ==========================================================
%       SUBROUTINE GAD_CALC_RHS(
%      I           bi,bj,iMin,iMax,jMin,jMax,k,kM1,kUp,kDown,
%      I           xA, yA, maskUp, uFld, vFld, wFld,
%      I           uTrans, vTrans, rTrans, rTransKp1,
%      I           diffKh, diffK4, KappaR, diffKr4, TracerN, TracAB,
%      I           deltaTLev, trIdentity,
%      I           advectionScheme, vertAdvecScheme,
%      I           calcAdvection, implicitAdvection, applyAB_onTracer,
%      I           trUseDiffKr4, trUseGMRedi, trUseKPP, trUseSmolHack,
%      O           fZon, fMer,
%      U           fVerT, gTracer,
%      I           myTime, myIter, myThid )
%
% C !DESCRIPTION:
% C Calculates the tendency of a tracer due to advection and diffusion.
% C It calculates the fluxes in each direction indepentently and then
% C sets the tendency to the divergence of these fluxes. The advective
% C fluxes are only calculated here when using the linear advection schemes
% C otherwise only the diffusive and parameterized fluxes are calculated.
% C
% C Contributions to the flux are calculated and added:
% C \begin{equation*}
% C {\bf F} = {\bf F}_{adv} + {\bf F}_{diff} +{\bf F}_{GM} + {\bf F}_{KPP}
% C \end{equation*}
% C
% C The tendency is the divergence of the fluxes:
% C \begin{equation*}
% C G_\theta = G_\theta + \nabla \cdot {\bf F}
% C \end{equation*}
% C
% C The tendency is assumed to contain data on entry.
%
% C !USES: ===============================================================
%       IMPLICIT NONE
% #include "SIZE.h"
% #include "EEPARAMS.h"
% #include "PARAMS.h"
% #include "GRID.h"
% #include "SURFACE.h"
% #include "GAD.h"
% #ifdef ALLOW_AUTODIFF
% # include "AUTODIFF_PARAMS.h"
% #endif /* ALLOW_AUTODIFF */
%
% C !INPUT PARAMETERS: ===================================================
% C bi, bj           :: tile indices
% C iMin, iMax       :: for called routines, to get valid output "gTracer"
% C jMin, jMax       ::                      over this range of indices
% C k                :: vertical index
% C kM1              :: =k-1 for k>1, =1 for k=1
% C kUp              :: index into 2 1/2D array, toggles between 1|2
% C kDown            :: index into 2 1/2D array, toggles between 2|1
% C xA, yA           :: areas of X and Y face of tracer cells
% C maskUp           :: 2-D array for mask at W points
% C uFld, vFld, wFld :: Local copy of velocity field (3 components)
% C uTrans, vTrans   :: 2-D arrays of volume transports at U,V points
% C rTrans           :: 2-D arrays of volume transports at W points
% C rTransKp1        :: 2-D array of volume trans at W pts, interf k+1
% C diffKh           :: horizontal diffusion coefficient
% C diffK4           :: horizontal bi-harmonic diffusion coefficient
% C KappaR           :: 2-D array for vertical diffusion coefficient, interf k
% C diffKr4          :: 1-D array for vertical bi-harmonic diffusion coefficient
% C TracerN          :: tracer field @ time-step n (Note: only used
% C                     if applying AB on tracer field rather than on tendency gTr)
% C TracAB           :: current tracer field (@ time-step n if applying AB on gTr
% C                     or extrapolated fwd in time to n+1/2 if applying AB on Tr)
% C trIdentity       :: tracer identifier (required for KPP,GM)
% C advectionScheme  :: advection scheme to use (Horizontal plane)
% C vertAdvecScheme  :: advection scheme to use (Vertical direction)
% C calcAdvection    :: =False if Advec computed with multiDim scheme
% C implicitAdvection:: =True if vertical Advec computed implicitly
% C applyAB_onTracer :: apply Adams-Bashforth on Tracer (rather than on gTr)
% C trUseDiffKr4     :: true if this tracer uses vertical bi-harmonic diffusion
% C trUseGMRedi      :: true if this tracer uses GM-Redi
% C trUseKPP         :: true if this tracer uses KPP
% C trUseSmolHack    :: true if this tracer uses Smolarkiewicz-Hack to remain > 0
% C myTime           :: current time
% C myIter           :: iteration number
% C myThid           :: thread number
%       INTEGER bi,bj,iMin,iMax,jMin,jMax
%       INTEGER k,kUp,kDown,kM1
%       _RS xA    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RS yA    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RS maskUp(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RL uFld  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RL vFld  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RL wFld  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RL uTrans(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RL vTrans(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RL rTrans(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RL rTransKp1(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RL diffKh, diffK4
%       _RL KappaR(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RL diffKr4(Nr)
%       _RL TracerN(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
%       _RL TracAB (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
%       _RL deltaTLev(Nr)
%       INTEGER trIdentity
%       INTEGER advectionScheme, vertAdvecScheme
%       LOGICAL calcAdvection
%       LOGICAL implicitAdvection, applyAB_onTracer
%       LOGICAL trUseDiffKr4, trUseGMRedi, trUseKPP, trUseSmolHack
%       _RL     myTime
%       INTEGER myIter, myThid
%
% C !OUTPUT PARAMETERS: ==================================================
% C gTracer          :: tendency array
% C fZon             :: zonal flux
% C fMer             :: meridional flux
% C fVerT            :: 2 1/2D arrays for vertical advective flux
%       _RL gTracer(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
%       _RL fZon  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RL fMer  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RL fVerT (1-OLx:sNx+OLx,1-OLy:sNy+OLy,2)
%
% C !FUNCTIONS:       ====================================================
% #ifdef ALLOW_DIAGNOSTICS
%       CHARACTER*4 GAD_DIAG_SUFX
%       EXTERNAL    GAD_DIAG_SUFX
% #endif /* ALLOW_DIAGNOSTICS */
%
% C !LOCAL VARIABLES: ====================================================
% C i,j              :: loop indices
% C df4              :: used for storing del^2 T for bi-harmonic term
% C af               :: advective flux
% C df               :: diffusive flux
% C localT           :: local copy of tracer field
% C locABT           :: local copy of (AB-extrapolated) tracer field
%       INTEGER i,j
%       _RS maskLocW(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RS maskLocS(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RL df4   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RL af    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RL df    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RL localT(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RL locABT(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
%       _RL advFac, rAdvFac


for j=1-OLy,sNy+OLy
    for i=1-OLx,sNx+OLx
        fZon(i,j)      = 0.0;
        fMer(i,j)      = 0.0;
        fVerT(i,j,kUp) = 0.0;
        df(i,j)        = 0.0;
        df4(i,j)       = 0.0;
    end
end

%--   Make local copy of tracer array
for j=1-OLy,sNy+OLy
    for i=1-OLx,sNx+OLx
        localT(i,j)=TracerN(i,j,k);
        locABT(i,j)=TracerN(i,j,k);
    end
end

% C--   Pre-calculate del^2 T if bi-harmonic coefficient is non-zero
%       IF (diffK4 .NE. 0.) THEN
%        CALL GAD_GRAD_X(bi,bj,k,xA,localT,fZon,myThid)
%        CALL GAD_GRAD_Y(bi,bj,k,yA,localT,fMer,myThid)
%        CALL GAD_DEL2(bi,bj,k,fZon,fMer,df4,myThid)
%       ENDIF

%--   Initialize net flux in X direction
for j=1-OLy,sNy+OLy
    for i=1-OLx,sNx+OLx
        fZon(i,j) = 0.0;
    end
end

%-    Advective flux in X
if (calcAdvection) 
    if (advectionScheme == ENUM_UPWIND_3RD ) 
        CALL GAD_U3_ADV_X( bi,bj,k, uTrans, maskLocW, locABT, af, myThid )
    end
end

for j=1-OLy,sNy+OLy
    for i=1-OLx,sNx+OLx
        fZon(i,j) = fZon(i,j) + af(i,j)
    end
end

%-    Diffusive flux in X
if (diffKh ~= 0) 
    CALL GAD_DIFF_X(bi,bj,k,xA,diffKh,localT,df,myThid)
else
    for j=1-OLy,sNy+OLy
        for i=1-OLx,sNx+OLx
            df(i,j) = 0.0;
        end
    end
end

% C-    Add bi-harmonic diffusive flux in X
%       IF (diffK4 .NE. 0.) THEN
%        CALL GAD_BIHARM_X(bi,bj,k,xA,df4,diffK4,df,myThid)
%       end
%

%     anelastic: advect.fluxes are scaled by rhoFac but hor.diff. flx are not
for j=1-OLy,sNy+OLy
    for i=1-OLx,sNx+OLx
        fZon(i,j) = fZon(i,j) + df(i,j)*rhoFacC(k)
    end
end

%--   Initialize net flux in Y direction
for j=1-OLy,sNy+OLy
    for i=1-OLx,sNx+OLx
        fMer(i,j) = 0.0;
    end
end

%-    Advective flux in Y
if (calcAdvection) 
    if (advectionScheme == ENUM_UPWIND_3RD ) 
        CALL GAD_U3_ADV_Y( bi,bj,k, vTrans, maskLocS, locABT, af, myThid )
    end
end

for j=1-OLy,sNy+OLy
    for i=1-OLx,sNx+OLx
        fMer(i,j) = fMer(i,j) + af(i,j)
    end
end

%-    Diffusive flux in Y
if (diffKh ~= 0) 
    CALL GAD_DIFF_Y(bi,bj,k,yA,diffKh,localT,df,myThid)
else
    for j=1-OLy,sNy+OLy
        for i=1-OLx,sNx+OLx
            df(i,j) = 0.0;
        end
    end
end

% C-    Add bi-harmonic flux in Y
%       IF (diffK4 .NE. 0.) THEN
%        CALL GAD_BIHARM_Y(bi,bj,k,yA,df4,diffK4,df,myThid)
%       ENDIF
%

%     anelastic: advect.fluxes are scaled by rhoFac but hor.diff. flx are not
for j=1-OLy,sNy+OLy
    for i=1-OLx,sNx+OLx
        fMer(i,j) = fMer(i,j) + df(i,j)*rhoFacC(k)
    end
end


% C--   Compute vertical flux fVerT(kUp) at interface k (between k-1 & k):
% C-    Advective flux in R

%--   Divergence of fluxes
% C     Anelastic: scale vertical fluxes by rhoFac and leave Horizontal fluxes unchanged
% C     for Stevens OBC: keep only vertical diffusive contribution on boundaries
for j=1-OLy,sNy+OLy-1
    for i=1-OLx,sNx+OLx-1
        gTracer(i,j,k) = gTracer(i,j,k) ...
            -_recip_hFacC(i,j,k,bi,bj)*recip_drF(k) ...
            *recip_rA(i,j,bi,bj)*recip_deepFac2C(k)*recip_rhoFacC(k) ...
            *( (fZon(i+1,j)-fZon(i,j))*maskInC(i,j,bi,bj) ...
            +(fMer(i,j+1)-fMer(i,j))*maskInC(i,j,bi,bj) ...
            +(fVerT(i,j,kDown)-fVerT(i,j,kUp))*rkSign ...
            -localT(i,j)*( (uTrans(i+1,j)-uTrans(i,j))*advFac ...
            +(vTrans(i,j+1)-vTrans(i,j))*advFac ...
            +(rTransKp1(i,j)-rTrans(i,j))*rAdvFac ...
            )*maskInC(i,j,bi,bj) ...
            )
    end
end


%       END
