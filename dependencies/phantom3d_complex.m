function [BILD_complex] = phantom3d_complex(mtx_reco,PhaseAmp)
%PHATNTOM3D_COMPLEX Generate a shepp-logan phantom with phase variations

if nargin<2;PhaseAmp=pi/2;end

BILD_magnitude=flipud(phantom3d('Modified Shepp-Logan',mtx_reco));
BILD_phase=generateSmoothPhase(mtx_reco,'PhaseAmp',PhaseAmp);
BILD_complex=BILD_magnitude.*exp(1i*BILD_phase);
end

