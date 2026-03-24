function [mtx_acq,mtx_reco,npts_arm,narms,npts,k,nexc] = getKeyWFVars(wf)
%GETKEYWFVARS Summary of this function goes here
%   Detailed explanation goes here
mtx_acq=wf.mtx;
mtx_reco=mtx_acq;
npts_arm=size(wf.ks,2);
narms=size(wf.phi,2);
npts=npts_arm*narms;
k=wf.k;
nexc=narms*mtx_acq(3);
end

