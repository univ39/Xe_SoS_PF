function E=getE(k,mtx_acq,mtx_reco,dim,shift)
if ~exist('shift','var');shift=[];end
if isempty(shift);shift=[0 0 0];end
switch dim
    case 1
        [XM]=single(-(mtx_reco/2-1):mtx_reco/2);   % 2D coordinates
        
        % scale from +-0.5 to +-pi
        Kx=single(k(:,:,1))./0.5*pi./mtx_reco.*mtx_acq;
        
        % Encoding matrix exp(i*k*r)
        E=exp(1i*(Kx(:)*(XM(:)).'));
    case 2
        if length(mtx_acq)>1;mtx_acq=mtx_acq(1);end
        if length(mtx_reco)>1;mtx_reco=mtx_reco(1);end
        [XM,YM]=ndgrid(single(-(mtx_reco/2-1):mtx_reco/2));   % 2D coordinates

        XM=XM+shift(1);
        YM=YM+shift(2);
        
        % scale from +-0.5 to +-pi
        Kx=single(k(:,:,1))./0.5*pi./mtx_reco.*mtx_acq;
        Ky=single(k(:,:,2))./0.5*pi./mtx_reco.*mtx_acq;
        
        % Encoding matrix exp(i*k*r)
        E=exp(1i*(Kx(:)*(XM(:)).'+Ky(:)*(YM(:)).'));
    case 3
        xvec=single(-(mtx_reco(1)/2-1):mtx_reco(1)/2);
        yvec=single(-(mtx_reco(2)/2-1):mtx_reco(2)/2);
        zvec=single(-(mtx_reco(3)/2-1):mtx_reco(3)/2);

        [XM,YM,ZM]=ndgrid(xvec,yvec,zvec);   % 2D coordinates
        
        % scale from +-0.5 to +-pi
        Kx=single(k(:,:,1))./0.5*pi./mtx_reco(1).*mtx_acq(1);
        Ky=single(k(:,:,2))./0.5*pi./mtx_reco(2).*mtx_acq(2);
        Kz=single(k(:,:,3))./0.5*pi./mtx_reco(3).*mtx_acq(3);
        
        % Encoding matrix exp(i*k*r)
        E=exp(-1i*(Kx(:)*(XM(:)).'+Ky(:)*(YM(:)).'+Kz(:)*(ZM(:)).'));
    otherwise, error('dim(=%g)~=1,2, or 3');

end