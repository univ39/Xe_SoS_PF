function [F_z,k,dim,E_z] = getFz(k,mtx_acq, mtx_reco,npts)
%GETFZ get the reconstruction matrix in Z
    z = k(:,:,3);[~,i,~]=unique(z);
    Kz = single(z(sort(i)))./0.5*pi./mtx_reco(3).*mtx_acq(3);
    if mod(mtx_reco(3),2)==1
        ZM = single(-(mtx_reco(3)-1)/2:(mtx_reco(3)-1)/2); 
    else
        ZM = single(-(mtx_reco(3)/2):(mtx_reco(3)/2-1));
    end
    E_z = exp(1i*(-Kz(:)*ZM(:).'));
    % if ~isempty(b1)
    %     cos_factor=cumprod(cosd(b1.FA*ones(size(E_z,1),1)));
    %     b1_factor=repmat(cos_factor,1,size(E_z,2));
    %     F_z = pinv(E_z.*b1_factor,0.5);
    % else
        F_z = pinv(E_z,0.5);
    % end


    clear z Kz ZM
    
    % modify trajectory
    k = k(:,1:npts/mtx_acq(3),1:2);dim=2;
end

