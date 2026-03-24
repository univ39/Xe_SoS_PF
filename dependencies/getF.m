function [F, pinv_time] = getF(E,thresh,mode,useGPU,verb)
%GETF Calculate reconstruction matrix from 
%            E encoding matrix
%       thresh regularization threshold         []
%         mode inversion mode (svd/qr/eig/cholesky)
%       useGPU true or false
%         verb verbosity
% Kylie Yeung 11/2024

if ~exist("thresh","var");thresh=[];end
if ~exist("mode","var");mode=[];end
if isempty(mode);mode='cholesky';end
if ~exist("useGPU","var");useGPU=[];end
if isempty(useGPU);useGPU=canUseGPU;end 
if ~exist("verb","var");verb=[];end
if isempty(verb);verb=1;end 

if useGPU
    E=gpuArray(E);
end


pinv_start=tic;

if strcmpi(mode,'svd')
    if isempty(thresh);thresh=0.95;end
    if verb;fprintf('Using %s with thresh = %.2f \n',mode,thresh);end
    [U,S,V] = svd(E,'econ');
    svdE = diag(S);    
    cumulative_energy = cumsum(svdE) / sum(svdE);
    cutoff = find(cumulative_energy >= thresh, 1);
    invS = 1./svdE;
    invS(cutoff:end) = 0;
    invS = diag(invS);
    F = V*invS*U';
elseif strcmpi(mode,'svdsketch')
    if isempty(thresh);thresh=0.95;end
    if verb;fprintf('Using %s with thresh = %.2f \n',mode,thresh);end
    [U,S,V] = svdsketch(E,eps(class(E))^(1/4),'MaxIterations',4);
    svdE = diag(S);    
    cumulative_energy = cumsum(svdE) / sum(svdE);
    cutoff = find(cumulative_energy >= thresh, 1);
    invS = 1./svdE;
    invS(cutoff:end) = 0;
    invS = diag(invS);
    F = V*invS*U';
elseif strcmpi(mode,'qr')
    if isempty(thresh);thresh=0.2;end
    if verb;fprintf('Using %s with thresh = %.2f \n',mode,thresh);end
    F = pinvUsingQR(E, thresh);
elseif strcmpi(mode,'eig')
    if isempty(thresh);thresh=0.01;end
    if verb;fprintf('Using %s with thresh = %.2f \n',mode,thresh);end
    if size(E,1)>size(E,2)
        EHE=E'*E;
        [V,D]=eig(EHE);
        diagD=diag(D);
        iCUT=nnz(cumsum(diagD)<thresh*sum(diagD));
        iEHE=V(:,iCUT:end)*diag(1./diagD(iCUT:end))*V(:,iCUT:end)';
        % iEHE=V(:,1:iCUT)*diag(1./diagD(1:iCUT))*V(:,1:iCUT)';
        F=iEHE*E';
    else
        EEH=E*E';
        [V,D]=eig(EEH);
        diagD=diag(D);
        iCUT=nnz(cumsum(diagD)<thresh*sum(diagD));
        iEEH=V(:,iCUT:end)*diag(1./diagD(iCUT:end))*V(:,iCUT:end)';
        F=E'*iEEH;
    end    
elseif strcmpi(mode,'cholesky')
    if isempty(thresh);thresh=0.05;end
    if verb;fprintf('Using %s with thresh = %.2f \n',mode,thresh);end
    if size(E,1)>size(E,2) %overdetermined
        lambda=thresh*norm(E,1);

        EHE=E'*E+lambda*eye(size(E,2));

        L = chol(EHE, 'lower');
        L_inv = inv(L);
        iEHE = L_inv' * L_inv;
        F=iEHE*E';
    else
        lambda=thresh*norm(E,1);

        EEH=E*E'+lambda*eye(size(E,1));
        L = chol(EEH, 'lower');
        L_inv = inv(L);
        iEEH = L_inv' * L_inv;
        F=E'*iEEH;
    end   
elseif strcmpi(mode,'block_cholesky')
    if isempty(thresh);thresh=0.05;end
    if verb;fprintf('Using %s with thresh = %.2f \n',mode,thresh);end
    if size(E,1)>size(E,2) %overdetermined
        lambda=thresh*norm(E,1);
        EHE=E'*E+lambda*eye(size(E,2));
        
        L = blockchol(EHE);
        L_inv = inv(L);
        iEHE = L_inv * L_inv';
        F=iEHE*E';
    else
        lambda=thresh*norm(E,1);

        EEH=E*E'+lambda*eye(size(E,1));
        L = blockchol(EEH);
        L_inv = inv(L);
        iEEH = L_inv * L_inv';
        F=E'*iEEH;
    end  
else
    error('Mode %s not found \n',mode)
end
pinv_time=toc(pinv_start);
F=gather(F);
end