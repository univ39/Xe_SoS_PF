function progressBar(ind,n,text)
% input iteration and total number required for progress
    if nargin<3;text=[];end
    persistent strCR
    if ind==1;strCR=-1;end

    percent = ind / n;
    numHashes = floor(percent * 50);
    bar = [repmat('#', 1, numHashes), repmat('-', 1, 50 - numHashes)];

    strOut=[text ' ' num2str(ind) ' of ' num2str(n) ' [' bar ']' '(' num2str(percent*100) '%%)'];

    if strCR == -1
        % Don't do carriage return during first run
        fprintf(strOut);
    else
        % Do it during all the other runs
        fprintf([strCR strOut]);
    end
    
    % Update carriage return
    strCR = repmat('\b',1,length(strOut)-1);

    if ind==n;fprintf(['\n']);end

    drawnow limitrate
end

