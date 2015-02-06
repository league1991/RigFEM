function W = solveWeight(dataFileName, weightFileName)

    eval(dataFileName);

    nFrame = size(q,1);
    pos = q;%repmat(initPos, [nFrame,1]) + q; 
    %plotPnt(pos(10,:))
    
    % prepare indices
    nIntPnt = size(intPntIdx,2);    
    nSurfPnt = size(surfPntIdx,2);
    
    % build coef matrix
    % A = [x; y; z], x = [frame1X; frame2X; ... ; frameNX]
    surfIdx = (surfPntIdx-1)*3;
    A = [pos(:, surfIdx+1); pos(:, surfIdx+2); pos(:, surfIdx+3)];
    
    % W = (nSurfPnt, nIntPnt)
    W = zeros(nSurfPnt, nIntPnt);
    for ithIntPnt = 1:nIntPnt
        idx = (intPntIdx(ithIntPnt)-1)*3;
        b = [pos(:,idx+1); pos(:,idx+2); pos(:,idx+3)];
        %x = lsqnonneg(A,b);
        x = nnls(A,b);
        W(:,ithIntPnt) = x;
    end
    imagesc(W);
    size(W);
	
	saveCmd = sprintf('saveMatrix(''%s'', ''weight'', W'');', weightFileName);
	eval(saveCmd);
end

function plotPnt(v)
    nPnt = length(v)/3;
    p = reshape(v,[3, nPnt]);
    scatter3(p(1,:),p(2,:),p(3,:),'.');
    axis equal
end
