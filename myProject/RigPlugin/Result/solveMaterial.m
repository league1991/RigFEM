function [x] = solveMaterial(dataFileName, inputFactor, lambda)
    eval(dataFileName);

    nFrame = size(reducedElementGF, 1) / nParam;
    nElement=size(reducedElementGF, 2);

    %lambda = 1000;

    reducedGF = sum(reducedElementGF, 2);
    factor    = inputFactor';%[20,1,1]';

    b0        = repmat(factor, [nFrame,1]) .* reducedGF;
    b         = [b0; ones(nElement,1)*lambda];

    AI        = sparse(ones(nElement,1) * lambda);
    AI        = diag(AI);
    A         = [reducedElementGF; AI];

    %x         = A\b;
    %x0        = reducedElementGF \ b0;
    [x,w,info]=nnls_mldivide(A, b);
    plot(x);
    
	saveCmd = sprintf('saveMatrix(''%sMaterial.m'', ''material'', x);', dataFileName);
	eval(saveCmd);
    
    %imagesc([x,x0]);

    %{
    b         = repmat(factor, [nFrame,1]) .* reducedGF;
    b         = [b; ones(nElement,1)*lambda];

    AI        = zeros(nElement, nElement);
    AI(1:nElement, 1:nElement) = lambda;
    A         = [reducedElementGF; AI];

    [x,w,info]=nnls(A, b);
    %}

end