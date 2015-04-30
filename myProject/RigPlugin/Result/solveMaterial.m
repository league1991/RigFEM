function [x] = solveMaterial(dataFileName, inputFactor, lambda, smoothFactor)
    eval(dataFileName);

    nFrame = size(reducedElementGF, 1) / nParam;
    nElement=size(reducedElementGF, 2);

    %lambda = 1000;

    reducedGF = sum(reducedElementGF, 2);
    factor    = inputFactor';

    b0        = repmat(factor, [nFrame,1]) .* reducedGF;
	b         = b0;
    % b         = [b0; ones(nElement,1)*lambda; zeros(nElement,1)];

	A         = sparse(reducedElementGF);
    % A         = [reducedElementGF; AI; L*smoothFactor];

	if lambda > 0
		AI = sparse(ones(nElement,1) * lambda);
		AI = diag(AI);
		b = [b; ones(nElement,1)*lambda];
		A = [A; AI];
	end
	if smoothFactor > 0
		b = [b; zeros(nElement,1)];
		A = [A; L*smoothFactor];
	end
	%imagesc(L);
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