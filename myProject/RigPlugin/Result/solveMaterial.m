function [x] = solveMaterial(dataFileName, inputFactor, lambda, smoothFactor)
    eval(dataFileName);

	% for general force fit 
	% reducedElementGF = [gf1 gf2 ... gfn]    	n = number of elements
	% gfi = [gfi1 gfi2 ... gfim]' 				m = frame number
	% gfij = [gfij1 gfij2 ... gfijl] 			l = number of parameters, gfijk = the projection into param k of general force of element i, in frame j
    
	% for general force derivative fit	
	% reducedElementGF = [gf1 gf2 ... gfi ... gfn]    	n = number of elements
	% gfi = [gfi1 gfi2 ... gfij ... gfim]' 				m = frame number
	% gfij = [gfij11 gfij21 ... gfijkl ... gfijpp] 		p = number of parameters, gfijkl = hessian(k,l) of element i, in frame j, note that they are arranged by column major order
    
	nRepeat = size(reducedElementGF, 1) / nParam;    % for GF fit, nRepeat = nFrame; for GF derivative fit nRepeat = nFrame * nParam
    nElement=size(reducedElementGF, 2);
	
	% normalize weights, so the fitting result will get converged when data increased
	nFrame = reducedElementGFCount; %size(reducedElementGF, 1) / 9;
	normalizeWeight = sqrt(nFrame);
	lambda = lambda * normalizeWeight;
	smoothFactor = smoothFactor * normalizeWeight;

    reducedGF = sum(reducedElementGF, 2);
    factor    = inputFactor';

    b0        = repmat(factor, [nRepeat,1]) .* reducedGF;
	b         = b0;
	A         = sparse(reducedElementGF);

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