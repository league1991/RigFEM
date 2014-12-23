tangentMat = sparse(399,399);
for i = 1:size(t,1)
    tangentMat(t(i,1), t(i,2)) = t(i,3);
end