function saveMatrix(fileName, matName, mat)
fid = fopen(fileName, 'w');
[max_row, max_col] = size(mat);
fprintf(fid, '%s =[ ', matName);
for row = 1:max_row
    for col = 1:max_col
		fprintf(fid, '%15g ', mat(row, col));
	end
	fprintf(fid, '\n', mat(row, max_col));
end
fprintf(fid, '];');
fclose(fid);