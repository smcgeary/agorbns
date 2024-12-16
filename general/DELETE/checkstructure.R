table = read.table("text.txt",header=FALSE,skip=2)[,1:3]

prob = rep(0,max(table[,1:2]))

for (row in seq(dim(table)[1])) {
	print(row)
	print(table[row,])
	print(table[row,1])
	print(table[row,3])
	print(prob[table[row,1]])
	print(prob[table[row,2]])
	prob[table[row,1]] <- prob[table[row,1]] + table[row,3]^2

	prob[table[row,2]] <- prob[table[row,2]] + table[row,3]^2
	print(prob[table[row,1]])
	print(prob[table[row,2]])

}

print(prob)