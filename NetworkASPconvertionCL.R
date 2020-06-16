suppressMessages(library(dplyr))
suppressMessages(library(plyr))
suppressMessages(library(BoolNet))

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputNetwork',  'i', 1, "character", "REQUIRED: Directed network with three collumns - TF, target, interaction ",
  'outdir',     'o',1, "character", 'Outdir path (default ./)'
), byrow=TRUE, ncol=5);

opt = getopt(spec)

if ( !is.null(opt$help) | is.null(opt$inputASPanswers)) {
  cat("Converting network to ASP readable format")
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

network <- read.csv(file = opt$inputNetwork)
asp <- data.frame()
nb <- union(network$TF,network$target)
asp <- paste('nbnode(',length(union(network$TF,network$target)),').', sep = '')
for(i in 1:length(nb)){
  asp <- rbind(asp,
  	paste('node("',nb[i],'").', 
  	sep = ''))
}
for(i in 1:dim(network)[1]){
  asp <- rbind(asp,
  	paste('in("',network$TF[i],'","',network$target[i],'",',network$interaction[i],").", 
  	sep = ''))
}
for(i in 1:length(nb)){
  asp <- rbind(asp,
  	paste('maxC("',nb[i],'",',dim(network[network$target==nb[i],])[1],').', 
  	sep = ''))
}

write.table(asp, file = opt$outdir, 
	quote = FALSE, 
	row.names = FALSE, 
	col.names = FALSE,
	append = TRUE)