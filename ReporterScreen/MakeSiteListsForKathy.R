source("general/general.R")

print(kMirnas)


MakeMirnaSiteSeqList <- function(mirna, experiment="equilibrium", n_constant=5,
                                 sitelist="paperfinal", buffer=FALSE) {
	if (mirna == "miR-1") {
		buffer <- TRUE
	}
	if (mirna == "miR-7-23nt") {
		experiment <- "equilibrium2_nb"
		mirna.str <- "miR-7"
	} else {
		mirna.str <- mirna
	}
	sXe <- SubfunctionCall(SitesXCounts)
	print(sXe)
	sites <- rownames(sXe)
	out <- sapply(sites, GetSiteSeq, mirna=mirna)
	out <- out[1:(length(out) - 1)]
	print(out)
	out
	file.path <- sprintf("variants/site_sequences/%s_%s_sites.txt", mirna.str, sitelist)
	write.table(out, file=file.path, quote=FALSE, col.names=FALSE, sep="\t")
}


sapply(kMirnas, MakeMirnaSiteSeqList)
MakeMirnaSiteSeqList(kMirnas)

