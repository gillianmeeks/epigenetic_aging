suppressMessages(library('plink2R'))
suppressMessages(library("optparse"))

option_list = list(
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
  make_option("--weights", action="store", default=NA, type='character',
              help="File listing molecular weight RDat files (must have columns WGT,ID,CHR,P0,P1) [required]"),
  make_option("--weights_dir", action="store", default=NA, type='character',
              help="Path to directory where weight files (WGT column) are stored [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Load in list of weights
# TODO : TEST FOR NO HEADER HERE
wgtlist = read.table(opt$weights,head=T,as.is=T)
chr_ = unique(wgtlist$CHR)

N = nrow(wgtlist)
out.tbl = data.frame("GENE" = character(N), "HSQ" = numeric(N) , "HSQ_PV" = numeric(N), "MODEL" = character(N) , "MODELCV.R2" = character(N), stringsAsFactors=FALSE)

## For each wgt file:
for ( w in 1:nrow(wgtlist) ) {
        #cat( unlist(wgtlist[w,]) , '\n' )
        # Load weights
        wgt.file = paste(opt$weights_dir,"/",wgtlist$WGT[w],sep='')
        out.tbl$GENE[w] <- unlist(strsplit(wgtlist$WGT[w],split=".wgt.RDat"))
        load(wgt.file)
        # Remove NAs (these should not be here)
        wgt.matrix[is.na(wgt.matrix)] = 0

        # which rows have rsq
        row.rsq = grep( "rsq" , rownames(cv.performance) )
        # which rows have p-values
        row.pval = grep( "pval" , rownames(cv.performance) )
        
        # Identify the best model
        mod.best = which.min(apply(cv.performance[row.pval,,drop=F],2,min,na.rm=T))

        if ( length(mod.best) == 0 ) {
                cat( "WARNING : " , unlist(wgtlist[w,]) , " did not have a predictive model ... skipping entirely\n" )
                FAIL.ctr = FAIL.ctr + 1
                next
        }
                if ( sum(wgt.matrix[, mod.best] != 0) == 0 ) {
                cat( "WARNING : " , unlist(wgtlist[w,])[[3]] , names(cv.performance)[ mod.best ] , "had", nrow(wgt.matrix) , "overlapping SNPs, but none with non-zero expression weights, try more SNPS or a different model\n")
		cur.FAIL = TRUE
        }

        # if this is a top1 model, clear out all the other weights
        if ( substr( (colnames(cv.performance))[ mod.best ],1,4) == "top1" ) wgt.matrix[ -which.max(wgt.matrix[,mod.best]^2)  , mod.best] = 0

        # populate the output
        if ( sum(names(wgtlist) == "PANEL") == 1 ) out.tbl$PANEL[w] = wgtlist$PANEL[w]
        if ( exists("hsq") ) {
                out.tbl$HSQ[w] = hsq[1]
        }
        if ( exists("hsq.pv") ) {
                out.tbl$HSQ.PV[w] = hsq.pv[1]
        }
        out.tbl$MODEL[w] = colnames( cv.performance )[ mod.best ]
        out.tbl$MODELCV.R2[w] = paste(format(cv.performance[row.rsq,mod.best],digits=2,trim=T),collapse=',')
}

write.table( format( out.tbl , digits=3 ) , quote=F , row.names=F , sep='\t' , file=opt$out)