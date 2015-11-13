####by Heather E. Wheeler 20150108####
## Modified by Kaanan Shah
## Nov 13, 2015
## script to run glmnet to get CV-R2 and betas for prediction model of gene expression
## using cis genotypes for each gene with GTEx data

date <- Sys.Date()
"%&%" = function(a,b) paste(a,b,sep="")
doTW<-function(tis,chrom,alpha) {
    tisout<-tis
    if (tisout == "Brain-Anteriorcingulatecortex(BA24)") {tisout = "Brain-Anteriorcingulatecortex-BA24"}
    if (tisout == "Brain-Caudate(basalganglia)") {tisout = "Brain-Caudate-basalganglia"}
    if (tisout == "Brain-FrontalCortex(BA9)") {tisout = "Brain-FrontalCortex-BA9"}
    if (tisout == "Brain-Nucleusaccumbens(basalganglia)") {tisout = "Brain-Nucleusaccumbens-basalganglia"}
    if (tisout == "Brain-Putamen(basalganglia)") {tisout = "Brain-Putamen-basalganglia"}
    if (tisout == "Skin-NotSunExposed(Suprapubic)") {tisout = "Skin-NotSunExposed-Suprapubic"}
    if (tisout == "Skin-SunExposed(Lowerleg)") {tisout = "Skin-SunExposed-Lowerleg"}
    
    ###############################################
    ### Directories & Variables
    
    my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
    exp.dir <- my.dir %&% "gtex-rnaseq/ind-tissues-from-nick/"
    annot.dir <- my.dir %&% "gtex-annot/"
    gt.dir <- my.dir %&% "gtex-genotypes/"
    en.dir <- "/group/im-lab/nas40t2/kaanan/PrediXcan/betas/Resultsbychromosome/"
    
    k <- 10 ### k-fold CV
    n <- 1 #number of k-fold CV replicates, remove nrep loop for this implementation
    #tis <- args[3]
    
    ##alpha = The elasticnet mixing parameter, with 0≤α≤ 1. The penalty is defined as
    #(1-α)/2||β||_2^2+α||β||_1.
    #alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.
    
    gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein.chr" %&% chrom ##genes split into 22 files
    snpset <- "hapmapSnpsCEU"
    
    ################################################
    ### Functions & Libraries
    
    library(glmnet)
    
    ################################################
    rpkmid <- exp.dir %&% "GTEx_PrediXmod." %&% tis %&% ".exp.adj.15PEERfactors.3PCs.gender.ID.list"
    rpkmgene <- exp.dir %&% "GTEx_PrediXmod." %&% tis %&% ".exp.adj.15PEERfactors.3PCs.gender.GENE.list"
    rpkmfile <- exp.dir %&% "GTEx_PrediXmod." %&% tis %&% ".exp.adj.15PEERfactors.3PCs.gender.IDxGENE.RDS"
    
    expid <- scan(rpkmid,"character")
    geneid <- scan(rpkmgene,"character")
    expdata <- readRDS(rpkmfile)
    #expdata <- matrix(expdata, ncol = length(geneid), byrow=TRUE)
    rownames(expdata) <- expid
    colnames(expdata) <- geneid
    
    t.expdata <- expdata 
    
    gencode <- read.table(gencodefile)
    rownames(gencode) <- gencode[,5]
    t.expdata <- t.expdata[,intersect(colnames(t.expdata),rownames(gencode))] ###pull gene expression data w/gene info
    
    expsamplelist <- rownames(t.expdata) ###samples with exp data###
    
    bimfile <- gt.dir %&% "GTEx_Analysis_2014-06-13.hapmapSnpsCEU.chr" %&% chrom %&% ".bim" ###get SNP position information###
    bim <- read.table(bimfile)
    rownames(bim) <- bim$V2
    
    famfile <- gt.dir %&% "GTEx_Analysis_2014-06-13.hapmapSnpsCEU.ID.list"
    fam <- scan(famfile,"character")
    samplelist <- intersect(fam,expsamplelist)
    
    exp.w.geno <- t.expdata[samplelist,] ###get expression of samples with genotypes###
    explist <- colnames(exp.w.geno)
    
    gtfile <- gt.dir %&% 'GTEx_Analysis_2014-06-13.hapmapSnpsCEU.chr' %&% chrom %&% '.mldose.gz'
    gtX <- read.table(gtfile)
    a<-gtX[,3:dim(gtX)[2]]
    gtX <- as.matrix(a)
    colnames(gtX) <- bim$V2
    rownames(gtX) <- fam
    X <- gtX[samplelist,]
    
    set.seed(42)
    groupid <- sample(1:10,length(samplelist),replace=TRUE) ##need to use same folds to compare alphas
    
    resultsarray <- array(0,c(length(explist),9))
    dimnames(resultsarray)[[1]] <- explist
    resultscol <- c("gene","alpha","cvm","lambda.iteration","lambda.min","n.snps","R2","pval","genename")
    dimnames(resultsarray)[[2]] <- resultscol
    workingbest <- en.dir %&% "working_TW_" %&% tisout %&% "_exp_" %&% k %&% "-foldCV_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% "_chr" %&% chrom %&% "_" %&% date %&% ".txt"
    write(resultscol,file=workingbest,ncolumns=9,sep="\t")
    
    weightcol = c("gene","SNP","refAllele","effectAllele","beta","genename")
    workingweight <- en.dir %&% "TW_" %&% tisout %&% "_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% "_weights_chr" %&% chrom %&% "_" %&% date %&% ".txt"
    write(weightcol,file=workingweight,ncol=6,sep="\t")
    for(i in 1:length(explist)){
        cat(i,"/",length(explist),"\n")
        gene <- explist[i]
        geneinfo <- gencode[gene,]
        chr <- geneinfo[1]
        c <- substr(chr$V1,4,5)
        start <- geneinfo$V3 - 1e6 ### 1Mb lower bound for cis-eQTLS
        end <- geneinfo$V4 + 1e6 ### 1Mb upper bound for cis-eQTLs
        chrsnps <- subset(bim,bim[,1]==c) ### pull snps on same chr
        cissnps <- subset(chrsnps,chrsnps[,4]>=start & chrsnps[,4]<=end) ### pull cis-SNP info
        cisgenos <- X[,intersect(colnames(X),cissnps[,2])] ### pull cis-SNP genotypes
        if(is.null(dim(cisgenos))){
            bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
        }else{
            minorsnps <- subset(colMeans(cisgenos), colMeans(cisgenos,na.rm=TRUE)>0) ###pull snps with at least 1 minor allele###
            minorsnps <- names(minorsnps)
            cisgenos <- cisgenos[,minorsnps]
            #cisgenos <- scale(cisgenos, center=T, scale=T) ## scaling of genotypes
            #cisgenos[is.na(cisgenos)] <- 0
            if(is.null(dim(cisgenos)) | dim(cisgenos)[2] == 0){###effectively skips genes with <2 cis-SNPs
                bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
            }else{
                
                exppheno <- exp.w.geno[,gene] ### pull expression data for gene
                exppheno <- scale(exppheno, center=T, scale=T)  ###need to scale for fastLmPure to work properly
                exppheno[is.na(exppheno)] <- 0
                rownames(exppheno) <- rownames(exp.w.geno)
                
                ##run Cross-Validation over alphalist
                fit <- cv.glmnet(cisgenos,exppheno,nfolds=k,alpha=alpha,keep=T,foldid=groupid,parallel=F) ##parallel=T is slower on tarbell, not sure why
                
                fit.df <- data.frame(fit$cvm,fit$lambda,1:length(fit$cvm)) ##pull info to find best lambda
                best.lam <- fit.df[which.min(fit.df[,1]),] # needs to be min or max depending on cv measure (MSE min, AUC max, ...)
                cvm.best = best.lam[,1]
                lambda.best = best.lam[,2]
                nrow.best = best.lam[,3] ##position of best lambda in cv.glmnet output
                
                ret <- as.data.frame(fit$glmnet.fit$beta[,nrow.best]) # get betas from best lambda
                ret[ret == 0.0] <- NA
                bestbetas = as.vector(ret[which(!is.na(ret)),]) # vector of non-zero betas
                names(bestbetas) = rownames(ret)[which(!is.na(ret))]
                
                pred.mat <- fit$fit.preval[,nrow.best] # pull out predictions at best lambda
                
            }
        }
        if(length(bestbetas) > 0){
            res <- summary(lm(exppheno~pred.mat))
            genename <- as.character(gencode[gene,6])
            rsq <- res$r.squared
            pval <- res$coef[2,4]
            
            resultsarray[gene,] <- c(gene, alpha, cvm.best, nrow.best, lambda.best, length(bestbetas), rsq, pval, genename)
            
            
            ### output best shrunken betas for PrediXcan
            bestbetalist <- names(bestbetas)
            bestbetainfo <- bim[bestbetalist,]
            betatable<-as.matrix(cbind(bestbetainfo,bestbetas))
            betafile<-cbind(gene,betatable[,2],betatable[,5],betatable[,6],betatable[,7],genename) ##output "gene","SNP","refAllele","effectAllele","beta"
            write(t(betafile),file=workingweight,ncolumns=6,append=T,sep="\t") # t() necessary for correct output from write() function
            
        }else{
            genename <- as.character(gencode[gene,6])
            resultsarray[gene,1] <- gene
            resultsarray[gene,2:8] <- c(NA,NA,NA,NA,0,NA,NA)
            resultsarray[gene,9] <- genename
        }
        write(resultsarray[gene,],file=workingbest,ncolumns=9,append=T,sep="\t")
    }
    write.table(resultsarray,file=en.dir %&% "TW_" %&% tisout %&% "_exp_" %&% k %&% "-foldCV_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% "_chr" %&% chrom %&% "_" %&% date %&% ".txt",quote=F,row.names=F,sep="\t")
}
