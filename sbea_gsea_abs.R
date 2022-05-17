
## taken from the Enrichment browser  package
## GSEA.Ranking metric changed from S2N to absolute S2N

sbea_abs <- function(   
  method = EnrichmentBrowser::sbeaMethods(), 
  se, 
  gs, 
  alpha = 0.05, 
  perm = 1000, 
  padj.method = "none",
  out.file = NULL,
  browse = FALSE,
  assay = "auto", 
  ...)
{   
  # get configuration
  GS.MIN.SIZE <- configEBrowser("GS.MIN.SIZE")
  GS.MAX.SIZE <- configEBrowser("GS.MAX.SIZE")
  FC.COL <-  configEBrowser("FC.COL")
  PVAL.COL <- configEBrowser("PVAL.COL")
  ADJP.COL <-  configEBrowser("ADJP.COL")
  
  # TODO: disentangle DE and EA analysis
  se <- .preprocSE(se)
  se <- .setAssay(method, se, perm, assay)
  
  # data type: ma or rseq?
  is.rseq <- metadata(se)$dataType == "rseq"
  
  # getting gene sets
  if(is(gs, "GeneSetCollection")) gs <- GSEABase::geneIds(gs)
  if(!is.list(gs)) gs <- getGenesets(gs)
  
  # restrict se and gs to intersecting genes
  igenes <- intersect(rownames(se), unique(unlist(gs)))
  if(!length(igenes)) stop("Expression dataset (se)", " and ", 
                           "gene sets (gs) have no gene IDs in common")
  se <- se[igenes,]
  gs <- lapply(gs, function(s) s[s %in% igenes]) 
  lens <- lengths(gs)
  gs <- gs[lens >= GS.MIN.SIZE & lens <= GS.MAX.SIZE]
  
  if(is.character(method))
  { 
    method <- match.arg(method)
    
    # needs conversion of gs list to adj. matrix?
    cmat.methods <- c("ora", "safe", "samgs", "ebm")
    if(method %in% cmat.methods) 
    {
      cmat <- .gs2cmat(gs)
      se <- se[rownames(cmat),]
    }        
    
    ## (1) data type independent (ora, mgsa, ebm)
    ## (2) dedicated RNA-seq mode (camera, roast, gsva)
    ## (3) need transformation (gsea, gsa, padog, safe, samgs, globaltest)
    if(method == "ora") 
    {
      call <- .stdArgs(match.call(), formals())
      exargs <- .matchArgs(.ora, call, list(mode = 1, cmat = cmat))
      exargs$se <- se
      gs.ps <- do.call(.ora, lapply(exargs, eval.parent, n = 2))
    }
    else if(method == "gsea") gs.ps <- .gsea(se, gs, perm)
    else if(method == "padog") gs.ps <- .padog(se, gs, perm)        
    else if(method == "safe") gs.ps <- .ora(2, se, cmat, perm, alpha)
    else if(method %in% c("roast", "camera"))
      gs.ps <- .roast.camera(method, se, gs, perm, rseq = is.rseq)
    else if(method == "gsva") gs.ps <- .gsva(se, gs, rseq = is.rseq)
    else if(method == "gsa") gs.ps <- .gsa(se, gs, perm)
    else if(method == "globaltest") gs.ps <- .globaltest(se, gs, perm)
    else if(method == "samgs") gs.ps <- .samgs(se, cmat, perm, out.file)
    else if(method == "mgsa") gs.ps <- .mgsa(se, gs, alpha, ...)
    else if(method == "ebm") gs.ps <- .ebm(se, cmat)
  }
  else if(is.function(method))
  { 
    call <- .stdArgs(match.call(), formals())
    exargs <- .matchArgs(method, call)
    exargs$se <- se
    exargs$gs <- gs
    gs.ps <- do.call(method, lapply(exargs, eval.parent, n = 2))
  }
  else stop(paste(method, "is not a valid method for sbea"))
  
  res.tbl <- .formatEAResult(gs.ps, padj.method, out.file)
  pcol <- ifelse(padj.method == "none", PVAL.COL, ADJP.COL) 
  res <- list(
    method = method, res.tbl = res.tbl,
    nr.sigs = sum(res.tbl[,pcol] < alpha),
    se = se, gs = gs, alpha = alpha)
  if(browse) eaBrowse(res)
  return(res)
}

#' @rdname eaBrowse
#' @export
.gsRanking <- function(res, signif.only=TRUE)
{
  if(signif.only)
  {
    nr.sigs <- res$nr.sigs
    if(nr.sigs) ranking <- res$res.tbl[seq_len(nr.sigs),]
    else return(NULL)
  }
  else ranking <- res$res.tbl
  return(ranking)
}

.setAssay <- function(method, se, perm, assay = "auto")
{
  # reorder assays
  if(length(assays(se)) > 1 && assay != "auto") se <- .reorderAssays(se, assay)
  
  # data type: ma or rseq?
  data.type <- .detectDataType(assay(se))
  metadata(se)$dataType <- data.type
  
  if(is.function(method)) return(se)
  stopifnot(is.character(method))
  
  # works on the rowData (FC, PVAL) or the assay itself?
  if(method == "ora" && perm == 0) method <- "ora0"
  fdat.methods <- c("ora0", "ebm", "mgsa")
  if(method %in% fdat.methods) return(se) 
  
  is.rseq <- data.type == "rseq"
  is.raw <- method %in% c("camera", "roast", "gsva")
  if(length(assays(se)) == 1)
  {
    if(!is.rseq || is.raw) return(se) 
    se <- EnrichmentBrowser::normalize(se, norm.method = "vst")
  }
  if(assay == "auto") assay <- ifelse(is.rseq && is.raw, "raw", "norm") 
  .reorderAssays(se, assay)    
}

.reorderAssays <- function(se, assay)
{
  ind <- match(assay, names(assays(se)))
  if(is.na(ind)) stop("Expression dataset (se) does not ",
                      "contain an assay named \"", assay, "\"")
  if(ind != 1)
  { 
    ind2 <- setdiff(seq_along(assays(se)), ind)
    assays(se) <- assays(se)[c(ind, ind2)]
  }
  return(se)
}

.formatEAResult <- function(res, padj.method, out.file)
{
  PVAL.COL <- configEBrowser("PVAL.COL")
  ADJP.COL <-  configEBrowser("ADJP.COL")
  
  res.tbl <- data.frame(signif(res, digits=3))
  sorting.df <- res.tbl[,ncol(res.tbl)]
  if(ncol(res.tbl) > 1) 
    sorting.df <- cbind(sorting.df, -res.tbl[,rev(seq_len(ncol(res.tbl)-1))])
  else colnames(res.tbl)[1] <- PVAL.COL 
  res.tbl <- res.tbl[do.call(order, as.data.frame(sorting.df)), , drop=FALSE]
  
  if(padj.method != "none")
    res.tbl[[ADJP.COL]] <- p.adjust(res.tbl[[PVAL.COL]], padj.method)
  
  res.tbl <- DataFrame(rownames(res.tbl), res.tbl)
  colnames(res.tbl)[1] <- configEBrowser("GS.COL")
  rownames(res.tbl) <- NULL
  
  if(!is.null(out.file))
  {
    write.table(res.tbl, 
                file=out.file, quote=FALSE, row.names=FALSE, sep="\t")
    message(paste("Gene set ranking written to", out.file)) 
  }
  return(res.tbl) 
}

.preprocSE <- function(se)
{
  FC.COL <-  configEBrowser("FC.COL")
  PVAL.COL <- configEBrowser("PVAL.COL")
  ADJP.COL <-  configEBrowser("ADJP.COL")
  
  if(is(se, "ExpressionSet")) se <- as(se, "SummarizedExperiment")
  
  if(!(FC.COL %in% colnames(rowData(se))))
    stop(paste("Required rowData column", FC.COL, "not found"))   
  if(!(ADJP.COL %in% colnames(rowData(se))))
    stop(paste("Required rowData column", ADJP.COL, "not found"))   
  
  # dealing with NA's
  se <- se[!is.na(rowData(se)[,FC.COL]),]
  se <- se[!is.na(rowData(se)[,ADJP.COL]),]    
  
  return(se)
}

.gs2cmat <- function(gs)
{
  f <- file()
  sink(file = f)
  cmat <- safe::getCmatrix(gs, as.matrix = TRUE)
  sink()
  close(f)
  return(cmat)
}

.gmt2cmat <- function(gs, features, min.size=0, max.size=Inf)
{
  if(is.character(gs)) gs <- getGenesets(gs)
  # transform gs gmt to cmat
  cmat <- sapply(gs, function(x) features %in% x)
  rownames(cmat) <- features
  
  # restrict to gene sets with valid size
  gs.sizes <- colSums(cmat)
  valid.size <- which((gs.sizes >= min.size) & (gs.sizes <= max.size))
  if(length(valid.size) == 0) stop("No gene set with valid size!")
  cmat <- cmat[, valid.size]
  
  # restrict to genes which are in sets with valid size
  has.set <- which(rowSums(cmat) > 0)
  cmat <- cmat[has.set,]
  
  return(cmat)
}

# deAna as local.stat for safe
.local.deAna <- function (X.mat, y.vec, args.local)
{
  return(function(data, ...) 
  {
    stat <- deAna(expr=data, grp=y.vec,
                  blk=args.local$blk,
                  de.method=args.local$de.method, 
                  stat.only=TRUE)
    return(stat)
  })
}

###
#
# ENRICHMENT METHODS
#
###

.rseqSBEA <- function(method, se, cmat, perm, alpha)
{
  assign("se", se, envir=.GlobalEnv)
  assign("local.deAna", local.deAna, envir=.GlobalEnv)
  de.method <- grep(".STAT$", colnames(rowData(se)), value=TRUE)
  de.method <- sub(".STAT$",  "", de.method)
  
  blk <- NULL
  blk.col <- configEBrowser("BLK.COL") 
  if(blk.col %in% colnames(colData(se))) blk <- colData(se)[,blk.col]
  
  args.local <- list(de.method=de.method, blk=blk)
  
  args.global <- list(one.sided=FALSE)
  if(method == "ora")
  {
    global <- "Fisher"
    nr.sigs <- sum(rowData(se)[, configEBrowser("ADJP.COL")] < alpha)
    args.global$genelist.length <- nr.sigs
  }
  else if(method == "safe") global <- "Wilcoxon" 
  else if(method == "gsea") global <- "Kolmogorov"
  else if(method %in% c("samgs", "gsa", "padog"))
  {
    global <- toupper(method)
    global.func <- paste("global", global, sep=".")
    assign(global.func, get(global.func), envir=.GlobalEnv)
    
    if(method == "padog") args.global$gf <- .getGeneFreqWeights(cmat)
  }
  
  x <- assay(se)
  y <- colData(se)[,configEBrowser("GRP.COL")]
  gs.ps <- safe::safe(X.mat=x, y.vec=y, C.mat=cmat,         
                      local="deAna", args.local=args.local,
                      global=global, args.global=args.global, 
                      Pi.mat=perm, alpha=alpha, error="none")
  
  res.tbl <- cbind(
    gs.ps@global.stat, 
    gs.ps@global.stat / colSums(cmat), 
    gs.ps@global.pval)
  
  colnames(res.tbl) <- c("GLOB.STAT", "NGLOB.STAT", configEBrowser("PVAL.COL"))
  
  return(res.tbl)
}

.isSig <- function(rdat, alpha=0.05, beta=1, sig.stat=c("p", "fc", "|", "&"))
{
  FC.COL <- configEBrowser("FC.COL")
  ADJP.COL <- configEBrowser("ADJP.COL")
  
  sig.stat <- sig.stat[1]
  if(grepl("p$", sig.stat))
  {
    if(sig.stat == "p") sig <- rdat[, ADJP.COL] < alpha
    else
    {
      perc <- as.integer(substring(sig.stat, 1, 2))
      p <- rdat[,ADJP.COL]
      names(p) <- rownames(rdat) 
      ordp <- sort(p)
      nr.sig <- round( length(p) * (perc / 100) )
      sigs <- names(ordp)[seq_len(nr.sig)]
      sig <- rownames(rdat) %in% sigs
    }
  }
  else if(grepl("fc$", sig.stat))
  { 
    if(sig.stat == "fc") sig <- abs(rdat[, FC.COL]) > beta
    else
    {
      perc <- as.integer(substring(sig.stat, 1, 2))
      fc <- rdat[,FC.COL]
      names(fc) <- rownames(rdat) 
      ordfc <- fc[order(abs(fc), decreasing=TRUE)]
      nr.sig <- round( length(fc) * (perc / 100) )
      sigs <- names(ordfc)[seq_len(nr.sig)]
      sig <- rownames(rdat) %in% sigs
    }
  }
  else 
  {
    psig <- rdat[, ADJP.COL] < alpha
    fcsig <- abs(rdat[, FC.COL]) > beta
    sig <- do.call(sig.stat, list(psig, fcsig))
  }
  return(sig)
}

# 1 HYPERGEOM ORA
.oraHypergeom <- function(rdat, cmat, 
                          alpha=0.05, beta=1, sig.stat=c("p", "fc", "|", "&"))
{
  # determine sig. diff. exp. genes of se, 
  # corresponds to sample size from urn
  isig <- .isSig(rdat, alpha, beta, sig.stat)
  nr.sigs <- sum(isig)
  
  # determine overlap of sig and set genes for each set
  sig.cmat <- cmat & isig
  
  # white balls observed when drawing nr.sigs balls from the urn
  ovlp.sizes <- colSums(sig.cmat)
  
  # white balls in the urn  (genes in gene set)
  gs.sizes <- colSums(cmat) 
  # black balls in the urn (genes not in gene set)
  uni.sizes <- nrow(rdat) - gs.sizes 
  
  # determine significance of overlap 
  # based on hypergeom. distribution
  gs.ps <- phyper(ovlp.sizes-1, gs.sizes, uni.sizes, nr.sigs, lower.tail=FALSE) 
  
  res.tbl <- cbind(gs.sizes, ovlp.sizes, gs.ps)
  colnames(res.tbl) <- c("NR.GENES", "NR.SIG.GENES", configEBrowser("PVAL.COL"))
  rownames(res.tbl) <- colnames(cmat)
  
  return(res.tbl)
}

# 2 RESAMPL ORA
# 3 SAFE
#
# wrapper to call safe functionality approriately for
# overrepresentation analysis (ORA)
#
# perm=0 will execute traditional hypergeom. ORA
#
# for perm > 0 use
#   mode=1 ... resampl ORA (fisher)
#   mode=2 ... safe default (wilcoxon)
#
#
.ora <- function(mode=2, se, cmat, perm=1000, alpha=0.05, 
                 padj="none", beta=1, sig.stat=c("p", "fc", "|", "&"))
{
  GRP.COL <- configEBrowser("GRP.COL")
  ADJP.COL <- configEBrowser("ADJP.COL")
  
  x <- assay(se)
  y <- colData(se)[, GRP.COL]
  
  # execute hypergeom ORA if no permutations
  rdat <- rowData(se)
  if(perm == 0) res.tbl <- .oraHypergeom(rdat, cmat, alpha, beta, sig.stat)
  # else do resampling using functionality of SAFE
  else{
    # use built-in p-adjusting?
    padj <- switch(padj,
                   BH = "FDR.BH",
                   fdr = "FDR.BH",
                   bonferroni = "FWER.Bonf",
                   holm = "FWER.Holm",
                   BY = "FDR.YB",
                   "none")
    
    # resampl ORA
    if(mode == 1){
      nr.sigs <- sum(.isSig(rdat, alpha, beta, sig.stat))
      args <- list(one.sided=FALSE, genelist.length=nr.sigs)
      
      gs.ps <- safe::safe(X.mat=x, y.vec=y, global="Fisher", C.mat=cmat, 
                          Pi.mat=perm, alpha=alpha, error=padj, args.global=args)
    } 
    # SAFE default                  
    else gs.ps <- safe::safe(X.mat=x, y.vec=y, 
                             C.mat=cmat, Pi.mat=perm, alpha=alpha, error=padj,  args.global = list(one.sided=F))
    pval <- if(padj == "none") gs.ps@global.pval else gs.ps@global.error
    res.tbl <- cbind(
      gs.ps@global.stat, 
      gs.ps@global.stat / colSums(cmat), 
      pval)
    colnames(res.tbl) <- c("GLOB.STAT", "NGLOB.STAT", configEBrowser("PVAL.COL"))
  }
  return(res.tbl)
}

# 4 GSEA
.gsea <- function(
  se, 
  gs.gmt, 
  perm=1000,
  padj="none", 
  out.file=NULL)
{        
  GRP.COL <- configEBrowser("GRP.COL")
  
  # npGSEA
  if(perm==0)
  {
    npGSEA <- pTwoSided <- NULL
    isAvailable("npGSEA", type="software")
    gsc <- .gsList2Collect(gs.gmt)
    res <- npGSEA(x=assay(se), y=se[[GRP.COL]], set=gsc)
    ps <- sapply(res, pTwoSided)
    names(ps) <- names(gs.gmt)
    return(ps)
  }
  
  # build class list
  cls <- list()
  cls$phen <- levels(as.factor(se[[GRP.COL]]))
  cls$class.v <- ifelse(se[[GRP.COL]] == cls$phen[1], 0, 1)
  
  if(is.null(out.file)) 
    out.dir <- configEBrowser("OUTDIR.DEFAULT") 
  else out.dir <- sub("\\.[a-z]+$", "_files", out.file)
  if(!file.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
  
  # use built-in p-adjusting?
  padj <- switch(padj,
                 BH = "fdr",
                 fdr = "fdr",
                 bonferroni = "fwer",
                 holm = "fwer",
                 BY = "fdr",
                 "none")
  
  # run GSEA
  res <- .GSEA_abs(input.ds=as.data.frame(assay(se)), 
              input.cls=cls, gs.db=gs.gmt, nperm=perm,
              padj.method=padj, output.directory=out.dir)
  
  gs.ps <- S4Vectors::as.matrix(res[,3:5])
  rownames(gs.ps) <- res[,1]
  
  return(gs.ps)
}

.samgs <- function(se, cmat, perm, out.file)
{
  GRP.COL <- configEBrowser("GRP.COL")
  
  if(is.null(out.file)) out.dir <- configEBrowser("OUTDIR.DEFAULT")
  else out.dir <- sub("\\.[a-z]+$", "_files", out.file)
  
  if(!file.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
  samt.file <- file.path(out.dir, "samt.RData")
  
  SAMGS(GS = as.data.frame(cmat), DATA = assay(se), 
        cl = as.factor(as.integer(se[[GRP.COL]])), 
        nbPermutations = perm, 
        tstat.file = samt.file)
}

# 5 EBM (_E_mpirical _B_rowns _M_ethod)
.ebm <- function(se, cmat)
{
  empiricalBrownsMethod <- NULL
  isAvailable("EmpiricalBrownsMethod", type="software")
  pcol <-  rowData(se)[, configEBrowser("ADJP.COL")]
  e <- assay(se)
  gs.ps <- apply(cmat, 2, function(s) empiricalBrownsMethod(e[s,], pcol[s]))
  return(gs.ps)
}


# 6 GSA
.gsa <- function(se, gs, perm=1000)
{  
  # setup
  GSA <- NULL
  isAvailable("GSA", type="software")
  
  minsize <- configEBrowser("GS.MIN.SIZE")
  maxsize <- configEBrowser("GS.MAX.SIZE")
  GRP.COL <- configEBrowser("GRP.COL")
  BLK.COL <- configEBrowser("BLK.COL")
  
  # prepare input
  x <- assay(se)
  genenames <- names(se)
  y <- se[[GRP.COL]] + 1
  
  # paired?
  blk <- NULL
  if(BLK.COL %in% colnames(colData(se))) blk <- se[[BLK.COL]] 
  paired <- !is.null(blk)
  resp.type <- ifelse(paired, "Two class paired", "Two class unpaired")
  
  # response vector y need to be differently coded for 2-class paired   
  if(paired)
  {
    y <- blk
    ublk <- unique(blk)
    for(i in seq_along(ublk)) y[blk==ublk[i]] <- c(i,-i)
    y <- as.integer(y)
  }
  
  # run GSA
  res <- GSA(x=x, y=y, nperms=perm, genesets=gs, resp.type=resp.type,
             genenames=genenames, minsize=minsize, maxsize=maxsize)
  
  # format output
  ps <- cbind(res$pvalues.lo, res$pvalues.hi)
  ps <- 2 * apply(ps, 1, min)
  scores <- res$GSA.scores
  res.tbl <- cbind(scores, ps)
  colnames(res.tbl) <- c("SCORE", configEBrowser("PVAL.COL"))
  rownames(res.tbl) <- names(gs)
  
  return(res.tbl)
}

# rseq: GSA maxmean stat as global.stat for safe
.global.GSA <- function(cmat, u, ...)
{
  # SparseM::as.matrix
  isAvailable("SparseM", type="software")
  am <- getMethod("as.matrix", signature="matrix.csr")
  tcmat <- t(am(cmat))
  
  return(
    function(u, cmat2=tcmat) 
    {
      ind.pos <- u > 0 
      
      upos <- u[ind.pos]
      lpos <- rowSums(cmat2[,ind.pos])
      vpos <- as.vector(cmat2[,ind.pos] %*% upos) / lpos
      vpos <- sapply(vpos, function(x) ifelse(is.na(x), 0, x))
      
      uneg <- abs(u[!ind.pos])
      lneg <- rowSums(cmat2[,!ind.pos])
      vneg <- as.vector(cmat2[,!ind.pos] %*% uneg) / lneg
      vneg <- sapply(vneg, function(x) ifelse(is.na(x), 0, x))
      
      mm <- apply(cbind(vpos, vneg), 1, max)
      return(mm)
    }
  )
}

# 7 PADOG
.padog <- function(se, gs, perm=1000)
{
  padog <- NULL
  isAvailable("PADOG", type="software")
  
  grp <- se[[configEBrowser("GRP.COL")]]
  grp <- ifelse(grp == 0, "c", "d") 
  
  blk <- NULL
  BLK.COL <- configEBrowser("BLK.COL")
  if(BLK.COL %in% colnames(colData(se))) 
    blk <- make.names(colData(se)[,BLK.COL]) 
  paired <- !is.null(blk)
  
  nmin <- configEBrowser("GS.MIN.SIZE")
  perm <- as.numeric(perm) 
  
  res <- padog(assay(se), group=grp, 
               paired=paired, block=blk, gslist=gs, Nmin=nmin, NI=perm)
  
  res.tbl <- res[, c("meanAbsT0", "padog0", "PmeanAbsT", "Ppadog")]
  colnames(res.tbl) <- c("MEAN.ABS.T0", 
                         "PADOG0", "P.MEAN.ABS.T",  configEBrowser("PVAL.COL"))
  rownames(res.tbl) <- as.vector(res[,"ID"]) 
  return(res.tbl)
}


# compute gene frequencies across genesets
.getGeneFreqWeights <- function(cmat)
{
  gf <- rowSums(cmat)
  if (!all(gf == 1)) 
  {
    q99 <- quantile(gf, 0.99)
    m3sd <- mean(gf) + 3 * sd(gf)
    if(q99 > m3sd) gf[gf > q99] <- q99
    gff <- function(x) 1 + ((max(x) - x)/(max(x) - min(x)))^0.5
    gf <- gff(gf)
  } 
  else 
  {
    gf <- rep(1, nrow(cmat))
    names(gf) <- rownames(cmat)
  }
  return(gf)
}

# rseq: PADOG weighted mean as global.stat for safe
.global.PADOG <- function(cmat, u, args.global)
{
  # SparseM::as.matrix
  isAvailable("SparseM", type="software")
  #pos <- grep("SparseM", search())
  am <- getMethod("as.matrix", signature="matrix.csr")#, where=pos)
  cmat <- t(am(cmat))
  gs.size <- rowSums(cmat) 
  
  return(
    function(u, cmat2=cmat, gf=args.global$gf, gs.size=rowSums(cmat)) 
    {
      wu <- abs(u) * gf
      return(as.vector(cmat2 %*% wu) / gs.size)
    }
  )
}

# 8a MGSA 
.mgsa <- function(se, gs, alpha=0.05, beta=1, sig.stat=c("p", "fc", "|", "&"))
{
  mgsa <- setsResults <- NULL
  isAvailable("mgsa", type = "software")
  
  # extract significant (DE) genes
  isig <- .isSig(rowData(se), alpha, beta, sig.stat)
  obs <- rownames(se)[isig]
  pop <- rownames(se)
  
  # run mgsa
  res <- mgsa(o=obs, sets=gs, population=pop)
  res <- setsResults(res)[,1:3]
  res[,3] <- 1 - res[,3]
  colnames(res)[3] <- configEBrowser("PVAL.COL")
  return(res)
}

# 8b GLOBALTEST
.globaltest <- function(se, gs, perm=1000)
{
  gt <- NULL
  isAvailable("globaltest", type="software")
  
  grp <- colData(se)[, configEBrowser("GRP.COL")]
  names(assays(se))[1] <- "exprs"
  se <- as(se, "ExpressionSet")
  res <- gt(grp, se, subsets=gs, permutations=perm)
  res <- res@result[,2:1]
  colnames(res) <- c("STAT", configEBrowser("PVAL.COL"))
  return(res)
}

# 9 ROAST
# 10 CAMERA
.roast.camera <- function(method=c("roast", "camera"), se, gs, perm=1000, rseq=FALSE)
{
  method <- match.arg(method)
  
  # design matrix
  grp <- colData(se)[, configEBrowser("GRP.COL")]
  blk <- NULL
  BLK.COL <- configEBrowser("BLK.COL")
  if(BLK.COL %in% colnames(colData(se))) blk <- colData(se)[,BLK.COL]
  
  group <- factor(grp)
  paired <- !is.null(blk)
  f <- "~" 
  if(paired) 
  {   
    block <- factor(blk)
    f <- paste0(f, "block + ") 
  }   
  f <- formula(paste0(f, "group"))
  design <- model.matrix(f)
  
  y <- assay(se)
  # rseq data
  if(rseq)
  {
    y <- edgeR::DGEList(counts=y,group=grp)
    y <- edgeR::calcNormFactors(y)
    y <- edgeR::estimateDisp(y, design)
  }
  
  # set gene sets
  gs.index <- limma::ids2indices(gs, rownames(se))
  
  # run roast / camera
  if(method == "roast")
    res <- limma::mroast(y, gs.index, design, 
                         nrot=perm, adjust.method="none", sort="none")
  else res <- limma::camera(y, gs.index, design, sort=FALSE)
  res <- res[,c("NGenes", "Direction", "PValue")]
  colnames(res) <- c("NR.GENES", "DIR", configEBrowser("PVAL.COL"))
  res[,"DIR"] <- ifelse(res[,"DIR"] == "Up", 1, -1)
  
  return(res)
}


# 11 GSVA
.gsva <- function(se, gs, rseq=FALSE)
{
  gsva <- NULL
  isAvailable("GSVA", type="software")
  
  # compute GSVA per sample enrichment scores
  kcdf <- ifelse(rseq, "Poisson", "Gaussian")
  es <- gsva(expr=assay(se), gset.idx.list=gs, kcdf=kcdf)
  
  # set design matrix
  grp <- colData(se)[, configEBrowser("GRP.COL")]
  blk <- NULL
  BLK.COL <- configEBrowser("BLK.COL")
  if(BLK.COL %in% colnames(colData(se))) blk <- colData(se)[,BLK.COL]
  
  group <- factor(grp)
  paired <- !is.null(blk)
  f <- "~"
  if(paired)
  {
    block <- factor(blk)
    f <- paste0(f, "block + ")
  }
  f <- formula(paste0(f, "group"))
  design <- model.matrix(f)  
  
  # fit the linear model to the GSVA enrichment scores
  fit <- limma::lmFit(es, design)
  fit <- limma::eBayes(fit)
  res <- limma::topTable(fit, number=nrow(es), coef="group1", sort.by="none", adjust.method="none")
  
  # process output
  res <- res[,c("t", "P.Value")]
  colnames(res) <- c("t.SCORE", configEBrowser("PVAL.COL"))
  
  return(res)
}


.detectDataType <- function(expr) 
  ifelse(all(.isWholenumber(expr), na.rm=TRUE), "rseq", "ma")

.isWholenumber <- function(x, tol=.Machine$double.eps^0.5) abs(x-round(x)) < tol


.GSEA_abs <- function(
  input.ds, 
  input.cls, 
  gs.db, 
  output.directory = "", 
  reshuffling.type = "sample.labels", 
  nperm = 1000, 
  padj.method=c("none", "fdr", "fwer"),
  weighted.score.type = 1, 
  #gs.size.threshold.min = 25, 
  #gs.size.threshold.max = 500, 
  random.seed = 123456) 
{
  
  if (.Platform$OS.type == "windows") 
  {
    memory.limit(6000000000)
    memory.limit()
  }
  
  # Read input data matrix
  set.seed(seed=random.seed, kind = NULL)
  adjust.param <- 0.5
  dataset <- input.ds
  gene.labels <- row.names(dataset)
  sample.names <- names(dataset)
  A <- data.matrix(dataset)
  cols <- ncol(A)
  rows <- nrow(A)
  
  # Read input class vector
  CLS <- input.cls
  class.labels <- CLS$class.v
  class.phen <- CLS$phen
  phen1 <- class.phen[2] # cases
  phen2 <- class.phen[1] # controls
  
  # sort samples according to phenotype
  col.index <- order(class.labels, decreasing=TRUE)
  class.labels <- class.labels[col.index]
  sample.names <- sample.names[col.index]
  A <- A[, col.index]
  colnames(A) <- sample.names
  
  # Read input gene set database
  #    temp <- readLines(gs.db, warn=FALSE)
  #    max.Ng <- length(temp)
  #    temp.size.G <- sapply(temp, 
  #        function(t) length(unlist(strsplit(t, "\t"))) - 2)
  
  #   max.size.G <- max(temp.size.G)      
  #    gs <- matrix(rep("null", max.Ng*max.size.G), nrow=max.Ng, ncol= max.size.G)
  #    temp.names <- temp.desc <- vector(length = max.Ng, mode = "character")
  #    gs.count <- 1
  #    for (i in seq_len(max.Ng)) 
  #    {
  #        spl <- unlist(strsplit(temp[[i]], "\t"))
  #        gene.set.size <- length(spl) - 2
  #        gs.line <- noquote(spl)
  #        gene.set.name <- gs.line[1] 
  #        gene.set.desc <- gs.line[2] 
  #        gene.set.tags <- 
  #            sapply(seq_len(gene.set.size), function(j) gs.line[j + 2])
  #        existing.set <- is.element(gene.set.tags, gene.labels)
  #        set.size <- sum(existing.set)
  #        if ((set.size >= gs.size.threshold.min) && 
  #              (set.size <= gs.size.threshold.max))
  #        {
  #            temp.size.G[gs.count] <- set.size
  #            gs[gs.count,] <- c(gene.set.tags[existing.set], 
  #                rep(NA, max.size.G - temp.size.G[gs.count]))
  #            temp.names[gs.count] <- gene.set.name
  #            temp.desc[gs.count] <- gene.set.desc
  #            gs.count <- gs.count + 1
  #        }
  #    } 
  
  Ng <- length(gs.db)
  gs.names <- names(gs.db)
  size.G <- sapply(gs.db, length) 
  gs <- matrix(NA, nrow=Ng, ncol=max(size.G))
  for(i in seq_len(Ng)) gs[i, seq_len(size.G[i])] <- gs.db[[i]]
  
  N <- nrow(A)
  Ns <- ncol(A)
  #all.gene.descs <- all.gene.symbols <-  gene.labels[i]
  
  Obs.indicator <- Obs.RES <- matrix(nrow= Ng, ncol=N)
  Obs.ES <- Obs.arg.ES <- Obs.ES.norm <- vector(length = Ng, mode = "numeric")
  
  # Compute observed and random permutation gene rankings
  obs.s2n <- vector(length=N, mode="numeric")
  signal.strength <- tag.frac <- gene.frac <- 
    coherence.ratio <- vector(length=Ng, mode="numeric")
  obs.phi.norm <- matrix(nrow = Ng, ncol = nperm)
  correl.matrix <- obs.correl.matrix <- order.matrix <- 
    obs.order.matrix <- matrix(nrow = N, ncol = nperm)
  
  nperm.per.call <- 100
  n.groups <- nperm %/% nperm.per.call
  n.rem <- nperm %% nperm.per.call
  n.perms <- c(rep(nperm.per.call, n.groups), n.rem)
  n.ends <- cumsum(n.perms)
  n.starts <- n.ends - n.perms + 1
  n.tot <- ifelse(n.rem == 0, n.groups, n.groups + 1)
  
  for (nk in seq_len(n.tot)) 
  {
    call.nperm <- n.perms[nk]
    message(paste("Permutations:", n.starts[nk], "--", n.ends[nk]))
    O <- .GSEA.GeneRanking(A, class.labels, gene.labels, call.nperm)
    order.matrix[,n.starts[nk]:n.ends[nk]] <- O$order.matrix
    obs.order.matrix[,n.starts[nk]:n.ends[nk]] <- O$obs.order.matrix
    correl.matrix[,n.starts[nk]:n.ends[nk]] <- O$s2n.matrix
    obs.correl.matrix[,n.starts[nk]:n.ends[nk]] <- O$obs.s2n.matrix
    rm(O)
  }
  
  message("Processing ...")
  # using median to assign enrichment scores
  obs.s2n <- apply(obs.correl.matrix, 1, median)    
  names(obs.s2n) <- gene.labels
  save(obs.s2n, file=file.path(output.directory, "gsea_s2n.RData"))  
  
  obs.index <- order(obs.s2n, decreasing=TRUE)            
  obs.s2n <- obs.s2n[obs.index]       
  #obs.gene.labels <- gene.labels[obs.index]       
  #obs.gene.descs <- all.gene.descs[obs.index]       
  #obs.gene.symbols <- all.gene.symbols[obs.index]       
  
  for (r in seq_len(nperm)) 
  {
    correl.matrix[, r] <- correl.matrix[order.matrix[,r], r]
    obs.correl.matrix[, r] <- obs.correl.matrix[obs.order.matrix[,r], r]
  }
  
  gene.list2 <- obs.index
  for (i in seq_len(Ng)) 
  {
    gene.set <- gs[i,!is.na(gs[i,])]
    gene.set2 <- match(gene.set, gene.labels)
    GSEA.results <- .GSEA.EnrichmentScore(
      gene.list=gene.list2, gene.set=gene.set2, 
      weighted.score.type=weighted.score.type, correl.vector = obs.s2n)
    Obs.ES[i] <- GSEA.results$ES
    Obs.arg.ES[i] <- GSEA.results$arg.ES
    Obs.RES[i,] <- GSEA.results$RES
    Obs.indicator[i,] <- GSEA.results$indicator
    if (Obs.ES[i] >= 0) 
    {  
      # compute signal strength
      tag.frac[i] <- sum(Obs.indicator[i,seq_len(Obs.arg.ES[i])])/size.G[i]
      gene.frac[i] <- Obs.arg.ES[i]/N
    } 
    else 
    {
      tag.frac[i] <- sum(Obs.indicator[i, Obs.arg.ES[i]:N])/size.G[i]
      gene.frac[i] <- (N - Obs.arg.ES[i] + 1)/N
    }
    signal.strength[i] <- tag.frac[i] * 
      (1 - gene.frac[i]) * (N / (N - size.G[i]))
  }
  
  # Compute enrichment for random permutations 
  phi <- phi.norm <- obs.phi <- matrix(nrow = Ng, ncol = nperm)
  if (reshuffling.type == "sample.labels") 
  { 
    # reshuffling phenotype labels
    for (i in seq_len(Ng)) 
    {
      gene.set <- gs[i,!is.na(gs[i,])]
      gene.set2 <- match(gene.set, gene.labels)
      for (r in seq_len(nperm)) 
      {
        gene.list2 <- order.matrix[,r]
        GSEA.results <- .GSEA.EnrichmentScore2(
          gene.list=gene.list2, gene.set=gene.set2, 
          weighted.score.type=weighted.score.type, 
          correl.vector=correl.matrix[, r])   
        phi[i, r] <- GSEA.results$ES
      }
      obs.gene.list2 <- obs.order.matrix[,1]
      GSEA.results <- .GSEA.EnrichmentScore2(gene.list=obs.gene.list2, 
                                             gene.set=gene.set2, weighted.score.type=weighted.score.type, 
                                             correl.vector=obs.correl.matrix[, nperm])
      obs.phi[i, ] <- GSEA.results$ES
    }
  }
  else if (reshuffling.type == "gene.labels") 
  { 
    # reshuffling gene labels
    for (i in seq_len(Ng)) 
    {
      gene.set <- gs[i,!is.na(gs[i,])]
      gene.set2 <- match(gene.set, gene.labels)
      for (r in seq_len(nperm)) 
      {
        reshuffled.gene.labels <- sample(1:rows)
        GSEA.results <- .GSEA.EnrichmentScore2(
          gene.list=reshuffled.gene.labels, gene.set=gene.set2, 
          weighted.score.type=weighted.score.type, 
          correl.vector=obs.s2n)   
        phi[i, r] <- GSEA.results$ES
      }
      obs.gene.list2 <- obs.order.matrix[,1]
      GSEA.results <- .GSEA.EnrichmentScore2(gene.list=obs.gene.list2, 
                                             gene.set=gene.set2, weighted.score.type=weighted.score.type, 
                                             correl.vector=obs.correl.matrix[, nperm])   
      obs.phi[i, ] <- GSEA.results$ES
    }
  }
  # Compute 3 types of p-values
  padj.method <- match.arg(padj.method)    
  
  
  #
  # Find nominal p-values
  #
  p.vals <- matrix(0, nrow = Ng, ncol = 2)
  for (i in seq_len(Ng)) 
  {
    ind <- phi[i,] >= 0
    pos.phi <- phi[i, ind]
    neg.phi <- phi[i, !ind] 
    ES.value <- Obs.ES[i]
    p.vals[i, 1] <- signif(ifelse(ES.value >= 0,
                                  sum(pos.phi >= ES.value)/length(pos.phi),
                                  sum(neg.phi <= ES.value)/length(neg.phi)), digits=5)
    
    # Rescaling normalization for each gene set null
    pos.m <- mean(pos.phi)
    neg.m <- mean(abs(neg.phi))
    pos.phi <- pos.phi/pos.m
    neg.phi <- neg.phi/neg.m
    for (j in seq_len(nperm))
    { 
      phi.norm[i, j] <- 
        phi[i, j] / ifelse(phi[i, j] >= 0, pos.m, neg.m)
      obs.phi.norm[i, j] <-
        obs.phi[i, j] / ifelse(obs.phi[i, j] >= 0, pos.m, neg.m)
    }
    Obs.ES.norm[i] <- Obs.ES[i] / ifelse(Obs.ES[i] >= 0, pos.m, neg.m)
  }
  
  #
  #
  #    
  
  # Compute FWER p-vals
  if(padj.method == "fwer")
  {
    max.ES.vals.p <- NULL
    max.ES.vals.n <- NULL
    for (j in seq_len(nperm)) 
    {
      ind <- phi.norm[,j] >= 0
      pos.phi <- phi.norm[ind, j]
      neg.phi <- phi[!ind, j] 
      
      if (length(pos.phi)) max.ES.vals.p <- c(max.ES.vals.p, max(pos.phi))
      if (length(neg.phi)) max.ES.vals.n <- c(max.ES.vals.n, min(neg.phi))
    }
    
    for (i in seq_len(Ng)) 
    {
      ES.value <- Obs.ES.norm[i]
      p.vals[i, 2] <- signif(ifelse(ES.value >= 0, 
                                    sum(max.ES.vals.p >= ES.value),
                                    sum(max.ES.vals.n <= ES.value)) / 
                               length(max.ES.vals.p), digits=5)
    }
    p.vals <- p.vals[,2]
  }
  
  
  # Compute FDRs
  if(padj.method == "fdr")
  { 
    NES <- phi.norm.mean <- obs.phi.norm.mean <- phi.norm.median <- 
      obs.phi.norm.median <- phi.norm.mean <- obs.phi.mean <- 
      FDR.mean <- FDR.median <- phi.norm.median.d <- 
      obs.phi.norm.median.d <- vector(length=Ng, mode="numeric")
    
    Obs.ES.index <- order(Obs.ES.norm, decreasing=TRUE)
    Orig.index <- seq(1, Ng)
    Orig.index <- Orig.index[Obs.ES.index]
    Orig.index <- order(Orig.index, decreasing=FALSE)
    Obs.ES.norm.sorted <- Obs.ES.norm[Obs.ES.index]
    gs.names.sorted <- gs.names[Obs.ES.index]
    
    NES <- Obs.ES.norm.sorted
    for (k in seq_len(Ng)) 
    {
      ES.value <- NES[k]
      count.col <- obs.count.col <- vector(length=nperm, mode="numeric")
      for (i in seq_len(nperm)) 
      {
        phi.vec <- phi.norm[,i]
        obs.phi.vec <- obs.phi.norm[,i]
        if (ES.value >= 0) 
        {
          count.col.norm <- sum(phi.vec >= 0)
          obs.count.col.norm <- sum(obs.phi.vec >= 0)
          count.col[i] <- ifelse(count.col.norm > 0, 
                                 sum(phi.vec >= ES.value)/count.col.norm, 0)
          obs.count.col[i] <- ifelse(obs.count.col.norm > 0, 
                                     sum(obs.phi.vec >= ES.value)/obs.count.col.norm, 0)
        } 
        else 
        {
          count.col.norm <- sum(phi.vec < 0)
          obs.count.col.norm <- sum(obs.phi.vec < 0)
          count.col[i] <- ifelse(count.col.norm > 0, 
                                 sum(phi.vec <= ES.value)/count.col.norm, 0)
          obs.count.col[i] <- ifelse(obs.count.col.norm > 0, 
                                     sum(obs.phi.vec <= ES.value)/obs.count.col.norm, 0)
        }
      }
      phi.norm.mean[k] <- mean(count.col)
      obs.phi.norm.mean[k] <- mean(obs.count.col)
      phi.norm.median[k] <- median(count.col)
      obs.phi.norm.median[k] <- median(obs.count.col)
      FDR.mean[k] <- ifelse(phi.norm.mean[k]/obs.phi.norm.mean[k] < 1, 
                            phi.norm.mean[k]/obs.phi.norm.mean[k], 1)
      FDR.median[k] <- ifelse(phi.norm.median[k]/obs.phi.norm.median[k] < 1, 
                              phi.norm.median[k]/obs.phi.norm.median[k], 1)
    }
    
    # adjust q-values
    adjust.FDR.q.val <- FALSE
    if (adjust.FDR.q.val) 
    {
      pos.nes <- sum(NES >= 0)
      min.FDR.mean <- FDR.mean[pos.nes]
      min.FDR.median <- FDR.median[pos.nes]
      for(k in seq(pos.nes - 1, 1, -1)) 
      {
        if(FDR.mean[k] < min.FDR.mean) min.FDR.mean <- FDR.mean[k]
        if(min.FDR.mean < FDR.mean[k]) FDR.mean[k] <- min.FDR.mean
      }
      neg.nes <- pos.nes + 1
      min.FDR.mean <- FDR.mean[neg.nes]
      min.FDR.median <- FDR.median[neg.nes]
      for (k in seq(neg.nes + 1, Ng)) 
      {
        if(FDR.mean[k] < min.FDR.mean) min.FDR.mean <- FDR.mean[k]
        if (min.FDR.mean < FDR.mean[k]) FDR.mean[k] <- min.FDR.mean
      }
    }   
    
    obs.phi.norm.mean.sorted <- obs.phi.norm.mean[Orig.index]
    phi.norm.mean.sorted <- phi.norm.mean[Orig.index]
    FDR.mean.sorted <- FDR.mean[Orig.index]
    FDR.median.sorted <- FDR.median[Orig.index]
    
    p.vals <- FDR.mean.sorted
  }
  
  
  #    #   Compute global statistic
  #    glob.p.vals <- vector(length=Ng, mode="numeric")
  #    NULL.pass <- OBS.pass <- vector(length=nperm, mode="numeric")
  #    for (k in seq_len(Ng)) 
  #    {
  #        if (NES[k] >= 0) 
  #        {
  #            for (i in seq_len(nperm)) 
  #            {
  #                NULL.pos <- sum(phi.norm[,i] >= 0)            
  #                NULL.pass[i] <- ifelse(NULL.pos > 0, 
  #                    sum(phi.norm[,i] >= NES[k])/NULL.pos, 0)
  #                OBS.pos <- sum(obs.phi.norm[,i] >= 0)
  #                OBS.pass[i] <- ifelse(OBS.pos > 0, 
  #                    sum(obs.phi.norm[,i] >= NES[k])/OBS.pos, 0)
  #            }
  #        } 
  #        else 
  #        {
  #            for (i in seq_len(nperm)) 
  #            {
  #                NULL.neg <- sum(phi.norm[,i] < 0)
  #                NULL.pass[i] <- ifelse(NULL.neg > 0, 
  #                    sum(phi.norm[,i] <= NES[k])/NULL.neg, 0)
  #                OBS.neg <- sum(obs.phi.norm[,i] < 0)
  #                OBS.pass[i] <- ifelse(OBS.neg > 0, 
  #                    sum(obs.phi.norm[,i] <= NES[k])/OBS.neg, 0)
  #            }
  #        }
  #        glob.p.vals[k] <- sum(NULL.pass >= mean(OBS.pass))/nperm
  #    }
  #    glob.p.vals.sorted <- glob.p.vals[Orig.index]
  
  # Produce results report
  Obs.ES <- signif(Obs.ES, digits=5)
  Obs.ES.norm <- signif(Obs.ES.norm, digits=5)
  p.vals <- signif(p.vals, digits=4)
  #    signal.strength <- signif(signal.strength, digits=3)
  #    tag.frac <- signif(tag.frac, digits=3)
  #    gene.frac <- signif(gene.frac, digits=3)
  #    FDR.mean.sorted <- signif(FDR.mean.sorted, digits=5)
  #    FDR.median.sorted <-  signif(FDR.median.sorted, digits=5)
  #    glob.p.vals.sorted <- signif(glob.p.vals.sorted, digits=5)
  
  report <- DataFrame(gs.names, size.G, Obs.ES, Obs.ES.norm, p.vals)
  #       p.vals[,1], FDR.mean.sorted, p.vals[,2], tag.frac, 
  #         gene.frac, signal.strength, FDR.median.sorted, glob.p.vals.sorted))
  colnames(report) <- c("GS", "SIZE", "ES", "NES", configEBrowser("PVAL.COL"))#, 
  rownames(report) <- NULL
  return(report)
  #         "FDR q-val", "FWER p-val", "Tag \\%", "Gene \\%", "Signal", 
  #         "FDR (median)", "glob.p.val")
  #    report2 <- report[order(Obs.ES.norm, decreasing=TRUE),]
  #    report3 <- report[order(Obs.ES.norm, decreasing=FALSE),]
  #    phen1.rows <- sum(Obs.ES.norm >= 0)
  #    phen2.rows <- length(Obs.ES.norm) - phen1.rows
  #    report.phen1 <- report2[seq_len(phen1.rows),]
  #    report.phen2 <- report3[seq_len(phen2.rows),]
  
  #    if (output.directory != "")  
  #    {
  #        if (phen1.rows > 0) 
  #        {
  #            filename <- paste(output.directory, doc.string, 
  #                ".SUMMARY.RESULTS.REPORT.", phen1,".txt", sep="", collapse="")
  #            write.table(report.phen1, 
  #                file = filename, quote=FALSE, row.names=FALSE, sep = "\t")
  #        }
  #        if (phen2.rows > 0) {
  #            filename <- paste(output.directory, doc.string, 
  #                ".SUMMARY.RESULTS.REPORT.", phen2,".txt", sep="", collapse="")
  #            write.table(report.phen2, 
  #                file = filename, quote=FALSE, row.names=FALSE, sep = "\t")
  #        }
  #    }
  #
  #    return(list(report1 = report.phen1, report2 = report.phen2))
  
}  # end of definition of GSEA.analysis

# Auxiliary functions and definitions 

.GSEA.GeneRanking <- function(A, class.labels, gene.labels, nperm) 
{ 
  
  # This function ranks the genes according to the signal to noise ratio for the 
  # actual phenotype and also random permutations and bootstrap subsamples of both
  # the observed and random phenotypes. It uses matrix operations to implement the
  # signal to noise calculation in stages and achieves fast execution speed. It 
  # supports two types of permutations: random (unbalanced) and balanced. It also
  # supports subsampling and bootstrap by using masking and multiple-count
  # variables.  When "fraction" is set to 1 (default) the there is no subsampling
  # or boostrapping and the matrix of observed signal to noise ratios will have 
  # the same value for all permutations. This is wasteful but allows to support 
  # all the multiple options with the same code. Notice that the second matrix for
  # the null distribution will still have the values for the random permutations 
  # (null distribution). This mode (fraction = 1.0) is the defaults, the 
  # recommended one and the one used in the examples.It is also the one that has 
  # be tested more thoroughly. The resampling and boostrapping options are 
  # intersting to obtain smooth estimates of the observed distribution but its is
  # left for the expert user who may want to perform some sanity checks before 
  # trusting the code.
  #
  # Inputs:
  #   A: Matrix of gene expression values (rows are genes, columns are samples) 
  #   class.labels: Phenotype of class disticntion of interest. 
  #       A vector of binary labels having first the 1's and then the 0's 
  #   gene.labels: gene labels. Vector of probe ids or 
  #       accession numbers for the rows of the expression matrix 
  #   nperm: Number of random permutations/bootstraps to perform 
  #
  # Outputs:
  #   s2n.matrix: Matrix with random permuted or bootstraps signal to noise ratios
  #       (rows are genes, columns are permutations or bootstrap subsamplings)
  #   obs.s2n.matrix: Matrix with observed signal to noise ratios 
  #       (rows are genes, columns are boostraps subsamplings. 
  #           If fraction is set to 1.0 then all the columns have the same values
  #   order.matrix: Matrix with the orderings that will 
  #       sort the columns of the obs.s2n.matrix in decreasing s2n order
  #   obs.order.matrix: Matrix with the orderings that will 
  #       sort the columns of the s2n.matrix in decreasing s2n order
  #
  A <- A + 0.00000001
  
  N <- nrow(A)
  Ns <- ncol(A)
  
  subset.mask <- reshuffled.class.labels1 <- reshuffled.class.labels2 <- 
    class.labels1 <- class.labels2 <- matrix(0, nrow=Ns, ncol=nperm)
  order.matrix <- obs.order.matrix <- s2n.matrix <- 
    obs.s2n.matrix <- matrix(0, nrow = N, ncol = nperm)
  #obs.gene.labels <- obs.gene.descs <- 
  #    obs.gene.symbols <- vector(length = N, mode="character")
  M1 <- M2 <- S1 <- S2 <- matrix(0, nrow = N, ncol = nperm)
  C <- split(class.labels, class.labels)
  class1.size <- length(C[[2]])
  class2.size <- length(C[[1]])
  class1.index <- seq_len(class1.size)
  class2.index <- (class1.size + 1):(class1.size + class2.size)
  
  for (r in seq_len(nperm)) 
  {
    class1.subset <- sample(class1.index, size = ceiling(class1.size))
    class2.subset <- sample(class2.index, size = ceiling(class2.size))
    subset.class1 <- as.integer(class1.index %in% class1.subset)
    subset.class2 <- as.integer(class2.index %in% class2.subset)
    
    subset.mask[, r] <- as.numeric(c(subset.class1, subset.class2))
    fraction.class1 <- class1.size/Ns
    fraction.class2 <- class2.size/Ns
    # random (unbalanced) permutation
    full.subset <- c(class1.subset, class2.subset)
    label1.subset <- sample(full.subset, size = Ns * fraction.class1)
    for (i in seq_len(Ns)) 
    {
      m1 <- sum(!is.na(match(label1.subset, i)))
      m2 <- sum(!is.na(match(full.subset, i)))
      reshuffled.class.labels1[i, r] <- m1
      reshuffled.class.labels2[i, r] <- m2 - m1
      if (i <= class1.size) 
      {
        class.labels1[i, r] <- m2
        class.labels2[i, r] <- 0
      } 
      else 
      {
        class.labels1[i, r] <- 0
        class.labels2[i, r] <- m2
      }
    }
  }
  
  # compute S2N for the random permutation matrix
  P <- reshuffled.class.labels1 * subset.mask
  n1 <- sum(P[,1])         
  M1 <- A %*% P
  M1 <- M1/n1      
  A2 <- A*A        
  S1 <- A2 %*% P   
  S1 <- S1/n1 - M1*M1    
  S1 <- sqrt(abs((n1/(n1-1)) * S1))   
  P <- reshuffled.class.labels2 * subset.mask
  n2 <- sum(P[,1])           
  M2 <- A %*% P           
  M2 <- M2/n2          
  A2 <- A*A           
  S2 <- A2 %*% P      
  S2 <- S2/n2 - M2*M2 
  S2 <- sqrt(abs((n2/(n2-1)) * S2))
  rm(P)
  rm(A2)
  
  # small sigma "fix" as used in GeneCluster
  S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
  S2 <- ifelse(S2 == 0, 0.2, S2)
  S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
  S1 <- ifelse(S1 == 0, 0.2, S1)
  M1 <- M1 - M2
  rm(M2)
  S1 <- S1 + S2
  rm(S2)
  s2n.matrix <- abs(M1/S1)
  order.matrix <- apply(s2n.matrix, 2, order, decreasing=TRUE)
  
  # compute S2N for the "observed" permutation matrix
  P <- class.labels1 * subset.mask
  n1 <- sum(P[,1])         
  M1 <- A %*% P
  M1 <- M1/n1      
  A2 <- A*A        
  S1 <- A2 %*% P   
  S1 <- S1/n1 - M1*M1    
  S1 <- sqrt(abs((n1/(n1-1)) * S1))   
  P <- class.labels2 * subset.mask
  n2 <- sum(P[,1])           
  M2 <- A %*% P           
  M2 <- M2/n2          
  A2 <- A*A           
  S2 <- A2 %*% P      
  S2 <- S2/n2 - M2*M2 
  S2 <- sqrt(abs((n2/(n2-1)) * S2))
  rm(P)
  rm(A2)
  
  # small sigma "fix" as used in GeneCluster
  S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
  S2 <- ifelse(S2 == 0, 0.2, S2)
  S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
  S1 <- ifelse(S1 == 0, 0.2, S1)
  M1 <- M1 - M2
  rm(M2)
  S1 <- S1 + S2
  rm(S2)
  obs.s2n.matrix <- abs(M1/S1)
  obs.order.matrix <- apply(obs.s2n.matrix, 2, order, decreasing=TRUE)
  
  return(list(s2n.matrix = s2n.matrix, 
              obs.s2n.matrix = obs.s2n.matrix, 
              order.matrix = order.matrix,
              obs.order.matrix = obs.order.matrix))
}

.GSEA.EnrichmentScore <- function(
  gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) 
{  
  #
  # Computes the weighted GSEA score of gene.set in gene.list. 
  # The weighted score type is the exponent of the correlation 
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). 
  # When the score type is 1 or 2 it is necessary to input the correlation vector
  # with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list 
  #       (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set (e.g. integers indicating 
  #       the location of those genes in the input dataset) 
  #   weighted.score.type: Type of score: 
  #       weight: 0 (unweighted = Kolmogorov-Smirnov), 
  #       1 (weighted), and 2 (over-weighted)  
  #   correl.vector: A vector with the coorelations (e.g. signal to noise scores)
  #       corresponding to the genes in the gene list 
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1) 
  #   arg.ES: Location in gene.list where the peak 
  #       running enrichment occurs (peak of the "mountain") 
  #   RES: Numerical vector containing the running 
  #       enrichment score for all locations in the gene list 
  #   tag.indicator: Binary vector indicating 
  #       the location of the gene sets (1's) in the gene list 
  
  # notice that the sign is 0 (no tag) or 1 (tag) 
  tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))
  no.tag.indicator <- 1 - tag.indicator 
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  if (weighted.score.type == 0) correl.vector <- rep(1, N)
  alpha <- weighted.score.type
  correl.vector <- abs(correl.vector**alpha)
  sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
  norm.tag    <- 1.0 / sum.correl.tag
  norm.no.tag <- 1.0 / Nm
  RES <- cumsum(tag.indicator * correl.vector * 
                  norm.tag - no.tag.indicator * norm.no.tag)      
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > -min.ES) 
  {
    ES <- signif(max.ES, digits = 5)
    arg.ES <- which.max(RES)
  } 
  else 
  {
    ES <- signif(min.ES, digits=5)
    arg.ES <- which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
}


.GSEA.EnrichmentScore2 <- function(
  gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) 
{  
  #
  # Computes the weighted GSEA score of gene.set in gene.list. It is the same 
  # calculation as in GSEA.EnrichmentScore but faster (x8) without producing the 
  # RES, arg.RES and tag.indicator outputs.
  # This call is intended to be used to asses the enrichment of random 
  # permutations rather than the observed one.
  # The weighted score type is the exponent of the correlation 
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). 
  # When the score type is 1 or 2 it is necessary to input the correlation vector
  # with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list 
  #       (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set 
  #       (e.g. integers indicating the location of those genes in the input dataset) 
  #   weighted.score.type: Type of score: 
  #       weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
  #   correl.vector: A vector with the coorelations (e.g. signal to noise scores)
  #       corresponding to the genes in the gene list 
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1) 
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  
  peak.res.vector <- valley.res.vector <- 
    tag.diff.vector <- vector(length=Nh, mode="numeric")
  
  loc.vector <- vector(length=N, mode="numeric")
  loc.vector[gene.list] <- seq_len(N)
  tag.loc.vector <- loc.vector[gene.set]
  tag.loc.vector <- sort(tag.loc.vector, decreasing =FALSE)
  
  if (weighted.score.type == 0) tag.correl.vector <- rep(1, Nh)
  else if (weighted.score.type == 1) 
  {
    tag.correl.vector <- correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } 
  else if (weighted.score.type == 2) 
  {
    tag.correl.vector <- 
      correl.vector[tag.loc.vector] * correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } 
  else 
  {
    tag.correl.vector <- correl.vector[tag.loc.vector] * weighted.score.type
    tag.correl.vector <- abs(tag.correl.vector)
  }
  
  norm.tag <- 1.0/sum(tag.correl.vector)
  tag.correl.vector <- tag.correl.vector * norm.tag
  norm.no.tag <- 1.0/Nm
  tag.diff.vector[1] <- (tag.loc.vector[1] - 1) 
  tag.diff.vector[2:Nh] <- 
    tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
  tag.diff.vector <- tag.diff.vector * norm.no.tag
  peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
  valley.res.vector <- peak.res.vector - tag.correl.vector
  max.ES <- max(peak.res.vector)
  min.ES <- min(valley.res.vector)
  ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
  
  return(list(ES = ES))
}
