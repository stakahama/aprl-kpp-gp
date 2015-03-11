#!/usr/bin/Rscript

###_* inputs ====================

root <- "apinene"
inputrun <- "run_01b"
outputrun <- "run_01b"
type <- c("FB","GAS")[1]

###_* functions ====================

###_* other values ====================

a0file <- "molefrac_init.txt"
cfactor <- 10^9       # atm to ppb

###_* data import ====================

###_ . vapor pressures
props <- read.csv(sprintf("%s_props_298.csv",root),row.names="compound")

###_ . compound indices
compounds <- read.csv("compound_indices_table.csv",
                      colClasses=c("compound"="character"),
                      na.string="",
                      row.names="compound")

if(type=="FB") {

  init.allspec <- 1e-5  # ppb
  
  Strip <- function(x,char="[ ;]*")
    gsub(sprintf("^%s|%s$",char,char),"",x)

  Convert <- function(x) {
    x <- Strip(x)
    setNames(x[2],x[1])
  }

###_ . [FB] input concentrations  
  inp <- sapply(strsplit(readLines(file.path(inputrun,"cgas_init.def")),
                         "[ ]?=[ ]?"),Convert)

###_ . [FB] create concentration vector  
  conc <- setNames(rep(init.allspec,nrow(props)),row.names(props))
  inp <- inp[intersect(names(inp),names(conc))]
  conc[names(inp)] <- as.numeric(inp)

###_ . [FB] calculate a0 (mole fractions)
  a0 <- conc/(props$p0*cfactor)
  a0 <- a0/sum(a0)

} else if(type=="GAS") {

  Partition <- function(p,p0) {
    ## p, p0: in units of ppb
    excess <- sapply(p-p0,max,0)
    excess/sum(excess)
  }

###_ . [GAS] input concentrations    
  inp <- unlist(tail(read.csv(file.path(inputrun,"gas",sprintf("%s_formatted.csv",root)),
                              row.names="TIME"),1))

###_ . [GAS] create concentration vector
  conc <- inp[row.names(props)]

###_ . [GAS] calculate a0 (mole fractions)
  a0 <- Partition(conc,props$p0*cfactor)
  
}

###_* export ====================

## create data frame
out <- data.frame(index=compounds[row.names(props),"index"],
                  molefrac=formatC(a0[row.names(props)],digits=5,format="g"))
out <- out[order(out$index),]

## write to file
outfile <- file.path(outputrun,a0file)
cat("#FBinit\n",file=outfile)
write.table(out,outfile,append=TRUE,sep="\t",quote=FALSE,
            col.names=FALSE,row.names=FALSE)
