#!/usr/bin/Rscript

###_* inputs ====================

## root <- "apinene"
## inputrun <- "run_01b"
## outputrun <- "run_01b"
## type <- c("FB","GAS")[1]

argv <- commandArgs(TRUE)
root <- argv[1]
type <- argv[2]
inputrun <- argv[3]

outputrun <- inputrun

## possible:
if( !type %in% c("purecomponent",
                 "infinitesink",
                 "equalcomponent",
                 "initialequilibrium",
                 "gasphasecomp",
                 "recycledseed",
                 "extrasolventinit") )
  stop("--- not valid 'type' ---")

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

if(type %in% c("initialequilibrium","gasphasecomp","extrasolventinit")) {

  init.allspec <- 1e-5  # ppb
  
###_ . input concentrations

  Strip <- function(x,char="[ ;]*")
    gsub(sprintf("^%s|%s$",char,char),"",x)

  Convert <- function(x) {
    x <- Strip(x)
    setNames(x[2],x[1])
  }
  
  inp <- sapply(strsplit(toupper(readLines(file.path(inputrun,"cgas_init.def"))),
                         "[ ]?=[ ]?"),Convert)

###_ . create concentration vector  
  conc <- setNames(rep(init.allspec,nrow(props)),row.names(props))
  inp <- inp[intersect(names(inp),names(conc))]
  conc[names(inp)] <- as.numeric(inp)

###_ . calculate a0 (mole fractions)

  if(type=="initialequilibrium") {
    a0 <- conc/(props$p0*cfactor)
    a0 <- a0/sum(a0)
  } else if(type=="gasphasecomp") {
    a0 <- conc/sum(conc)
  } else if(type=="extrasolventinit") {
    Navogadro <- 6.02214129e23
    conversion <- c(Navogadro*1e-12,    # microg=>g, m^3=>cm^3
                    1/0.08206*1/298.15) # ppb to microg/m^3 (x MW)
    masstable <- read.table("mcm_apinene_mass.txt",skip=18,
                            col.names=c("compound","SMILES","InChI","molwt"),
                            row.names="compound")
    COA <- scan(file.path(inputrun,"input_partitioning.txt"),0,n=1)
    mean.molwt <- sum(conc*masstable[names(conc),"molwt"])/sum(conc)
    Mconc <- conc*masstable[names(conc),"molwt"]*conversion[2]
    Nconc <- Mconc/masstable[names(conc),"molwt"]*conversion[1]
    Nsolvent <- COA/mean.molwt*conversion[1]-sum(Nconc)
    a0 <- Nconc/(Nsolvent+Nconc)
  }

} else if(type=="recycledseed") {

  Partition <- function(p,p0) {
    ## p, p0: in units of ppb
    excess <- sapply(p-p0,max,0)
    excess/sum(excess)
  }

###_ . input concentrations    
  inp <- unlist(tail(read.csv(file.path(inputrun,"gas",sprintf("%s_formatted.csv",root)),
                              row.names="TIME"),1))

###_ . create concentration vector
  conc <- inp[row.names(props)]

###_ . calculate a0 (mole fractions)
  a0 <- Partition(conc,props$p0*cfactor)
  
} else if(type=="purecomponent") {

  a0 <- setNames(rep(1,nrow(props)),row.names(props))
  
} else if(type=="infinitesink") {

  a0 <- setNames(rep(0,nrow(props)),row.names(props))  

} else if(type=="equalcomponent") {

  a0 <- setNames(rep(1/nrow(props),nrow(props)),row.names(props))
  
}

###_* export ====================

## create data frame
out <- data.frame(index=compounds[row.names(props),"index"],
                  molefrac=formatC(a0[row.names(props)],digits=5,format="g"))
out <- out[order(out$index),]

## write to file
outfile <- file.path(outputrun,a0file)
cat(paste0("#",type,"\n"),file=outfile)
write.table(out,outfile,append=TRUE,sep="\t",quote=FALSE,
            col.names=FALSE,row.names=FALSE)
