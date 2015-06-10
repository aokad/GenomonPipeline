#!/usr/local/package/r/3.0.1_icc_mkl/bin/R
#
# Load library
library( "cummeRbund" )

#
# Get arguments
args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]
input_dir   <- args[2]
output_dir  <- args[3]

if (FALSE) {
    sample_name <- 'd_100_1'
    input_dir <- '/home/eigos/Data/GenomonProj/genomon/results/20150413/20150409_small_rna/cufflinks/d_100_1'
    output_dir <- '/home/eigos/Data/GenomonProj/genomon/results/20150413/20150409_small_rna/cummeRbund/d_100_1'
}

#
# Read data
cuff<-readCufflinks( input_dir )

png(filename=paste( output_dir, paste( sample_name, "disp.png", sep="_"), sep="/" ) )
disp<-dispersionPlot(genes(cuff))
disp
dev.off()

#    Error in exists(name, envir = env, mode = mode) :argument "env" is missing, with no default
#png(filename=paste( output_dir, paste( sample_name, "dens.png", sep="_"), sep="/" ) )
#dens<-csDensity(genes(cuff))
#dens
#dev.off()

png(filename=paste( output_dir, paste( sample_name, "box.png", sep="_"), sep="/" ) )
b<-csBoxplot(genes(cuff))
b
dev.off()

png(filename=paste( output_dir, paste( sample_name, "sm.png", sep="_"), sep="/" ) )
sm<-csScatterMatrix(genes(cuff))
sm
dev.off()

png(filename=paste( output_dir, paste( sample_name, "dend.png", sep="_"), sep="/" ) )
dend<-csDendro(genes(cuff))
dend
dev.off()

png(filename=paste( output_dir, paste( sample_name, "volcano.png", sep="_"), sep="/" ) )
v<-csVolcanoMatrix(genes(cuff))
v
dev.off()
