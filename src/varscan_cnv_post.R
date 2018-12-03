#!/usr/bin/env Rscript

# installation
install.packages('optparse', repos = "http://cran.us.r-project.org")
source("http://bioconductor.org/biocLite.R");
biocLite("DNAcopy");

# load libraries
library ('optparse');
library(DNAcopy);

### command line parameters
option_list <- list(make_option('--in', type = 'character', default = NULL, help = 'varscan input', metavar = 'character'),
                    make_option('--out', type = 'character', default = NULL, help = 'output', metavar = 'character')
                    );

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

# varscan post processing
cn <- read.table(opt$'in' ,header=F)
# raw: chrom   chr_start       chr_stop        num_positions   normal_depth    tumor_depth     log2_ratio      gc_content
# adjusted: chrom   chr_start       chr_stop        num_positions   normal_depth    tumor_depth     adjusted_log_ratio      gc_content      region_call     raw_ratio
CNA.object <-CNA(genomdat = cn[,7], chrom = cn[,1], maploc = cn[,2], data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
segs <- segment(CNA.smoothed, verbose=0, min.width=2)
segs2 = segs$output
write.table(segs2[,2:6], file=opt$out, row.names=F, col.names=F, quote=F, sep="\t")
