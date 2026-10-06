#!/usr/bin/env Rscript

# Local PCA (lostruct) on one chromosome.
#
# Usage:
#   Rscript local_pca.R <indir> <chrom> <outdir> [window_size] [snp|bp] [cpus]
#
#   indir        directory containing <chrom>.bcf (plus its .csi index)
#   chrom        chromosome / scaffold name, also the BCF file prefix
#   outdir       where results are written
#   window_size  default 50000
#   snp|bp       window unit, default bp
#   cpus         default: 1
#

# a single library contains all the functionality we need
library(lostruct)

# grab positional arguments
args   = commandArgs(trailingOnly = TRUE)
nargs  = length(args) 
if (nargs < 3 || nargs > 6){
  stop("Usage: local_pca.R <indir> <chrom> <outdir> [window_size] [snp|bp] [cpus]",
       call. = FALSE)
}

# assume the arguments are in proper order
indir    = args[1]
chrom    = args[2]
outdir   = args[3]
win_size = if (nargs >= 4) as.integer(args[4]) else 50000 # 50 Kbp window default
win_type = if (nargs >= 5) args[5] else "bp"
cpus     = if (nargs == 6) max(as.integer(args[6]), 1) else 1

# check arguments
if (is.na(win_size) || win_size < 1){
  stop("window_size must be a positive integer", call. = FALSE)
}
if (!win_type %in% c("snp", "bp")){
  stop("window type must be 'snp' or 'bp'", call. = FALSE) 
}
if (is.na(cpus)){
  stop("Invalid argument passed for cpus", call. = FALSEå)
}

# ensure the vcf file exist
vcf = file.path(indir, paste0(chrom, ".bcf"))
if (!file.exists(vcf)){
  stop("Could not locate a .bcf file for ", vcf, call. = FALSE)
}

# silently create the output directory if it doesn't exist
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# start computing
sites  = vcf_positions(vcf)
win_fn = vcf_windower(file = vcf, size = win_size, type = win_type, sites = sites)

# compute the pcs
outpcs  = file.path(outdir, paste0(chrom, ".pcs.rds"))
pcs     = eigen_windows(win_fn, k = 2, mc.cores = cpus) # only compute 2 PCs per window
saveRDS(pcs, file = outpcs)                           # save this object

outpcs_dist = file.path(outdir, paste0(chrom, ".pcsdist.rds"))
pcdist      = pc_dist(pcs, npc = 2, mc.cores = cpus) # pairwise distances between windows
saveRDS(pcdist, file = outpcs_dist)
na.inds     = is.na(pcdist[, 1])                     # remove failed windows

# now to find the differentiating windows
outfile                 = file.path(outdir, paste0(chrom, "_windows.tsv"))
win_pos                 = as.data.frame(region(win_fn)())
colnames(win_pos)[1:3]  = c("chr", "start", "end")
win_pos$window          = seq_len(nrow(win_pos))
win_pos$failed          = na.inds
write.table(win_pos, outfile, sep = "\t", quote = FALSE, row.names = FALSE)

# collect the samples to place onto the figure
good    = which(!na.inds)
N_MDS   = 4 
if (length(good) <= N_MDS) {
  message("Only ", length(good), " usable windows on ", chrom, "; skipping MDS")
  quit(save = "no", status = 0)
}
outfile                   = file.path(outdir, paste0(chrom, "_mds_windows.tsv"))
mds                       = cmdscale(pcdist[good, good], eig = TRUE, k = N_MDS)
mds_tab                   = cbind(win_pos[good, c("window", "chr", "start", "end")], mds$points)
colnames(mds_tab)[-(1:4)] = paste0("MDS", seq_len(N_MDS))
write.table(mds_tab, outfile, sep = "\t", quote = FALSE, row.names = FALSE)

# ---- plot ---------------------------------------------------------------
mid_mb = (mds_tab$start + mds_tab$end) / 2 / 1e6
outpng = file.path(outdir, paste0(chrom, "_mds.png"))
png(outpng, width = 1600, height = 2000, res = 150)

par(mfrow = c(3, 1), mar = c(4.5, 4.5, 2, 1))
plot(mid_mb, mds_tab$MDS1, pch = 16, cex = 0.5, xlab = "Position (Mb)", ylab = "MDS1", main = chrom)
plot(mid_mb, mds_tab$MDS2, pch = 16, cex = 0.5, xlab = "Position (Mb)", ylab = "MDS2")
plot(mds_tab$MDS1, mds_tab$MDS2, pch = 16, cex = 0.5, xlab = "MDS1", ylab = "MDS2")
invisible(dev.off())


