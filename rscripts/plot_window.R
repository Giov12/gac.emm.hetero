#!/usr/bin/env Rscript

# plot a PCA window with the option to include the flanking windows
#
# Usage:
#   Rscript plot_window.R <chrom_outdir> <chrom> <position> <samples.txt> <groups.tsv> <out.png> [flank]
#
#   chrom_outdir  folder holding <chrom>.pcs.rds and <chrom>_windows.tsv
#   chrom         chromosome / scaffold name
#   position      bp position
#   samples.txt   one sample name per line, in BCF order
#   popmap.tsv    no header, two tab-separated columns: sample<TAB>pop
#   out.png       output figure
#   flank         neighbouring windows to show on each side (default 0)

usage  = "Usage: plot_window.R <chrom_outdir> <chrom> <position> <samples.txt> <groups.tsv> <out.png> [flank]"
args   = commandArgs(trailingOnly = TRUE)
n_args = length(args)

if (n_args < 6 || n_args > 7){
  stop(usage, call. = FALSE)
}

res_dir      = args[1]
chrom        = args[2]
pos          = as.numeric(args[3])
samples_file = args[4]
groups_file  = args[5]
outpng       = args[6]
flank        = if (n_args == 7) max(as.integer(args[7]), 0) else 0

if (is.na(pos)){
  stop("position must be a number", call. = FALSE)
}
if (is.na(flank)){ 
  stop("Invalid number for flanking option", call. = FALSE)
}

pcs_file = file.path(res_dir, paste0(chrom, ".pcs.rds"))
win_file = file.path(res_dir, paste0(chrom, "_windows.tsv"))

for (f in c(pcs_file, win_file, samples_file, groups_file)){
  if (!file.exists(f)){
    stop("Cannot find ", f, call. = FALSE)
  }
}

pcs     = readRDS(pcs_file)
wins    = read.delim(win_file, stringsAsFactors = FALSE)
samples = readLines(samples_file)
samples = samples[nzchar(samples)]
groups  = read.delim(groups_file, header = FALSE, col.names = c("sample", "group"),
                     stringsAsFactors = FALSE)

if (nrow(pcs) != nrow(wins)){
  stop("pcs has ", nrow(pcs), " rows but windows.tsv has ", nrow(wins), call. = FALSE)
}

# one style per group: color by population, shape by type
style = list(
  Paulina     = list(col = "#128a5e", pch = 16),
  Barabara    = list(col = "#1c5cab", pch = 16),
  BearPaw_ENA = list(col = "#b84e22", pch = 17),
  BearPaw_TNP = list(col = "#b84e22", pch = 15)
)

grp = groups$group[match(samples, groups$sample)]
if (anyNA(grp)){
  stop("These samples have no group: ", paste(samples[is.na(grp)], collapse = ", "), call. = FALSE)
}

bad = setdiff(unique(grp), names(style))
if (length(bad) > 0){
  stop("Unknown group(s): ", paste(bad, collapse = ", "), call. = FALSE)
}

cols = sapply(grp, function(g) style[[g]]$col)
pchs = sapply(grp, function(g) style[[g]]$pch)
n    = length(samples)

get_pcs = function(i){
  # eigenvectors for window i: PC1 for every sample & PC2 for every sample
  nm = colnames(pcs)
  if (!is.null(nm) && any(grepl("^PC_1", nm))){
    pc1 = as.numeric(pcs[i, grep("^PC_1", nm)])
    pc2 = as.numeric(pcs[i, grep("^PC_2", nm)])
  } else {
    vec = as.numeric(pcs[i, 4:ncol(pcs)])  # skip total, lam_1, lam_2 (k = 2)
    pc1 = vec[1:n]
    pc2 = vec[(n + 1):(2 * n)]
  }
  if (length(pc1) != n || length(pc2) != n){
    stop("Found ", length(pc1), " PC1 values but ", n, " samples. dim(pcs) = ",
         paste(dim(pcs), collapse = " x "), "; first columns: ",
         paste(head(colnames(pcs), 6), collapse = ", "), call. = FALSE)
  }
  list(pc1 = pc1, pc2 = pc2)
}

# move Barabara to the right to orientate all samples to its respective axis
orient = function(pc1){
  if (mean(pc1[grp == "Barabara"]) < mean(pc1[grp == "Paulina"])) -pc1 else pc1
}

# which windows to draw based on overlap
hit = which(wins$start <= pos & pos <= wins$end)

if (length(hit) == 0){
  stop("No window contains ", chrom, ":", pos, call. = FALSE)
}

center = hit[1] # get the center window & N flanking windows
idx    = seq(max(1, center - flank), min(nrow(wins), center + flank))

fmt = function(x){
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

nindices = length(idx)
png(outpng, width = 3 * nindices + 1.5, height = 3.5, units = "in", res = 300)
layout(matrix(1:(nindices + 1), nrow = 1), widths = c(rep(3, nindices), 1.5))
par(mar = c(4.5, 4.5, 5, 1))

for (i in idx){
  ttl = paste0(chrom, "\n", fmt(wins$start[i]), "-", fmt(wins$end[i]))
  
  if (i == center){
    ttl = paste0(ttl, "\n(contains position)")
  }

  p = get_pcs(i)
  if (isTRUE(wins$failed[i]) || all(is.na(p$pc1))){
    plot.new()
    title(main = ttl, cex.main = 0.9)
    text(0.5, 0.5, "PCA failed\n(too few SNPs)")
    next
  }

  plot(orient(p$pc1), p$pc2, col = cols, pch = pchs, cex = 0.8, xlab = "PC1", ylab = "PC2", main = ttl, cex.main = 0.9)
}
# add the legend to the extra margin space
# last panel holds only the legend
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = names(style),
       col = sapply(style, function(s) s$col),
       pch = sapply(style, function(s) s$pch),
       cex = 0.9, bty = "n")

invisible(dev.off())

message("Saved ", outpng)
