#!/usr/bin/Rscript

library(tibble, quietly = T)
library(data.table, quietly = T)
library(readr)
library(dplyr, warn.conflicts = F)
library(stringi, quietly = T)
library(parallel)
library(Rmpfr, warn.conflicts = F, quietly = T)

options(error = function() {traceback(2, max.lines=10); if(!interactive()) quit(save="no", status=1, runLast=T)})

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

cmToInches = function (cm) { # Thanks grDevices...
    cm * 25 / 64
}

# 
# Parsing arguments
# 

args = commandArgs(trailingOnly = T)

indexFiles = args[1] %>% stri_split_fixed(pattern = ",", simplify = T)
coverageFiles = args[2] %>% stri_split_fixed(pattern = ",", simplify = T)
samples = args[3] %>% stri_split_fixed(pattern = ",", simplify = T)
infix = args[4]
condition = args[5]
axisSubsample = as.numeric(args[6])
do.plots = args[7]
cores = as.numeric(args[8])
pLimits = args[9] %>% stri_split_fixed(pattern = ",", simplify = T)

indexFiles = stri_split_fixed(indexFiles, pattern = ",") %>% unlist %>% trimws()
coverageFiles = stri_split_fixed(coverageFiles, pattern = ",") %>% unlist %>% trimws()
samples = stri_split_fixed(samples, pattern = ",") %>% unlist %>% trimws()
if (condition == "no-condition-argument-passed") {
    condition = ""
}

if (F) {
    cat("indexFiles", indexFiles, "\n")
    cat("coverageFiles:", coverageFiles, "\n")
    cat("samples:", samples, "\n")
    cat("infix: ", infix, "\n")
    cat("condition: ", condition, "\n")
    cat("axisSubsample: ", axisSubsample, "\n")
    cat("do.plots: ", do.plots, "\n")
    cat("cores: ", cores, "\n")
    cat("pLimits: ", pLimits, "\n")
}

resultsTable = do.call(bind_rows, 
    lapply(1:length(samples), function (i) {
        sample = samples[[i]]
        coverageFile = coverageFiles[[i]]
        indexFile = indexFiles[[i]]
        
        print(sample)
        
        index = read_tsv(indexFile, col_names = c("chr", "start", "end", "name"))
        index$start = index$start
        
        svs = index$name %>% stri_sort(numeric = T)

        #
        # Since base by base files are large, lets load only one SV at a time
        #
    
        do.call(bind_rows,
                mclapply(svs, function (svName) {
                    # svName = "SV_0_DEL"
                    print(paste(sample, svName))
                    
                    svRelevantRegion = index[index$name == svName, ]
                    
                    svNameShort = stri_replace(svName, regex = "(SV_.*)_.*", replacement = "$1")
                    print(paste0("Reading coverage for ", sample, " : ", svName))
                    fReadCommand = paste0("tabix ", coverageFile, " ", 
                                          svRelevantRegion$chr, ":", 
                                          format(svRelevantRegion$start, scientific = F), 
                                          "-",
                                          format(svRelevantRegion$end, scientific = F),
                                          " | grep -E \"(", svNameShort,"[+_-]|[^\t]*\t0\t)\"" )
                    coverage = tibble(
                        fread(
                            cmd = fReadCommand,
                            col.names = c("chr", "start", "end", "name", "coverage"),
                            showProgress = T)
                    )
                    print(paste0("Finished reading coverage for ", sample, " : ", svName))
                    
                    coverage$start = coverage$start + 1 # BED start coords are 0-based, lets change them to 1-based
    
                    svData = coverage[coverage$name == svName, ]
                    svInData = coverage[stri_detect_fixed(coverage$name, "_in") &
                                        !stri_detect_fixed(coverage$name, "_plot"), ]

                    regionData = coverage[stri_detect_regex(coverage$name, "[+-]") &
                                              !stri_detect_fixed(coverage$name, "_plot"), ]
                    
                    # Let's perform the t-test
                    test = t.test(svInData$coverage, regionData$coverage)
                    t = test$statistic[[1]]
                    df = test$parameter[[1]]
                    
                    probLog = pt(-abs(t), df, log.p = T)
                    pValue = 2 * exp(mpfr(probLog, 100))
                    
                    asterisks = (
                        if (pValue >= mpfr(pLimits[1], 100) | is.na(pValue)) '-'
                        else if (pValue > mpfr(pLimits[2], 100)) '*' 
                        else if (pValue > mpfr(pLimits[3], 100)) '**' 
                        else '***'
                    )
    
                    # plotting: a side-effect
                    if (do.plots) {
                        sv = svData
                        meanCoverageChr = coverage[coverage$name == sv$chr, "coverage"][[1]]
    
                        points = coverage[stri_detect_fixed(coverage$name, pattern = "_plot"), ]
                        points$color = ifelse(rownames(points) %in% rownames(points)[grep("in", points$name)], "limegreen",
                                              ifelse(points$coverage > meanCoverageChr, "steelblue", "maroon"))
                        points$length = points$end - points$start + 1
                        # svPointsRows = stri_detect_fixed(points$name, "in") #This is stored separately for the segment function
                        svPointsRows = grep("in", points$name) #This is stored separately for the segment function
                        moo = par(no.readonly = T)
                        filename = paste0(getwd(), "/test")
                        filename = paste0(getwd(), "/", sample, ".", infix, ".", "coverageAround.", 
                                          paste0("chr"[stri_detect_fixed(sv$chr, "chr", negate = T)], sv$chr), 
                                          "_", sv$start, "-", sv$end)
                        png(paste0(filename, ".png"), width=12.5, height=10, units = "cm", 
                            res=600, pointsize = 8)
                        jpeg(paste0(filename, ".jpg"), width=12.5, height=10, units = "cm", 
                             quality = 100,
                            res=600, pointsize = 8)
                         for (device in dev.list()) {
                            par(oma = c(2,2,0,0), mar = c(5.6, 4, 4.1, 4)) ; {yMax = (if (sv$coverage > mean(regionData$coverage)) sv$coverage * 1.2
                                                                                      else quantile(regionData$coverage, 0.95, names = F))
                            if (yMax < 10) {
                                if (yMax %% 5) yMax = yMax + (5 - yMax %% 5);
                            } else if (yMax %% 10) yMax = yMax + (10 - yMax %% 10)}
                            # Lets initialize an empty plot to put the gray lines
                            barplot(rep(NA, nrow(points)), ylim = c(0, yMax), axes = FALSE,
                                    space = 0,
                                    ylab = "Coverage",
                                    xlab = paste0(stri_trans_totitle(paste0("chr"[stri_detect_fixed(sv$chr, "chr", negate = T)], sv$chr)), " coordinates"),
                                    main = paste0(sample, " - ", paste0(
                                        "chr"[stri_detect_fixed(sv$chr, "chr", negate = T)], sv$chr), ":", sv$start, "-", sv$end, 
                                        ifelse(condition != "", paste0(" (", condition, ")"), "")),
                                    )
                            grid(nx = NA, ny = NULL, col = "gray", lty = 1)
                            mp = barplot(pmin(points$coverage, yMax),
                                         # border = NA,
                                         ylim = c(0, yMax),
                                         space = 0,
                                         col = points$color,
                                         border = points$color,
                                         xaxt = "n", # dont plot x axis (we plot it later)
                                         las = 1,
                                         add = T
                            )
                            
                            # Avg. coverage of chromosome
                            lines(x = c(0, length(mp)), y = rep(meanCoverageChr, 2),
                                  col = "chocolate4", lwd = 1.5, lty = "solid")
                            # Avg. coverage of surroundings
                            lines(x = c(0, length(mp)), y = rep(mean(points[grep("in", points$name, invert = T), ]$coverage), 2),
                                  col = "darkorange2", lwd = 1.5, lty = "dotted")
                            # Avg. coverage of SV
                            segments(x0 = svPointsRows[1] - 1, y0 = sv$coverage, x1 = svPointsRows[length(svPointsRows)], y1 = sv$coverage,
                                     lty = "solid", col = "red")
                            # Asterisks
                            text(x = (svPointsRows[length(svPointsRows)] + (svPointsRows[1] - 1) )/2,
                                 y = sv$coverage + mean(points$coverage) * 0.1, labels = asterisks) # The Mann-Whitneys's pvalue asterisks
                            axisPoints =
                                axis(side = 1, at = c(mp[seq(1, length(mp), length(mp)/axisSubsample)], mp[length(mp)]), labels = F, col = "grey", col.ticks = "grey")
                            text(c(mp[seq(1, length(mp), length(mp)/axisSubsample)], mp[length(mp)]), par("usr")[3] * 4.5,
                                 labels = c(points$start[seq(1, length(mp), length(mp)/axisSubsample)], points$start[nrow(points)]), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.7)
                            legend("bottom", xpd = NA,
                                   inset = c(0, -0.315),
                                   ncol = 2,
                                   legend = c("Cover. chr.", "Cover. region (no SV)", "Cover. SV", "Cover. > chr. mean", "Cover. < chr. mean", "SV"),
                                   col=c("chocolate4", "darkorange2", "red", "steelblue", "maroon", "limegreen"),
                                   lty = c("solid", "dotted", "solid", NA, NA, NA), pch = c(NA, NA, NA, 15, 15, 15), cex = 0.7, box.lty=0)
                            par(moo)
                            dev.off()}
                    }
                    # return value: a row for the tibble
                    pValueFormatted = formatMpfr(pValue)
                    meanCoefficient = svData$coverage / mean(regionData$coverage)
                    concordance = (if (stri_detect_fixed(tolower(svName), pattern = c("del", "dup"), max_count = 1)[1]) {
                        if (stri_detect_fixed(tolower(svName), pattern = "del")) {
                            if (meanCoefficient <= 0.6) {
                                "yes"
                            }
                            else {"no"}
                        }
                        else { # duplications
                            if (meanCoefficient >= 1.4) {
                                "yes"
                            }
                            else {"no"}
                        }
                    } else {"NA"})
                    significant = (if (concordance == "yes") {asterisks} 
                                   else if (concordance == "no") {"-"}
                                   else {"NA"})
                    list(sample = sample, sv = svName, coords = paste0(svData$chr, ":", svData$start, "-", svData$end),
                         meanCoverageSV = svData$coverage,
                         meanCoverageRegion = mean(regionData$coverage),
                         meanCoefficient = meanCoefficient,
                         probLog = probLog,
                         pValue = pValueFormatted,
                         pValueExp = stri_replace_first_regex(pValueFormatted, ".*e(-?.*)", replacement = "$1"),
                         asterisks = asterisks,
                         concordance = concordance,
                         significance = significant
                    )
                },
                mc.cores = cores)
        )
})
)

readr::write_csv(resultsTable, paste0("table.", infix, ".csv"), col_names = T)

