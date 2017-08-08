# Load Phenoflow
library("Phenoflow")

# Set plot settings
my.settings <- list(
  strip.background=list(col="transparent"),
  strip.border=list(col="transparent", cex=5),
  gate=list(col="black", fill="lightblue", alpha=0.2,border=NA,lwd=3),
  panel.background=list(col="lightgray"),
  background=list(col="white"))

### Create a PolygonGate for extracting the single-cell information
### Inland data
sqrcut1 <- matrix(c(8.5,8.5,16,16,
                    3,8.05,19,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate_inland <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")
### Lake Michigan/Muskegon data
sqrcut1 <- matrix(c(8.5,8.5,16,16,
                    3,8,16,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate_MI <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

### HNA and LNA gates
sqrcut1 <- matrix(c(asinh(11000),asinh(11000),16,16,
                    3,9.75,16,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
rGate_HNA <- polygonGate(.gate=sqrcut1, filterId = "HNA")
sqrcut1 <- matrix(c(8.5,8.5,asinh(11000),asinh(11000)
                    ,3, 8,9.75,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
rGate_LNA <- polygonGate(.gate=sqrcut1, filterId = "LNA")
filters <- filters(list(rGate_LNA, rGate_HNA))
flist <- list(filters)

# Parameters of interest
param=c("FL1-H", "FL3-H","SSC-H","FSC-H")

# Go through directorires in ./data and analyze folders with .fcs files in them.
# The data can be found on flowRepository:
# Lake Michigan 2013/2015 + Muskegon Lake 2014/2015: FR-FCM-ZYZN
# Inland Lakes: FR-FCM-ZY9J

teller <- 1
for(dir in list.dirs("./data")){
  files <- list.files(dir, pattern = "*.fcs")
  if(length(files) > 1){
    flowData <- read.flowSet(path = dir, 
                             transformation = FALSE, pattern=".fcs")
    flowData <- flowData[, param]
    flowData <- transform(flowData,`FL1-H`=asinh(`FL1-H`), 
                                      `SSC-H`=asinh(`SSC-H`), 
                                      `FL3-H`=asinh(`FL3-H`), 
                                      `FSC-H`=asinh(`FSC-H`))
    ### Denoise data
    if(dir=="./data/Inland_Lakes"){
      ###  Gating quality check for inland lakes
      names(flist) <- flowCore::sampleNames(flowData[1])
      png(file = paste(dir, "_filtered.png", sep = ""), width = 7, height = 6, res = 500, 
          units = "in", pointsize = 10)
      print(xyplot(`FL3-H`~`FL1-H`, data=flowData[1],
                   filter=flist,
                   xbins=200,nbin=128, par.strip.text=list(col="black", font=3,cex=1.85), 
                   smooth=FALSE, xlim=c(6,16),ylim=c(0,16),xlab=list(label="Green fluorescence intensity (FL1-H)",cex=2),
                   ylab=list(label="Red fluorescence intensity (FL3-H)",cex=2),
                   par.settings=my.settings,
                   scales=list(x=list(at=seq(from=6, to=16, by=2),cex=1),
                               y=list(at=seq(from=0, to=16, by=4),cex=1)),
                   strip.text = "NULL",
                   margin=TRUE,
                   binTrans="log"
      )
      )
      dev.off()
      flowData <- Subset(flowData, polyGate_inland)
      TC <- filter(flowData, polyGate_inland)
      TotalCount <- summary(TC);TotalCount <- toTable(TotalCount)
      HC <- filter(flowData, rGate_HNA) 
      HNACount <- summary(HC);HNACount <- toTable(HNACount)
      LC <- filter(flowData, rGate_LNA) 
      LNACount <- summary(LC);LNACount <- toTable(LNACount)
    } else {
      ###  Gating quality check for muskegon/michigan
      names(flist) <- flowCore::sampleNames(flowData[1])
      png(file = paste(dir, "_filtered.png", sep = ""), width = 7, height = 6, res = 500, 
          units = "in", pointsize = 10)
      print(xyplot(`FL3-H`~`FL1-H`, data=flowData[1],
                   filter=flist,
                   xbins=200,nbin=128, par.strip.text=list(col="black", font=3,cex=1.85), 
                   smooth=FALSE, xlim=c(6,16),ylim=c(0,16),xlab=list(label="Green fluorescence intensity (FL1-H)",cex=2),
                   ylab=list(label="Red fluorescence intensity (FL3-H)",cex=2),
                   par.settings=my.settings,
                   scales=list(x=list(at=seq(from=6, to=16, by=2),cex=1),
                               y=list(at=seq(from=0, to=16, by=4),cex=1)),
                   margin=TRUE,
                   binTrans="log"
        )
      )
      dev.off()
      flowData <- Subset(flowData, polyGate_MI)
      TC <- filter(flowData, polyGate_MI)
      TotalCount <- summary(TC);TotalCount <-  flowCore::toTable(TotalCount)
      HC <- filter(flowData, rGate_HNA) 
      HNACount <- summary(HC);HNACount <- toTable(HNACount)
      LC <- filter(flowData, rGate_LNA) 
      LNACount <- summary(LC);LNACount <- toTable(LNACount)
    }
  
    ### Transform back to initial space before exporting filtered data
    flowData_unfiltered <- transform(flowData,`FL1-H`=sinh(`FL1-H`), 
                                    `SSC-H`=sinh(`SSC-H`), 
                                    `FL3-H`=sinh(`FL3-H`), 
                                    `FSC-H`=sinh(`FSC-H`))
    ### Export fcs files
    if(!dir.exists(paste(dir,"_filtered", sep = ""))) dir.create(paste(dir,"_filtered", sep = ""))
    setwd(paste(dir,"_filtered", sep = ""))
    for(index in 1:length(flowData_unfiltered)){
      write.FCS(flowData_unfiltered[[index]], filename = paste(flowCore::sampleNames(flowData_unfiltered[index]),
                                                      "_filtered", sep = ""))
    }
    setwd("../..")
    
    ### Extract the volumes
    vol.temp <- c()
    for(i in 1:length(flowData)){
      vol.temp[i] <- as.numeric(flowData[[i]]@description$`$VOL`)/1000
    }
    
    Count_results <- data.frame(Samples = sampleNames(flowData), 
                                Total.cells = TotalCount$true/vol.temp,
                                HNA.cells = HNACount$true/vol.temp,
                                LNA.cells = LNACount$true/vol.temp)
    ### Store results
    if(teller == 1){
      df_results <- Count_results
    } else {
      df_results <- rbind(df_results, Count_results)
    }
    teller <- teller + 1
  }
}

# Load file with dilutions of experiments
dilutions_subset <- read.csv2("./data/dilutions.csv")
dilutions_subset <- dilutions_subset[dilutions_subset$dilution == 2, ]
dilutions_subset$Sample <- paste(dilutions_subset$Sample, ".fcs", sep ="")
# First start with giving a dilution of 10 to all samples
# Then change the ones necessary to 2
df_results <- data.frame(df_results, dilution = rep(10, nrow(df_results)))
df_results$dilution[df_results$Samples %in% dilutions_subset$Sample] <- 2

# Multiply all values by their dilution factor

df_results[, c(2,3,4)] <- df_results[, c(2,3,4)] * df_results$dilution

# Export final result table to CSV
write.csv(df_results[, c(1:4)], file = "recalculated_HNA_LNA.csv", quote = FALSE,
           row.names = FALSE)

df_results <- read.csv("recalculated_HNA_LNA.csv", stringsAsFactors = FALSE)

tmp <- strsplit(df_results$Samples, "_rep*")
df_results$Samples <- do.call(rbind, lapply(tmp, rbind))[,1]

df_final_results <- data.frame(
  aggregate(Total.cells~Samples, data = df_results, FUN = mean),
  HNA.cells= aggregate(HNA.cells~Samples, data = df_results, FUN = mean)[,2],
  LNA.cells= aggregate(LNA.cells~Samples, data = df_results, FUN = mean)[,2],
  sd.Total.cells = aggregate(Total.cells~Samples, data = df_results, FUN = sd)[,2],
  sd.HNA.cells = aggregate(HNA.cells~Samples, data = df_results, FUN = sd)[,2],
  sd.LNA.cells = aggregate(LNA.cells~Samples, data = df_results, FUN = sd)[,2]
)

# Export final mean result table to CSV
write.csv(df_final_results, file = "recalculated_mean_HNA_LNA.csv", quote = FALSE,
          row.names = FALSE)

