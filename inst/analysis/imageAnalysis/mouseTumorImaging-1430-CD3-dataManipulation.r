# Load R packages
library(openxlsx)
library(plyr)
library(ggplot2)
library(gridExtra)
library(spatstat)
library(lsmeans)

fileDir <- system.file("mouse/raw", 
  package="IMvigor210CoreBiologies")

# Read cd3-coordinate and perimeter data
barcode <- c(paste("33K", 5:9, sep=""), paste("33K", LETTERS, sep=""), paste("33L", 0:9, sep=""), paste("33L", LETTERS[1:18], sep=""))
barcode <- barcode[-c(6, 45)]
cd3.coord <- per.coord <- data.frame()
for(i in 1:length(barcode))
{
  # Read ROI and first-pass cell coordinates
  wkbk <- loadWorkbook(file.path(fileDir,
    paste("20170515_MSR10166_cd3_1430_Lx_coord_",barcode[i],".xlsx",sep="")))
  s <- length(sheets(wkbk))
  dat.1L0 = rbind(readWorkbook(wkbk, sheet = 3, colNames = FALSE), readWorkbook(wkbk, sheet = 4, colNames = FALSE))
  dat.1L1 = rbind(readWorkbook(wkbk, sheet = 5, colNames = FALSE), readWorkbook(wkbk, sheet = 6, colNames = FALSE))
  dat.1L2 = rbind(readWorkbook(wkbk, sheet = 7, colNames = FALSE), readWorkbook(wkbk, sheet = 8, colNames = FALSE))
  dat.1per = readWorkbook(wkbk, sheet = 9, colNames = FALSE)

  cd3.coord <- rbind(cd3.coord, rbind(data.frame(Barcode = barcode[i], sub = 1, level = "L0", dat.1L0),
                                      data.frame(Barcode = barcode[i], sub = 1, level = "L1", dat.1L1),
                                      data.frame(Barcode = barcode[i], sub = 1, level = "L2", dat.1L2)))
  per.coord <- rbind(per.coord, rbind(data.frame(Barcode = barcode[i], sub = 1, dat.1per)))

  sec <- 2
  while(s > 9)
  {
    dat.2L0 = rbind(readWorkbook(wkbk, sheet = (sec - 1)*9 + 3, colNames = FALSE), readWorkbook(wkbk, sheet = (sec - 1)*9 + 4, colNames = FALSE))
    dat.2L1 = rbind(readWorkbook(wkbk, sheet = (sec - 1)*9 + 5, colNames = FALSE), readWorkbook(wkbk, sheet = (sec - 1)*9 + 6, colNames = FALSE))
    dat.2L2 = rbind(readWorkbook(wkbk, sheet = (sec - 1)*9 + 7, colNames = FALSE), readWorkbook(wkbk, sheet = (sec - 1)*9 + 8, colNames = FALSE))
    dat.2per = readWorkbook(wkbk, sheet = (sec - 1)*9 + 9, colNames = FALSE)

    cd3.coord <- rbind(cd3.coord, rbind(data.frame(Barcode = barcode[i], sub = sec, level = "L0", dat.2L0),
                                        data.frame(Barcode = barcode[i], sub = sec, level = "L1", dat.2L1),
                                        data.frame(Barcode = barcode[i], sub = sec, level = "L2", dat.2L2)))
    per.coord <- rbind(per.coord, rbind(data.frame(Barcode = barcode[i], sub = sec, dat.2per)))
    s <- s - 9
    sec <- sec + 1
  }

  print(paste(i, barcode[i]))
}

# Read metadata and map treatment groups
wkbk <- loadWorkbook(file.path(fileDir,
  paste("./20170411_MSR10166_cd3_1430_subset_var.xlsx",sep="")))
mice.info <- readWorkbook(wkbk, sheet = 8, colNames = TRUE)
mice.info.small <- mice.info[,c("barcode.text","Antibody","Group","Treatment")]
mice.info.small$Barcode <- sapply(strsplit(mice.info.small$barcode, "-"), function(x) x[1])
mice.info.small <- mice.info.small[which(mice.info.small$Antibody == "CD3*" & !is.na(mice.info.small$Treatment) & !(mice.info.small$Barcode %in% c("33KA","33LD"))),]
cd3.coord <- join(cd3.coord, mice.info.small[,c("Barcode","Group")])
per.coord <- join(per.coord, mice.info.small[,c("Barcode","Group")])

## NOTE - Coordinate system is as used in MATLAB, i.e. YX instead of XY for perimeter
# Transform coordinate system for perimeter
temp <- per.coord$X2
per.coord$X2 <- per.coord$X1
per.coord$X1 <- temp
rm(temp)

# Remove (0, 0)
cd3.coord <- cd3.coord[-which(cd3.coord$X1 == 0 & cd3.coord$X2 == 0),]

# Check for and fix holes
temp <- per.coord
for(i in 1:length(barcode))
{
  Dat.per <- per.coord[which(per.coord$Barcode == barcode[i]),]
  Dat.cd3 <- cd3.coord[which(cd3.coord$Barcode == barcode[i]),]
  num.newSubs <- 0
  for(k in 1:length(unique(Dat.per$sub)))
  {
    dat.per <- Dat.per[which(Dat.per$sub == unique(Dat.per$sub)[k]),]
    dis <- sqrt(abs(diff(dat.per$X1))^2 + abs(diff(dat.per$X2))^2)
    dis.1 <- sqrt(abs(dat.per$X1[-1] - dat.per$X1[1])^2 + abs(dat.per$X2[-1] - dat.per$X2[1])^2)
    index.gaps <- which(dis > 100)
    print(paste(barcode[i],unique(Dat.per$sub)[k],length(index.gaps)))
    if(length(index.gaps) > 0)
    {
      stopifnot(!(any(c(0, dim(dat.per)[1]) %in% index.gaps)))
      index.gaps <- c(0, index.gaps, dim(dat.per)[1])
      sub.corrected <- data.frame()
      for(j in 1:(length(index.gaps) - 1))
      {
        ind <- (index.gaps[j] + 1):(index.gaps[j+1])
        sub.corrected <- rbind(sub.corrected, data.frame(sub = sum(num.newSubs) + j, X1 = dat.per$X1[ind], X2 = dat.per$X2[ind]))
      }
      stopifnot(dim(dat.per)[1] == dim(sub.corrected)[1])
      dat.per[,c("sub","X1","X2")] <- sub.corrected
      num.newSubs <- c(num.newSubs, length(index.gaps) - 1)
    } else {
      dat.per$sub <- sum(num.newSubs) + 1
      num.newSubs <- c(num.newSubs, 1)
    }
    per.coord[which(temp$Barcode == barcode[i] & temp$sub == unique(Dat.per$sub)[k]),] <- dat.per
  }
}
rm(temp)

# Create a hyperframe of the data
i = 1
dat.cd3 <- cd3.coord[which(cd3.coord$Barcode == barcode[i]),]
dat.per <- per.coord[which(per.coord$Barcode == barcode[i]),]
s <- length(unique(dat.per$sub))
w <- list(); length(w) <- s
for(j in 1:s)
{
  dat.per.s <- dat.per[which(dat.per$sub == j),]
  w[[j]] <- list(x=dat.per.s$X1, y=dat.per.s$X2)
}
W <- owin(poly = w)
p <- ppp(dat.cd3$X1, dat.cd3$X2, window=W)
dat57 <- hyperframe(Barcode=mice.info.small$Barcode[which(mice.info.small$Barcode == barcode[i])],
                    Group=mice.info.small$Group[which(mice.info.small$Barcode == barcode[i])],
                    window = W, y = p)
for(i in 2:length(barcode))
{
  dat.cd3 <- cd3.coord[which(cd3.coord$Barcode == barcode[i]),]
  dat.per <- per.coord[which(per.coord$Barcode == barcode[i]),]
  s <- length(unique(dat.per$sub))
  w <- list(); length(w) <- s
  for(j in 1:s)
  {
    dat.per.s <- dat.per[which(dat.per$sub == j),]
    w[[j]] <- list(x=dat.per.s$X1, y=dat.per.s$X2)
  }
  W <- owin(poly = w)
  p <- ppp(dat.cd3$X1, dat.cd3$X2, window=W)
  dat57 <- rbind(dat57,
                 hyperframe(Barcode=mice.info.small$Barcode[which(mice.info.small$Barcode == barcode[i])],
                      Group=mice.info.small$Group[which(mice.info.small$Barcode == barcode[i])],
                      window = W, y = p))
  print(i)
}

## dat57 included in package
#save(dat57, file = paste("dat57.RData",sep=""))
