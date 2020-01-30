packages = c("eeptools", "lme4", "car", "beeswarm", "MASS", "survival", "matrixStats", "zoo", "data.table", "lmerTest", "caret", "lmtest")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) install.packages(setdiff(packages, rownames(installed.packages())))  
for (i in packages) library(i, character.only = T)



pvalsig = function (x) ifelse(x < 0.001, "***", ifelse(x < 0.01, "**", ifelse(x < 0.05, "*", round(x, 4))))

datefun = function (x) {
  x = gsub("/", "-", x)
  date1 = sub("-.*", "", x)
  date2 = sub("-.*", "", sub("[[:digit:]]*-", "", x))
  date3 = sub(".*-", "", x)
  if (max(nchar(date1)) == 4) {
    as.Date(x)
  } else if (max(nchar(date3)) == 4) {
    if (max(as.numeric(date2), na.rm = T) > 12 & max(as.numeric(date1), na.rm = T) <= 12) {
      as.Date(paste0(date3, "/", date1, "/", date2))
    } else if (max(as.numeric(date1), na.rm = T) > 12 & max(as.numeric(date2)) <= 12){
      as.Date(paste0(date3, "/", date2, "/", date1))
    } else as.Date(paste0(date3, "/", date2, "/", date1)) # Might mistake days for months
  } else if (max(as.numeric(date1), na.rm = T) > 31 | max(as.numeric(date3), na.rm = T) > as.integer(format(Sys.Date(),"%Y"))-2000) {
    as.Date(paste0(20,x))
  } else if (max(as.numeric(date3), na.rm = T) > 31 | max(as.numeric(date1), na.rm = T) > as.integer(format(Sys.Date(),"%Y"))-2000) {
    if (max(as.numeric(date2), na.rm = T) > 12 & max(as.numeric(date1), na.rm = T) <= 12) {
      as.Date(paste0(20, date3, "/", date1, "/", date2))
    } else if (max(as.numeric(date1), na.rm = T) > 12 & max(as.numeric(date2), na.rm = T) <= 12){
      as.Date(paste0(20, date3, "/", date2, "/", date1))
    } else as.Date(paste0(20, date3, "/", date2, "/", date1)) # Might mistake days for months
  } else if (max(as.numeric(date2), na.rm = T) > 12) {
    as.Date(paste0(20, date3, "/", date1, "/", date2))
  } else as.Date(paste0(20, date3, "/", date2, "/", date1)) # Might mistake between days, months and years
}

boxcoxt = function(y) {
  lambda = BoxCoxTrans(y, na.rm = T)$lambda
  if (round(lambda, 1) == 0) {
    y = log(y)
  } else y = (y^lambda - 1)/lambda
}

stresiduals = function (x) {
  s = sqrt(abs(deviance(x))/df.residual(x))
  hii = hatvalues(x)
  residuals(x)/(s * sqrt(1 - hii))
}

colfun = function (x) {
  if (x == 2) {
    c("black", "red")
  } else if (x == 3) {
    c("black", "red", "brown")
  } else if (x == 4) {
    c("black", "red", "grey60", "brown")
  } else if (x == 5) {
    c("black", "red", "grey60", "brown", "tomato")
  } else if (x == 6) {
    c("black", "red", "grey60", "brown", "grey30", "tomato")
  } else if (x == 7) {
    c("black", "red", "grey60", "brown", "grey30", "tomato", "orange")
  } else if (x == 8) {
    c("black", "red", "grey60", "brown", "grey30", "tomato", "grey75", "orange")
  } else if (x > 8) {
    rainbow(x, start = 0, end = 2/6)
  }
}


if (exists("project")) if (project == "mir29") colfun = function (x) {
  if (x == 2) {
    c("blue", "red")
  } else if (x == 3) {
    c("blue", "red", "green4")
  }
}

plotsigaxes = data.frame(rbind(
  c(2,1,2,0),
  c(3,1,2,0),c(3,1,3,1),c(3,2,3,0),
  c(4,1,2,0),c(4,1,3,1),c(4,2,4,0),c(4,3,4,1),
  c(5,1,2,0),c(5,2,3,0),c(5,3,4,0),c(5,4,5,0),c(5,1,3,1),c(5,3,5,1),c(5,2,4,2),
  c(6,1,2,0),c(6,3,4,0),c(6,5,6,0),c(6,2,4,1),c(6,4,6,1),c(6,2,6,2),
  c(8,1,2,0),c(8,3,4,0),c(8,5,6,0),c(8,7,8,0),c(8,2,4,1),c(8,4,6,1),c(8,6,8,1),c(8,2,6,2),c(8,2,8,3),c(8,4,8,4)
))
# if (exists("study")) if (study == "old.pro") plotsigaxes = data.frame(rbind(c(4,1,2,0),c(4,2,3,0),c(4,3,4,0),c(4,1,4,1),c(4,1,3,2),c(4,2,4,3)))
names(plotsigaxes) = c("Ngroups", "P1", "P2", "Height")
mainplotsigaxes = plotsigaxes
for (n in 4:8) {
  psa = cbind.data.frame(n, t(combn(1:n, 2)))
  names(psa) = c("Ngroups", "P1", "P2")
  psa$Height = -1
  # psa$Height = 0
  # for (i in 2:nrow(psa)) while (any(sapply(seq(psa[i,2]+1, psa[i,3]), function (x) any(x == psa[1:(i-1),][psa[1:(i-1),"Height"] == psa[i,"Height"], 3])))) {
  # psa[i,"Height"] = psa[i,"Height"] + 1
  # }
  # plotsigaxes = rbind.data.frame(plotsigaxes, psa)
  # plotsigaxes = plotsigaxes[!duplicated(plotsigaxes[,1:3]),]
  
  psa2 = rbind.data.frame(plotsigaxes[plotsigaxes$Ngroups == n,], psa)
  psa2 = psa2[!duplicated(psa2[,1:3]),]
  
  for (i in which(psa2$Height < 0)){
    psa2[i,"Height"] = 0
    # while (any(psa2$Height == 0) & any(sapply(seq(psa2[i,2]+1, psa2[i,3]), function (x) any(x <= psa2[1:(i-1),][psa2[1:(i-1),"Height"] == psa2[i,"Height"], 3])))) {
    while (any(psa2$Height == 0) & any(sapply(seq(psa2[i,2]+1, psa2[i,3]), function (x) any(x == unlist(mapply(seq, psa2[index(psa2) < i & psa2[,"Height"] == psa2[i,"Height"], 2] + 1, psa2[index(psa2) < i & psa2[,"Height"] == psa2[i,"Height"], 3])))))) {
      psa2[i,"Height"] = psa2[i,"Height"] + 1
    }
    plotsigaxes = rbind.data.frame(plotsigaxes, psa2[i,])
  }
}
allplotsigaxes = plotsigaxes[!duplicated(plotsigaxes[,1:3]),]
