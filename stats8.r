time0 = proc.time()

project = "Heart progerin" # mir29   Heart progerin   Heart lamin  ISO Challenge  PCTX Treatment  Zmpste-Rankl  HGPS Amanda  DBU Alberto

## Import data
if (.Platform$OS.type == "unix") setwd("/Volumes/Victor/") else setwd("S:/LAB_VA/LAB/Victor/")
# setwd("/Users/Victor/Downloads/")
# setwd("S:/LAB_VA/LAB/Alvaro Macias/Results/Esperanza de vida/Challenge de Isoproterenol/Electros/")
# setwd("S:/LAB_VA/LAB/Alvaro Macias/Results/Esperanza de vida/Paclitaxel 0.5mg_kg 3_semana/ECGs/")
# setwd("/Volumes/LAB/Alberto M/Alberto/Ratones endotelio y SMC/ECG DBU/")
if(project == "PCTX Treatment") setwd("/Volumes/LAB_VA/LAB/Alvaro Macias/Results/Tratamientos/Paclitaxel 0.5mg_kg 3_semana/ECGs sucesivos tras PCTX agudo - REAL/")
baseroute = paste0(project, " project/")

source("Methodology/Scripts/utility.r")
# source("S:/LAB_VA/LAB/Victor/Methodology/Scripts/utility.r")
if(project == "PCTX Treatment") source("/Volumes/Victor/Methodology/utility.r")

rconfig = data.frame(fread(paste0(baseroute, "Design/rconfig.txt"), encoding = "Latin-1", stringsAsFactors = F, na.strings = c("N/A", "NA", "#N/A", "")))
rconfig = rconfig[rconfig$Run == 1,-1]
invisible(sapply(names(rconfig), function (x) if (class(rconfig[,x]) == "integer") rconfig[,x] <<- as.logical(rconfig[,x])))

exp = 1
for (exp in 1:nrow(rconfig)) {
  for(opt in 1:ncol(rconfig)) assign(names(rconfig)[opt], rconfig[exp,opt])
  if (torpor_sel != "All") torpor_sel = as.numeric(torpor_sel)
  tryCatch({
    telemetryexp = F
    if (grepl("Telemetry", technique) | technique == "Metabolic.cages") telemetryexp = T
    designdata = data.frame(fread(paste0(baseroute, "Design/", study, ".txt"), encoding = "Latin-1", stringsAsFactors = F, na.strings = c("N/A", "NA", "#N/A", "")))
    designdata$Id = toupper(designdata$Id)
    if (technique == "Telemetry") channel = "CH2"
    rawdata = data.frame(fread(paste0(baseroute, "Raw data/", technique, ifelse(any(c("ECG", "Telemetry") %in% technique), paste0("/analyzed ", channel, "/"), "/"), "rawdata.txt"), encoding = "Latin-1", stringsAsFactors = F, na.strings = c("N/A", "NA", "#N/A", "", "nan")))
    rawdata$Id = toupper(rawdata$Id)
    expdata = merge(designdata, rawdata, by = c("Id"), sort = F)
    # expdata[is.na(expdata$Treatment),"Treatment"] = ""
    expdata = expdata[!is.na(expdata$Date),]
    expdata$DOB = datefun(expdata$DOB)
    if (!is.null(expdata$DOD)) expdata$DOD = datefun(expdata$DOD)
    expdata$Date = datefun(expdata$Date)
    expdata$Start.Date = datefun(expdata$Start.Date)
    expdata$End.Date = datefun(expdata$End.Date)
    expdata = expdata[expdata$Date >= expdata$Start.Date & expdata$Date <= expdata$End.Date,]
    if (!is.null(expdata$DOD)) {
      expdata$Age = as.numeric(age_calc(expdata$DOB,expdata$DOD, "days"))/7
    } else expdata$Age = as.numeric(age_calc(expdata$DOB,expdata$Date, "days"))/7
    expdata$Age0 = expdata$Age
    for (i in unique(expdata$Id)) expdata$Age0[expdata$Id == i] = min(expdata[expdata$Id == i,"Age"], na.rm = T)
    expdata$Day = as.numeric(age_calc(sort(unique(expdata$Date))[1],expdata$Date, "days"))
    for (i in sort(unique(expdata$Condition.Id))) expdata$Condition.Id[expdata$Condition.Id == i] = which(sort(unique(expdata$Condition.Id)) == i)
    # expdata$Condition.Id = sapply(1:length(expdata$Condition.Id), function (x) which(sort(unique(expdata$Condition.Id)) == expdata$Condition.Id[x]))
    followups = data.frame("Follow.up" = 1:length(unique(expdata$Day)), "Day" = sort(unique(expdata$Day)))
    expdata = merge(expdata, followups, by = "Day", all.x = T)
    expdata = expdata[,c(2:ncol(expdata),1)]
    expdata$Follow.up.Id = paste(expdata$Follow.up, expdata$Condition.Id, sep = ".")
    # if (any(names(rawdata) == "Time") & !telemetryexp & technique != "qPCR") {
    # if (any(names(rawdata) == "Time")) {
    if (any(names(rawdata) == "Time") & !telemetryexp) {
      tfollowups = data.frame("TFollow.up" = 1:length(unique(expdata$Time)), "Time" = sort(unique(expdata$Time)))
      expdata = merge(expdata, tfollowups, by = "Time", all.x = T, sort = F)[, union(names(expdata), names(tfollowups))]
      expdata$Time_Id = paste(expdata$TFollow.up, expdata$Condition.Id, sep = ".")
    }
    for (i in sort(unique(expdata$Condition))) expdata[,i] = expdata$Condition == i
    
    
    ## Determine plot types 
    checklevels = function (y) all(sapply(unique(expdata$Condition), function (x) length(unique(expdata[expdata$Condition == x, y]))) > 1)
    
    if (!any(duplicated(expdata$Id))) ind_plots = F
    age_plots = F; if (checklevels("Age")) age_plots = T
    fu_plots = F; if ((followup_boxplots | fl_followup_boxplots) & checklevels("Follow.up")) fu_plots = T
    time_plots = F; if (checklevels(grep("^Time", names(rawdata), value = T))) time_plots = T
    for (i in c("beeboxplots", "followupboxplots", "timeboxplots", "flfollowupboxplots", "fltimeboxplots", "loessplots", "indplots")) assign(i, NULL)
    boxbeeswarm_plots = strsplit(ifelse(boxplots, ifelse(beeswarm_plots, paste("beeswarm", "boxbeeswarm", "boxplot"), "boxplot"), "beeswarm"), split = " ")
    if (telemetryexp)  onlymeanplots = T
    if (onlymeanplots) allmean = "mean" else allmean = c("all", "mean")
    if (boxplots | beeswarm_plots) beeboxplots = as.vector(sapply(allmean, function (x) sapply (boxbeeswarm_plots, function (y) paste(x, y))))
    if (fu_plots) {
      if (followup_boxplots) followupboxplots = as.vector(sapply (boxbeeswarm_plots, function (x) paste("followup", x)))
      if (fl_followup_boxplots) flfollowupboxplots = as.vector(sapply (boxbeeswarm_plots, function (x) paste("first-last followup", x)))
    }
    if (time_plots) {
      if (followup_boxplots) timeboxplots = as.vector(sapply (boxbeeswarm_plots, function (x) paste("time", x)))
      if (fl_followup_boxplots) fltimeboxplots = as.vector(sapply (boxbeeswarm_plots, function (x) paste("first-last time", x)))
    }
    if (loess_plots) {
      if (fu_plots)  loessplots = c(loessplots, "loess vs day", "loess vs followup")
      if (time_plots)  loessplots = c(loessplots, "loess vs time")
      if (time_plots & !telemetryexp)  loessplots = c(loessplots, "loess vs time followup")
      if (age_plots)  loessplots = c(loessplots, "loess vs age")
      if (weight_factor)  loessplots = c(loessplots, "loess vs weight")
      if (temp_factor)  loessplots = c(loessplots, "loess vs temperature")
      if (hr_factor)  loessplots = c(loessplots, "loess vs heart rate")
      if (rr_factor)  loessplots = c(loessplots, "loess vs RR")
    }
    if (ind_plots) {
      if (fu_plots)  indplots = c(indplots, "ind vs day")
      if (time_plots) indplots = c(indplots, "ind vs time")
      if (age_plots)  indplots = c(indplots, "ind vs age")
      if (temp_factor)  indplots = c(indplots, "ind vs temperature")
      if (hr_factor)  indplots = c(indplots, "ind vs heart rate")
      if (rr_factor)  indplots = c(indplots, "ind vs RR")
    }
    
    plot_types = c(ifelse(boxplots | beeswarm_plots, "beeswarm and boxplots", NA),
                   ifelse(loess_plots & (fu_plots | age_plots | time_plots | temp_factor | hr_factor | rr_factor), "loess plots", NA),
                   ifelse(ind_plots & (fu_plots | age_plots | time_plots | temp_factor | hr_factor | rr_factor), "ind plots", NA))
    plot_types = plot_types[!is.na(plot_types)]
    
    
    ## Create directories
    dir.create(paste0(baseroute, "Results"), showWarnings = F)
    dir.create(paste0(baseroute, "Results/", study), showWarnings = F)
    route = paste0(baseroute, "Results/", study, "/", technique, 
                   ifelse(any(c("ECG", "Telemetry") %in% technique), paste0(" ", channel), ""), 
                   ifelse(technique == "Pressure", paste0(" -", daysremoved, " days ", outrounds, " outrounds", ifelse(removehighpressure, paste0(" max ", pressurethreshold, " mmHg"), ""), ifelse(daymedians, " daymedians", ""), ifelse(outpaired, " pairedouts", "")), ""),
                   ifelse(telemetryexp & cycle_sel != "All", paste0(" ", cycle_sel), ""),
                   ifelse(telemetryexp & torpor_sel != "All", paste0(" torpor ", torpor_sel), ""),
                   ifelse(telemetryexp & phase_sel != "All", paste0(" ", phase_sel), ""),
                   ifelse(weight_factor, " weight", ""),
                   ifelse(temp_factor, " temp", ""),
                   ifelse(hr_factor, " HR", ""),
                   ifelse(rr_factor, " RR", ""),
                   ifelse(sex_sel != "All", paste0(" sex ", sex_sel), "")
    )
    dir.create(route, showWarnings = F)
    for (i in plot_types) dir.create(paste0(route, "/", i), showWarnings = F)
    for (i in c(beeboxplots, followupboxplots, timeboxplots, flfollowupboxplots, fltimeboxplots)) dir.create(paste0(route, "/beeswarm and boxplots/", i), showWarnings = F)
    for (i in loessplots) dir.create(paste0(route, "/loess plots/", i), showWarnings = F)
    for (i in indplots) dir.create(paste0(route, "/ind plots/", i), showWarnings = F)
    dir.create(paste0(route, "/Diagnostic plots/"), showWarnings = F)
    
    ## Response variable and plot y labels
    vars = names(rawdata)
    exclvars = c("Id", "Date", "Time", "Absolute.time", "Event", "Cycle", "Torpor", "Phase")
    exclfacs = c()
    if (temp_factor) exclfacs = c(exclfacs, grep("temperature", vars, ignore.case = T, value = T))
    if (weight_factor) exclfacs = c(exclfacs, grep("weight", vars, ignore.case = T, value = T)[1]) # Include first weight only
    if (hr_factor) exclfacs = c(exclfacs, grep("heart.rate|hr", vars, ignore.case = T, value = T))
    if (rr_factor) exclfacs = c(exclfacs, grep("Delta.rr|^rr", vars, ignore.case = T, value = T))
    exclvars = c(exclvars, exclfacs)
    
    var_ommit = function (x) {
      var_subjects = suppressWarnings(min(summary(factor(expdata$Condition[!is.na(expdata[,x])]))))
      var_levels = suppressWarnings(length(levels(factor(expdata$Condition[!is.na(expdata[,x])]))))
      var_subjects == 1 | var_subjects == Inf | var_levels == 1
    }
    var_rm = which(vars %in% exclvars)
    var_rm = c(var_rm, which(sapply(vars, var_ommit) == T))
    vars = vars[-var_rm]
    
    collabels = as.character(fread(paste0(baseroute, "Raw data/", technique, ifelse(any(c("ECG", "Telemetry") %in% technique), paste0("/analyzed ", channel, "/"), "/"), "rawdata.txt"), encoding = "Latin-1", header = F, stringsAsFactors = F, dec = ".", na.strings = "N/A")[1,])
    collabels = gsub("\\([[:print:]]*C)", "( C)", collabels)
    ylabels = collabels[-var_rm]
    ylabels = gsub("\\([[:print:]]*C)", "( C)", ylabels)
    if (sex_sel != "All") expdata = expdata[expdata$Sex == sex_sel,]
    if (telemetryexp & cycle_sel != "All") expdata = expdata[expdata$Cycle == cycle_sel,]
    if (telemetryexp & torpor_sel != "All") expdata = expdata[expdata$Torpor == torpor_sel,]
    if (telemetryexp & technique != "Metabolic.cages" & phase_sel != "All") expdata = expdata[expdata$Phase == phase_sel,]
    
    idmeandata = expdata[!duplicated(expdata$Id), c(names(designdata), "Age0", exclfacs, vars, sort(unique(expdata$Condition)))]
    if (length(vars) > 1) {
      idmeandata[,vars] = t(sapply(unique(expdata$Id), function (x) colMeans(expdata[expdata$Id == x,vars], na.rm = T)))
    } else  idmeandata[,vars] = sapply(unique(expdata$Id), function (x) colMeans(data.frame(expdata[expdata$Id == x,vars]), na.rm = T))
    
    ## Stats and plots for each response variable
    allstats = data.frame(matrix(ncol = 1))[-1,]
    allint = data.frame(matrix(ncol = 1))[-1,]
    i = vars[1]
    i = vars[2]
    i = vars[length(vars)]
    for (i in vars){
      tryCatch({
        ## Determine groups, factors (fixed effects), random effects and labels for plots
        groups = unique(expdata[!is.na(expdata[,i]), "Condition"][order(expdata[!is.na(expdata[,i]), "Condition.Id"])])
        xlabel = groups # review
        colors = colfun(length(groups))
        if ("Color" %in% names(expdata)) colors = unique(expdata$Color[order(expdata$Condition.Id)])
        varlabel = ylabels[which(vars == i)]
        if (grepl("\\( C)", varlabel)) {
          if (grepl("Delta", varlabel)) {
            varlabel = as.expression(bquote(Delta*.(gsub("Delta", "", gsub(" C\\)", "", varlabel)))*degree*"C)"))
          } else varlabel = as.expression(bquote(.(gsub(" C\\)", "", varlabel))*degree*"C)"))
        } else if (grepl("Delta", varlabel)) {
          varlabel = as.expression(bquote(Delta*.(gsub("Delta", "", varlabel))))
        } else varlabel = as.expression(bquote(.(varlabel)))
        
        checkvar = function (x) grep("Delta", grep(x, names(expdata), value = T, ignore.case = T), value = T, invert = T)
        checklevels = function (y) all(sapply(groups, function (x) length(unique(expdata[expdata$Condition == x & !is.na(expdata[,i]), y]))) > 1)
        
        allfactors = c()
        if (any(grepl("^Time", names(rawdata)))) allfactors = c(allfactors, "^Time$")
        if (age_factor) allfactors = c(allfactors, "^Age$")
        if (sex_factor) allfactors = c(allfactors, "Sex")
        if (temp_factor) allfactors = c(allfactors, "Temperature.*C")
        if (weight_factor) allfactors = c(allfactors, "Weight..g")
        if (hr_factor) allfactors = c(allfactors, "heart.rate..b|hr..b")
        if (rr_factor) allfactors = c(allfactors, "RR..ms")
        if (telemetryexp) allfactors = c(allfactors, "Torpor", "Cycle")
        if (telemetryexp & technique != "Metabolic.cages") allfactors = c(allfactors, "Phase")
        if (date_factor & date_fixed) allfactors = c(allfactors, "^Date")
        
        factors = c()
        for (fac in allfactors) if (checklevels(checkvar(fac))) factors = c(factors, checkvar(fac)[1])
        
        randeff = c()
        if (any(duplicated(expdata[!is.na(expdata[,i]), "Id"]))) randeff = c(randeff, "Id")
        # if (!any(grepl("^Time", factors)) & checklevels("Follow.up")) randeff = c(randeff, "Follow.up")
        if (checklevels("Follow.up")) randeff = c(randeff, "Follow.up")
        if (date_factor & !date_fixed & checklevels("Date")) randeff = c(randeff, "Date")
        
        ## Deal with 0 values for transformations
        # expdata[!is.na(expdata[,i]) & expdata[,i] == 0, i] = NA
        # expdata[!is.na(expdata[,i]) & expdata[,i] == 0, i] = min(abs(expdata[!is.na(expdata[,i]) & expdata[,i] != 0, i]))/100
        # expdata[!is.na(expdata[,i]) & abs(expdata[,i]) == min(abs(expdata[!is.na(expdata[,i]),i])), i] = NA
        ## Stats
        lmstats = data.frame(matrix(ncol = 1))[-1,]
        lmint = data.frame(matrix(ncol = 1))[-1,]
        # for (vs in 1:(length(groups)-1)) {
        for (vs in 1:length(groups)) {
          if (length(unique(expdata[,i])) == 2 & length(randeff) == 0) {
            lm = glm(paste(i, " ~", paste(c(factors, groups[-vs]), collapse = " + ")), expdata, family = "binomial")
            Model = "Glm binomial"
            if (length(residuals(lm)) >= 5000) ShapiroWilk = 1 else ShapiroWilk = shapiro.test(residuals(lm))[[2]]
            BreuschPagan = bptest(lm)[[4]]
            lmcoef = summary(lm)$coef
          } else if (length(unique(expdata[,i])) == 2 & length(randeff) > 0) {
            lm = glmer(paste(i, " ~", paste(c(factors, groups[-vs]), collapse = " + "), paste0(" + (1|",paste0(randeff, collapse = ") + (1|"),")")), expdata, family = "binomial")
            Model = "Glmer binomial"
            if (length(residuals(lm)) >= 5000) ShapiroWilk = 1 else ShapiroWilk = shapiro.test(residuals(lm))[[2]]
            BreuschPagan = 1
            lmcoef = summary(lm)$coef
          } else if (length(randeff) == 0) {
            lm = lm(paste(i, "~", paste(c(factors, groups[-vs]), collapse = " + ")), expdata)
            Model = "Lm"
            if (length(residuals(lm)) >= 5000) ShapiroWilk = 1 else ShapiroWilk = shapiro.test(residuals(lm))[[2]]
            BreuschPagan = bptest(lm)[[4]]
            if ((ShapiroWilk < 0.05 | BreuschPagan < 0.05) & all(expdata[!is.na(expdata[,i]),i] != 0) & (all(expdata[!is.na(expdata[,i]),i] < 0) | all(expdata[!is.na(expdata[,i]),i] > 0))) {
              lm = lm(paste("log(abs(", i, ")) ~", paste(c(factors, groups[-vs]), collapse = " + ")), expdata)
              Model = "Lm log"
              if (length(residuals(lm)) >= 5000) ShapiroWilk = 1 else ShapiroWilk = shapiro.test(residuals(lm))[[2]]
              BreuschPagan = bptest(lm)[[4]]
              if (ShapiroWilk < 0.05 | BreuschPagan < 0.05) {
                lm = lm(paste("boxcoxt(", i, ") ~", paste(c(factors, groups[-vs]), collapse = " + ")), expdata)
                Model = "Lm BoxCox"
                if (length(residuals(lm)) >= 5000) ShapiroWilk = 1 else ShapiroWilk = shapiro.test(residuals(lm))[[2]]
                BreuschPagan = bptest(lm)[[4]]
              }
            }
            lmcoef = summary(lm)$coef
            
          } else if (length(randeff) > 0) {
            lm = lmer(paste(i, " ~", paste(c(factors, groups[-vs]), collapse = " + "), paste0(" + (1|",paste0(randeff, collapse = ") + (1|"),")")), expdata, REML = F)
            Model = "Lme"
            if (length(residuals(lm)) >= 5000) ShapiroWilk = 1 else ShapiroWilk = shapiro.test(residuals(lm))[[2]]
            BreuschPagan = 1
            if ((ShapiroWilk < 0.05 | BreuschPagan < 0.05) & all(expdata[!is.na(expdata[,i]),i] != 0) & (all(expdata[!is.na(expdata[,i]),i] < 0) | all(expdata[!is.na(expdata[,i]),i] > 0))) {
              lm = lmer(paste("log(abs(", i, ")) ~", paste(c(factors, groups[-vs]), collapse = " + "), paste0(" + (1|",paste0(randeff, collapse = ") + (1|"),")")), expdata, REML = F)
              Model = "Lme log"
              if (length(residuals(lm)) >= 5000) ShapiroWilk = 1 else ShapiroWilk = shapiro.test(residuals(lm))[[2]]
              BreuschPagan = 1
              if (ShapiroWilk < 0.05 | BreuschPagan < 0.05) {
                lm = lmer(paste("boxcoxt(", i, ") ~", paste(c(factors, groups[-vs]), collapse = " + "), paste0(" + (1|",paste0(randeff, collapse = ") + (1|"),")")), expdata, REML = F)
                Model = "Lme BoxCox"
                if (length(residuals(lm)) >= 5000) ShapiroWilk = 1 else ShapiroWilk = shapiro.test(residuals(lm))[[2]]
                BreuschPagan = 1
              }
            }
            lmcoef = summary(lm)$coef
          }
          Call = as.character(terms(lm))
          if (length(Call) > 1) Call = paste(Call[c(2,1,3)], collapse = " ")
          if (length(randeff) > 0) Call = paste(Call, paste0("+ (1|",paste0(randeff, collapse = ") + (1|"),")"), collapse = "")
          Warnings = c()
          if (ShapiroWilk < 0.05) Warnings = c(Warnings, "Non-normal")
          if (ShapiroWilk == 1) Warnings = c(Warnings, "Normality not evaluated")
          if (BreuschPagan < 0.05) Warnings = c(Warnings, "Non-equal variance")
          if (BreuschPagan == 1) Warnings = c(Warnings, "Variance constance not evaluated")
          Warnings[-1] = tolower(Warnings[-1])
          Warnings = paste(Warnings, collapse = ", ")
          assign(paste0("lmstats", vs), cbind.data.frame("Comparison" = c(groups[vs], factors, paste0(groups[-vs], sep = paste0("vs", groups[vs]))),
                                                         lmcoef, "2.5 %" = ifelse(confints, as.matrix(confint(lm)[rownames(lmcoef),1]), 0), "97.5 %" = ifelse(confints, as.matrix(confint(lm)[rownames(lmcoef),2]), 0), 
                                                         ShapiroWilk, BreuschPagan, Model, Call, Warnings, stringsAsFactors = F))
          lmstats = rbind.data.frame(lmstats,get(paste0("lmstats",vs))[-1,], make.row.names = F, stringsAsFactors = F)
          lmint = rbind.data.frame(lmint,get(paste0("lmstats",vs))[1,], make.row.names = F, stringsAsFactors = F)
        }
        names(lmstats) = gsub(" ", "", gsub("Pr.*", "p.value", gsub("t value", "t.value", names(lmstats))))
        names(lmint) = gsub("Comparison", "Group", gsub(" ", "", gsub("Pr.*", "p.value", gsub("t value", "t.value", names(lmstats)))))
        # names(lmstats)[3:5] = c("Std.Error", "t.value", "p.value")
        # names(lmint)[c(1,3:5)] = c("Group", "Std.Error", "t.value", "p.value")
        if(any(names(lmstats) == "df")) {
          lmstats = lmstats[,-which(names(lmstats) == "df")]
          lmint = lmint[,-which(names(lmint) == "df")]
        }
        lmstats = lmstats[!duplicated(lmstats$Comparison),]
        lmstats[,unlist(lapply(lmstats, is.numeric))] = round(lmstats[,unlist(lapply(lmstats, is.numeric))],5)
        lmstats = lmstats[!duplicated(abs(lmstats[,c("Estimate", "Std.Error", "p.value")])),]
        if (length(unique(expdata$Sex[!is.na(expdata[,i])])) > 1) lmstats[lmstats$Comparison == "Sex","Comparison"] = "MalevsFemale"
        lmstats$Sig = sapply(lmstats$p.value, function (x) ifelse(x < 0.001, "***", ifelse(x < 0.01, "**", ifelse(x < 0.05, "*", round(x, 4)))))
        lmstats$Sign = ifelse(lmstats$p.value < 0.05, sign(lmstats$Estimate), NA)
        
        allstats = rbind.data.frame(allstats, cbind.data.frame("Parameter" = ylabels[which(vars == i)], lmstats))
        allint = rbind.data.frame(allint, cbind.data.frame("Parameter" = ylabels[which(vars == i)], lmint))
        
        
          ## Stats for boxplots
          plotstats = lmstats[,c("Comparison","p.value","Sig")]
          if (length(factors) > 0) plotstats = plotstats[-1:-length(factors),]
          plotstats$G1 = sapply(gsub("vs[[:alnum:][:punct:]]*", "", plotstats$Comparison), function(x) which(groups == x))
          plotstats$G2 = sapply(gsub("[[:alnum:][:punct:]]*vs", "", plotstats$Comparison), function(x) which(groups == x))
          plotstats$P1 = sapply(1:nrow(plotstats), function (x) min(plotstats$G1[x],plotstats$G2[x]))
          plotstats$P2 = sapply(1:nrow(plotstats), function (x) max(plotstats$G1[x],plotstats$G2[x]))
          plotstats$Mean = rowMeans(data.frame(plotstats$P1, plotstats$P2))
          plotstats$Length = plotstats$P2 - plotstats$P1
          if (onlymainplotsigaxes) plotsigaxes = mainplotsigaxes else plotsigaxes = allplotsigaxes
          plotstats = merge(plotstats, plotsigaxes[plotsigaxes$Ngroups == length(groups),], by = c("P1","P2"), all.x = T)
          heightsigaxes = data.frame("Height" = sort(unique(plotstats[plotstats$p.value < 0.05, "Height"])))
          if (nrow(heightsigaxes) > 0) {
            heightsigaxes$Heightsig = seq(1, length(sort(unique(plotstats[plotstats$p.value < 0.05, "Height"]))), 1) - 1
            plotstats = merge(plotstats, heightsigaxes, by = c("Height"), all.x = T)
          } else plotstats$Heightsig = plotstats$Height  
          if (onlyplotsigcomps == F) {
            plotaxes = plotstats
          } else {
            plotaxes = plotstats[plotstats$p.value < 0.05 & !is.na(plotstats$Heightsig),]
            plotaxes$Height = plotaxes$Heightsig
          }
          
          ## Beeswarm and boxplots
          if (boxplots | beeswarm_plots) {
            for (k in beeboxplots){
              pdf(paste0(route, "/beeswarm and boxplots/", k, "/", which(vars == i), " ", i, " ", k, ".pdf"), (1.6 + 1.7*round(length(groups)/2))/0.9, 4/0.9, pointsize = 24, useDingbats = F)
              par(bty = "l", cex.axis = 0.75, mar = c(0.5,3,3.5,1), mgp = c(1.5,0.25,0), tck = - 0.03)
              if (strwidth(varlabel, units = "inches") > 2.5) {
                varlabel2 = varlabel
                varlabel = ""
              }
              if (k == "all beeswarm") {
                beeswarm(get(i) ~ Condition.Id, expdata, pch = 20, ylim = as.vector(summary(expdata[,i]))[c(1,6)], xaxt = "n", xlab ="", ylab = varlabel, add = F, col = colors, corral = "wrap")
                bxplot(get(i) ~ Condition.Id, expdata, add = T, width = 0.3)
              } else if (k == "mean beeswarm"){
                beeswarm(get(i) ~ Condition.Id, idmeandata, pch = 20, ylim = as.vector(summary(idmeandata[,i]))[c(1,6)], xaxt = "n", xlab ="", ylab = varlabel, add = F, col = colors, corral = "wrap")
                bxplot(get(i) ~ Condition.Id, idmeandata, add = T, width = 0.3)
              } else if (k == "all boxplot"){
                boxplot(get(i) ~ Condition.Id, expdata, pch = 20, ylim = as.vector(summary(expdata[,i]))[c(1,6)], xaxt = "n", xlab ="", ylab = varlabel, add = F, border = colors, boxlwd = 2.5, notch = boxnotch)
              } else if (k == "mean boxplot"){
                boxplot(get(i) ~ Condition.Id, idmeandata, pch = 20, ylim = as.vector(summary(idmeandata[,i]))[c(1,6)], xaxt = "n", xlab ="", ylab = varlabel, add = F, border = colors, boxlwd = 2.5, notch = boxnotch)
              } else if (k == "all boxbeeswarm"){
                boxplot(get(i) ~ Condition.Id, expdata, pch = 20, ylim = as.vector(summary(expdata[,i]))[c(1,6)], xaxt = "n", xlab ="", ylab = varlabel, add = F, border = colors, boxlwd = 2.5, boxlwd = 2.5, notch = boxnotch)
                beeswarm(get(i) ~ Condition.Id, expdata, pch = 20, ylim = as.vector(summary(expdata[,i]))[c(1,6)], xaxt = "n", xlab ="", ylab = varlabel, add = T, col = adjustcolor(colors, alpha.f = 0.2), corral = "wrap")
              } else if (k == "mean boxbeeswarm"){
                boxplot(get(i) ~ Condition.Id, idmeandata, pch = 20, ylim = as.vector(summary(idmeandata[,i]))[c(1,6)], xaxt = "n", xlab ="", ylab = varlabel, add = F, border = colors, boxlwd = 2.5, notch = boxnotch)
                beeswarm(get(i) ~ Condition.Id, idmeandata, pch = 20, ylim = as.vector(summary(idmeandata[,i]))[c(1,6)], xaxt = "n", xlab ="", ylab = varlabel, add = T, col = adjustcolor(colors, alpha.f = 0.2), corral = "wrap")
              }
              # axis(side = 1, at = 1:length(groups), labels = xlabel, tck = 0, lty = 0, mgp = c(3,1.5,0))
              # axis(side = 1, at = sort(unique(expdata$Condition.Id)), labels = xlabel, tck = 0, lty = 0, mgp = c(3,1.5,0))
              if (as.character(varlabel) == "") {
                mtext(side = 2, line = 1.5, varlabel2, adj = 0, padj = 0)
                varlabel = varlabel2
              }
              if (nrow(plotaxes) > 0) for (j in 1:nrow(plotaxes)) {
                axis(side = 3, line = plotaxes[j,"Height"]*0.75, at = plotaxes[j,c("P1","P2")] + c(0.03,-0.03), tck = 0, labels = F, lwd = 1.5)
                if (onlyplotsigcomps) {
                  # mtext(side = 3, line = plotaxes[j,"Height"]*0.75 - 0.3, at = plotaxes[j,"Mean"], plotaxes[j,"Sig"], cex = 0.8)
                  mtext(side = 3, line = plotaxes[j,"Height"]*0.75 - 0.4, at = plotaxes[j,"Mean"], plotaxes[j,"Sig"], cex = 1.2)
                } else mtext(side = 3, line = plotaxes[j,"Height"]*0.75 - 0.1, at = plotaxes[j,"Mean"], plotaxes[j,"Sig"], cex = 0.6)
              }
              dev.off()
            }
          }
          
          ## Follow-up and time beeswarm and boxplots
          timefuboxplots = NULL
          if (followup_boxplots & length(which(randeff == "Follow.up")) > 0) timefuboxplots = c(timefuboxplots, "followup")
          if (followup_boxplots & (length(which(randeff == "Time")) > 0 | length(which(factors == "Time")) > 0)) timefuboxplots = c(timefuboxplots, "time")
          for (l in timefuboxplots) {
            timeorfu = ifelse(l == "time", "TFollow.up", "Follow.up")
            # fupid = as.numeric(as.vector(sapply(1:length(unique(expdata[,timeorfu])), function (x) paste(x, 1:length(groups), sep = "."))))
            fupid = as.numeric(as.vector(sapply(unique(expdata[,timeorfu]), function (x) paste(x, 1:length(groups), sep = "."))))
            abscissa = ifelse(l == "time", "Time_Id", "Follow.up.Id")
            for (k in get(paste0(l,"boxplots"))) {
              pdf(paste0(route, "/beeswarm and boxplots/", k, "/", which(vars == i), " ", i, " ", k, ".pdf"), 7, 7, pointsize = 24, useDingbats = F)
              par(bty = "l", cex.axis = 0.7)
              if (k == paste(l,"beeswarm")) {
                beeswarm(get(i) ~ factor(get(abscissa), levels = fupid), expdata, pch = 20, ylim = as.vector(summary(expdata[,i]))[c(1,6)], xaxt = "n", xlab ="", ylab = varlabel, add = F, col = colors, corral = "wrap")
                bxplot(get(i) ~ factor(get(abscissa), levels = fupid), expdata, add = T, width = 0.3)
              } else if (k == paste(l,"boxplot")){
                boxplot(get(i) ~ factor(get(abscissa), levels = fupid), expdata, pch = 20, ylim = as.vector(summary(expdata[,i]))[c(1,6)], xaxt = "n", xlab ="", ylab = varlabel, add = F, border = colors, boxlwd = 2.5, notch = boxnotch)
              } else if (k == paste(l,"boxbeeswarm")){
                boxplot(get(i) ~ factor(get(abscissa), levels = fupid), expdata, pch = 20, ylim = as.vector(summary(expdata[,i]))[c(1,6)], xaxt = "n", xlab ="", ylab = varlabel, add = F, border = colors, boxlwd = 2.5, notch = boxnotch)
                beeswarm(get(i) ~ factor(get(abscissa), levels = fupid), expdata, pch = 20, ylim = as.vector(summary(expdata[,i]))[c(1,6)], xaxt = "n", xlab ="", ylab = varlabel, add = T, col = adjustcolor(colors, alpha.f = 0.2), corral = "wrap")
              }
              axis(side = 1, line = -1, at = 1:length(fupid), labels = rep(xlabel, length(unique(expdata[,timeorfu]))), tck = 0, lty = 0)
              # axis(side = 1, line = -1, at = 1:length(unique(expdata$Follow.up.Id)), labels = rep(xlabel, length(unique(expdata$Follow.up))), tck = 0, lty = 0)
              for (j in 1:length(unique(expdata[,timeorfu]))) {
                axis(side = 1, line = 1.5, at = c((j-1)*length(unique(expdata$Condition.Id)) + 1,j*length(unique(expdata$Condition.Id))), tck = 0, labels = F, lwd = 1)
                mtext(side = 1, line = 1.5, at = mean(c((j-1)*length(unique(expdata$Condition.Id)) + 1,j*length(unique(expdata$Condition.Id)))), paste("Follow-up", j), cex = 0.8)
              }
              dev.off()
            }
          }
          
          ## First-last follow-up beeswarm and boxplots
          fltimefuboxplots = NULL
          if (fl_followup_boxplots & length(which(randeff == "Follow.up")) > 0) fltimefuboxplots = c(fltimefuboxplots, "followup")
          if (fl_followup_boxplots & (length(which(randeff == "Time")) > 0 | length(which(factors == "Time")) > 0)) fltimefuboxplots = c(fltimefuboxplots, "time")
          for (l in fltimefuboxplots) {
            timeorfu = ifelse(l == "time", "Time", "Follow.up")
            abscissa = ifelse(l == "time", "Time_Id", "Follow.up.Id")
            ffu = min(expdata[!is.na(expdata[,i]),timeorfu])
            lfu = max(expdata[!is.na(expdata[,i]),timeorfu])
            for (k in get(paste0("fl",l,"boxplots"))){
              pdf(paste0(route, "/beeswarm and boxplots/", k, "/", which(vars == i), " ", i, " ", k, ".pdf"), 7, 7, pointsize = 24, useDingbats = F)
              par(bty = "l", cex.axis = 0.7)
              if (k == paste("first-last", l,"beeswarm")) {
                beeswarm(get(i) ~ get(abscissa), expdata[expdata[,timeorfu] == ffu | expdata[,timeorfu] == lfu,], pch = 20, ylim = as.vector(summary(expdata[,i]))[c(1,6)], xaxt = "n", xlab ="", ylab = varlabel, add = F, col = colors, corral = "wrap")
                bxplot(get(i) ~ get(abscissa), expdata[expdata[,timeorfu] == ffu | expdata[,timeorfu] == lfu,], add = T, width = 0.3)
              } else if (k == paste("first-last", l,"boxplot")){
                boxplot(get(i) ~ get(abscissa), expdata[expdata[,timeorfu] == ffu | expdata[,timeorfu] == lfu,], pch = 20, ylim = as.vector(summary(expdata[,i]))[c(1,6)], xaxt = "n", xlab ="", ylab = varlabel, add = F, border = colors, boxlwd = 2.5, notch = boxnotch)
              } else if (k == paste("first-last", l,"boxbeeswarm")){
                boxplot(get(i) ~ get(abscissa), expdata[expdata[,timeorfu] == ffu | expdata[,timeorfu] == lfu,], pch = 20, ylim = as.vector(summary(expdata[,i]))[c(1,6)], xaxt = "n", xlab ="", ylab = varlabel, add = F, border = colors, boxlwd = 2.5, notch = boxnotch)
                beeswarm(get(i) ~ get(abscissa), expdata[expdata[,timeorfu] == ffu | expdata[,timeorfu] == lfu,], pch = 20, ylim = as.vector(summary(expdata[,i]))[c(1,6)], xaxt = "n", xlab ="", ylab = varlabel, add = T, col = adjustcolor(colors, alpha.f = 0.2), corral = "wrap")
              }
              axis(side = 1, line = -1, at = 1:length(rep(xlabel, 2)), labels = rep(xlabel, 2), tck = 0, lty = 0)
              axis(side = 1, line = 1.5, at = c(1,length(unique(expdata$Condition.Id))), tck = 0, labels = F, lwd = 1)
              axis(side = 1, line = 1.5, at = c(length(unique(expdata$Condition.Id)) + 1,2*length(unique(expdata$Condition.Id))), tck = 0, labels = F, lwd = 1)
              mtext(side = 1, line = 1.5, at = mean(c(1,length(unique(expdata$Condition.Id)))), "First follow-up", cex = 0.8)
              mtext(side = 1, line = 1.5, at = mean(c(length(unique(expdata$Condition.Id)) + 1,2*length(unique(expdata$Condition.Id)))), "Last follow-up", cex = 0.8)
              dev.off()
            }
          }
          
          ## Loess plots
          if (!checklevels(checkvar("^Age$"))) loessplots = loessplots[-which(loessplots =="loess vs age")]
          if ((technique == "Metabolic.cages") & loess_plots) {
            plotdata = idmeandata
            plotdata$Age = plotdata$Age0
            loessplots = "loess vs age"
            if (weight_factor)  loessplots = c(loessplots, "loess vs weight")
          } else plotdata = expdata
          for (k in loessplots) {
            if (k == "loess vs day") { abscissa = "Day"; abslab = "Time (days)"
            } else if (k == "loess vs followup") { abscissa = "Follow.up"; abslab = "Follow-up"
            } else if (k == "loess vs time followup") { abscissa = "TFollow.up"; abslab = "Follow-up"
            } else if (k == "loess vs time") { abscissa = checkvar("^Time$"); abslab = "Time (min)"
            } else if (k == "loess vs weight") { abscissa = checkvar("Weight..g"); abslab = "Weight (g)"
            } else if (k == "loess vs temperature") { abscissa = checkvar("Temperature.*C"); abslab = "Temperature (C)"
            } else if (k == "loess vs heart rate") { abscissa = checkvar("heart.rate..b|hr..b"); abslab = "Heart rate (bpm)"
            } else if (k == "loess vs RR") { abscissa = checkvar("RR..ms"); abslab = "RR (ms)"
            } else { abscissa = "Age"; abslab = "Age (weeks)"
            }
            
            # names(plotdata)
            plotdata = plotdata[order(plotdata[,abscissa]),]
            pdf(paste0(route, "/loess plots/", k, "/", which(vars == i), " ", i, " ", k, ".pdf"), 10/0.9, 4/0.9, pointsize = 24, useDingbats = F)
            par(bty = "l", cex = 3/4, cex.axis = 0.75*4/3, cex.lab = 4/3, mar = c(3,3,1,1)*4/3, mgp = c(1.5,0.25,0)*4/3, tck = - 0.03)
            # par(bty = "l", cex.axis = 0.75, mar = c(3,3,1,1), mgp = c(1.5,0.25,0), tck = - 0.05)
            if (strwidth(varlabel, units = "inches") > 2.5) {
              varlabel2 = varlabel
              varlabel = ""
            }
            if (k == "loess vs followup" | k == "loess vs time followup") {
              plot(get(i) ~ get(abscissa), plotdata, xaxt="n", ylim = as.vector(summary(plotdata[,i]))[c(1,6)], xaxs = "r", yaxs = "r", xlab = abslab, ylab = varlabel, col = 0)
              axis(1, at = unique(plotdata[,abscissa]),labels = as.integer(unique(plotdata[,abscissa])))
            } else plot(get(i) ~ get(abscissa), plotdata, ylim = as.vector(summary(plotdata[,i]))[c(1,6)], xaxs = "r", yaxs = "r", xlab = abslab, ylab = varlabel, col = 0)
            for (j in 1:length(groups)) assign(paste0("loess",j), predict(loess(get(i) ~ get(abscissa), plotdata[plotdata$Condition.Id == j,], span = 1.25, na.action = na.exclude), se = T))
            for (j in 1:length(groups)) polygon(c(plotdata[!is.na(plotdata[,i]) & !is.na(plotdata[,abscissa]) & plotdata$Condition.Id == j,abscissa],
                                                  rev(plotdata[!is.na(plotdata[,i]) & !is.na(plotdata[,abscissa]) & plotdata$Condition.Id == j,abscissa])
            ), c(get(paste0("loess",j))$fit + qt(0.975,get(paste0("loess",j))$df)*get(paste0("loess",j))$se,
                 rev(get(paste0("loess",j))$fit - qt(0.975,get(paste0("loess",j))$df)*get(paste0("loess",j))$se)
            ), col = adjustcolor(colors[j], alpha.f = 0.3), border = NA)
            for (j in 1:length(groups)) points(get(i) ~ get(abscissa), plotdata[plotdata$Condition.Id == j,], pch = 20, col = adjustcolor(colors[j], alpha.f = 0.2))
            for (j in 1:length(groups)) lines(plotdata[!is.na(plotdata[,i]) & !is.na(plotdata[,abscissa]) & plotdata$Condition.Id == j, abscissa], get(paste0("loess",j))$fit, col = colors[j], lwd = 3)
            if (as.character(varlabel) == "") {
              mtext(side = 2, line = 1.5*4/3, varlabel2, adj = 1)
              varlabel = varlabel2
            }
            dev.off()
          }
          
          ## Individual plots
          for (k in indplots) {
            abscissa = ifelse(k == "ind vs time", "Time", ifelse(k == "ind vs day", "Day", "Age"))
            expdata = expdata[order(expdata[,abscissa]),]
            expdata = expdata[order(expdata$Id),]
            expdata = expdata[order(expdata$Condition.Id),]
            plots_per_col = ifelse(length(unique(expdata$Id)) == 3, 1, round(sqrt(length(unique(expdata$Id))),0))
            plots_per_row = ifelse(length(unique(expdata$Id)) == 3, 3, ceiling(length(unique(expdata$Id))/plots_per_col))
            pdf(paste0(route, "/ind plots/", k, "/", which(vars == i), " ", i, " ", k, ".pdf"), plots_per_row*2, plots_per_col*2, pointsize = 12, useDingbats = F)
            par(mfrow = c(plots_per_col,plots_per_row), bty = "l", cex.axis = 0.7, cex.main = 0.7, mar = c(2,2,0.5,0.5), mgp = c(2,1,0), oma = c(1,1,1,1))
            for (j in unique(expdata$Id)) {
              tryCatch({
                plot(get(i) ~ get(abscissa), expdata, ylim = as.vector(summary(expdata[,i]))[c(1,6)], xaxs = "r", yaxs = "r", col = 0)
                legend("top", legend = j, bty = "n")
                points(get(i) ~ get(abscissa), expdata[expdata$Id == j,], pch = 20, col = adjustcolor(colors[unique(expdata[expdata$Id == j, "Condition.Id"])], alpha.f = 0.5))
                # lines(get(i) ~ get(abscissa), expdata[expdata$Id == j & !is.na(expdata[,i]),], col = colors[unique(expdata[expdata$Id == j, "Condition.Id"])], lwd = 3)
                loessind = predict(loess(get(i) ~ get(abscissa), expdata[expdata$Id == j & !is.na(expdata[,i]),], span = 1.25, na.action = na.exclude), se = T)
                lines(expdata[expdata$Id == j & !is.na(expdata[,i]), abscissa], loessind$fit, col = colors[unique(expdata[expdata$Id == j, "Condition.Id"])], lwd = 3)
              }, error = function (e) cat())
            }
            title(paste0(ylabels[which(vars == i)], " vs ", ifelse (k == "ind vs time", "Time (min)", ifelse(k == "ind vs day", "Time (days)", "Age (weeks)"))), outer = T, cex.main = 1)
            dev.off()
          }
          
          ## Diagnostic plots
          pdf(paste0(route, "/Diagnostic plots/", which(vars == i), " ", i, " diagnostic plots.pdf"), useDingbats = F)
          par(mfrow = c(2,2), pch = 20)
          plot(residuals(lm) ~ fitted(lm), xlab = "Fitted values", ylab = "Residuals")
          mtext("Residuals vs Fitted", line = 0.2, side = 3)
          abline(h = 0, col = "gray", lty = 3)
          text(x = fitted(lm)[order(residuals(lm), decreasing = T)][1:3], y = sort(residuals(lm), decreasing = T)[1:3], labels = names(sort(residuals(lm), decreasing = T)[1:3]), cex = 0.75, pos = 2)
          lines(lowess(residuals(lm) ~ fitted(lm)),  col = 2)
          
          qq = qqnorm(residuals(lm), main = "", ylab = "Standardized residuals")
          mtext("Normal Q-Q", line = 0.2, side = 3)
          qqline(residuals(lm), lty = 3, col = "gray50")
          text(x = qq$x[order(residuals(lm), decreasing = T)][1:3], y = sort(residuals(lm), decreasing = T)[1:3], labels = names(sort(residuals(lm), decreasing = T)[1:3]), cex = 0.75, pos = 2)
          
          plot(sqrt(abs(stresiduals(lm))) ~ fitted(lm), xlab = "Fitted values", ylab = as.expression(substitute(sqrt(abs("Standardized residuals")))))
          mtext("Scale-Location", line = 0.2, side = 3)
          text(x = fitted(lm)[order(residuals(lm), decreasing = T)][1:3], y = sort(sqrt(abs(stresiduals(lm))), decreasing = T)[1:3], labels = names(sort(sqrt(abs(stresiduals(lm))), decreasing = T)[1:3]), cex = 0.75, pos = 2)
          lines(lowess(sqrt(abs(stresiduals(lm))) ~ fitted(lm)),  col = 2)
          
          plot(stresiduals(lm) ~ hatvalues(lm), xlab = "Leverage", ylab = "Standardized residuals", xlim = c(0,max(hatvalues(lm))))
          mtext("Residuals vs Leverage", line = 0.2, side = 3)
          text(x = hatvalues(lm)[order(stresiduals(lm), decreasing = T)][1:3], y = sort(stresiduals(lm), decreasing = T)[1:3], labels = names(sort(stresiduals(lm), decreasing = T)[1:3]), cex = 0.75, pos = 2)
          lines(lowess(stresiduals(lm) ~ hatvalues(lm)),  col = 2)
          abline(h = 0, v = 0, lty = 3, col = "gray")
          legend("bottomleft", legend = "Cook's distance", lty = 2, col = 2, bty = "n")
          cook.levels = c(0.5, 1)
          hh = seq.int(min(range(hatvalues(lm), na.rm = T)[1], range(hatvalues(lm), na.rm = T)[2]/100), par("usr")[2], length.out = 101)
          for (crit in cook.levels) {
            cl.h = sqrt(crit * nrow(summary(lm)$coef) * (1 - hh)/hh)
            lines(hh, cl.h, lty = 2, col = 2)
            lines(hh, -cl.h, lty = 2, col = 2)
          }
          aty = sqrt(cook.levels) * sqrt(nrow(summary(lm)$coef) * (1 - min(0.99, par("usr")[2]))/min(0.99, par("usr")[2]))
          axis(4, at = c(-rev(aty), aty), labels = paste(c(rev(cook.levels), cook.levels)), mgp = c(0.25, 0.25, 0), las = 2, tck = 0, cex.axis = 0.75, col.axis = 2)
          dev.off()
          
          
        }, error = function (e) cat("ERROR :", study, " ", technique, " ", i, " ", conditionMessage(e), "\n"))
      graphics.off()
      }
      
      ## Legend
      expdata = expdata[order(expdata$Condition.Id),]
      legtext = expdata[!duplicated(expdata$Condition.Id),"Genotype"]
      if (any(names(expdata) == "Treatment")) legtext = paste(legtext, expdata[!duplicated(expdata$Condition.Id),"Treatment"])
      if (any(names(expdata) == "Age.Group")) legtext = paste(legtext, expdata[!duplicated(expdata$Condition.Id),"Age.Group"])
      # legtext = paste(expdata[!duplicated(expdata$Condition.Id),"Genotype"], expdata[!duplicated(expdata$Condition.Id),"Treatment"], expdata[!duplicated(expdata$Condition.Id),"Age.group"]) # Review
      pdf(paste0(route, "/Legend.pdf"), max(strwidth(legtext, units = "inches")) + 2/5, (length(legtext) + 1)/5, useDingbats = F)
      par(mar = rep(0,4))
      plot(1, type = "n", axes = F, xlab = "", ylab = "")
      legend("center", pch = 16, cex = 1, legend = legtext, col = colors, bty = "n")
      dev.off(); dev.off()
      
      ## Sample size
      if (any(names(rawdata) == "Time") & !telemetryexp) {
        sampledata = as.data.frame(rbind(addmargins(table(expdata[expdata$TFollow.up == 1, c("Date", "Condition")])),
                                         N = addmargins(table(idmeandata[,"Condition"]))))
      } else if (telemetryexp) {
        sampledata = as.data.frame(rbind(N = addmargins(table(idmeandata[,"Condition"]))))
      } else sampledata = as.data.frame(rbind(addmargins(table(expdata[, c("Date", "Condition")])),
                                              N = addmargins(table(idmeandata[,"Condition"]))))
      
      
      ## Export data
      if (confints == F) allstats = allstats[, -which(names(allstats) == "2.5%"):-which(names(allstats) == "97.5%")]
      filterstats = allstats[allstats$p.value < 0.05,]
      
      fwrite(allstats, file = paste0(route, "/expstats.xls"), row.names = F, sep = "\t")
      fwrite(filterstats, file = paste0(route, "/fexpstats.xls"), row.names = F, sep = "\t")
      fwrite(allint, file = paste0(route, "/expintercepts.xls"), row.names = F, sep = "\t")
      
      
      names(expdata)[names(expdata) %in% names(rawdata)] = collabels[names(rawdata) %in% names(expdata)]
      fwrite(expdata, file = paste0(route, "/expdata.xls"), row.names = F, sep = "\t")
      names(idmeandata)[names(idmeandata) %in% names(rawdata)] = collabels[names(rawdata) %in% names(idmeandata)]
      fwrite(idmeandata, file = paste0(route, "/idmeanexpdata.xls"), row.names = F, sep = "\t")
      fwrite(sampledata, file = paste0(route, "/sampledata.txt"), row.names = T, sep = "\t")
      
      if (!file.exists(paste0(baseroute, "Results/processed.log.txt"))) fwrite(cbind.data.frame(t(names(rconfig)), "Processing.time"), file = paste0(baseroute, "Results/processed.log.txt"), row.names = F, col.names = F, sep = "\t", append = F)
      fwrite(cbind.data.frame(rconfig[exp,], Sys.time()), file = paste0(baseroute, "Results/processed.log.txt"), row.names = F, col.names = F, sep = "\t", append = T)
      setTxtProgressBar(txtProgressBar(style = 3), exp/nrow(rconfig))
    }, error = function (e) cat("ERROR :", study, " ", technique, " ", conditionMessage(e), "\n"))
  }
  time1 = (proc.time() - time0)[[3]]
  paste0(round(round(round(time1/60)/60)/24), "d ", round(round(time1/60)/60)%%24, "h ", round(time1/60)%%60, "min ", round(time1%%60), "s")
  