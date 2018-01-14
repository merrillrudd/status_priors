
rm(list=ls())

########################################
## packages
########################################
library(ggplot2)
library(dplyr)
library(reshape2)
# library(nlme)

########################################
## directories
########################################
wd <- "C:\\merrill\\status_priors"

fig_dir <- file.path(wd, "figures")
dir.create(fig_dir, showWarnings=FALSE)

data_dir <- file.path(wd, "data")

########################################
## figure theme
########################################
theme_lsd <- function (base_size = 14, base_family = "") 
{
    theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(axis.title.x = element_text(margin = margin(10,0,0,0)),
          #axis.title.x = element_text(vjust = -1.5),
          #axis.title.y = element_text(margin = margin(0,20,0,0)),
          #axis.title.y = element_text(vjust = -0.1),
          axis.text = element_text(size = rel(0.8)),
          axis.ticks = element_line(colour = "black"), 
          legend.key = element_rect(colour = "grey80"),
          panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey50"),
          panel.grid.major = element_line(colour = "grey90", size = 0.2),
          panel.grid.minor = element_line(colour = "grey98", size = 0.5),
          strip.background = element_rect(fill = "grey80", colour = "grey50", size = 0.2))
}

########################################
## combine information
########################################
timeseries <- read.csv(file.path(data_dir, "RLSADB_v4.25_timeseries.csv"), stringsAsFactors=FALSE, header=TRUE)
	timeseries <- melt(timeseries, measure.vars=c("BdivBmsypref", "UdivUmsypref"))
	timeseries <- timeseries[-which(is.na(timeseries$value)),]
stockinfo <- read.csv(file.path(data_dir, "RLSADB_v4.25_stockinfo.csv"), stringsAsFactors=FALSE, header=TRUE)
	colnames(stockinfo)[1] <- "stockid"
	stockinfo <- stockinfo %>% dplyr::select("stockid", "region", "scientificname")
	stockinfo <- stockinfo[which(stockinfo$stockid %in% timeseries$stockid),]

bioparams <- read.csv(file.path(data_dir, "RLSADB_v4.25_bioparams.csv"), stringsAsFactors=FALSE)
	bioparams <- bioparams[which(bioparams$stockid %in% timeseries$stockid),] %>%
				dplyr::mutate(M = as.numeric(M))

taxonomy <- read.csv(file.path(data_dir, "RLSADB_v4.25_taxonomy.csv"), stringsAsFactors=FALSE)
	taxonomy <- taxonomy[which(taxonomy$scientificname %in% stockinfo$scientificname),]

data <- timeseries %>% 
		dplyr::select(stockid, stocklong, year, variable, value) %>%
		dplyr::mutate(variable = ifelse(variable=="BdivBmsypref","B/Bmsy", ifelse(variable=="UdivUmsypref","U/Umsy",NA))) %>%
		dplyr::full_join(stockinfo, by="stockid") %>%
		dplyr::full_join(bioparams %>% dplyr::select(stockid, M), by="stockid") %>%
		dplyr::full_join(taxonomy %>% dplyr::select(scientificname, FisheryType, taxGroup), by="scientificname") %>%
		dplyr::mutate(value = ifelse(value > 10, 10, value)) %>%
		na.omit()
data$region <- as.factor(data$region)
data$FisheryType <- as.factor(data$FisheryType)
data$taxGroup <- as.factor(data$taxGroup)

region_vec <- unique(data$region)
fishtype_vec <- unique(data$FisheryType)
taxgroup_vec <- unique(data$taxGroup)

data <- data %>% 
		dplyr::mutate(region2 = as.factor(match(region, region_vec))) %>%
		dplyr::mutate(FisheryType2 = as.factor(match(FisheryType, fishtype_vec))) %>%
		dplyr::mutate(taxGroup2 = as.factor(match(taxGroup, taxgroup_vec))) %>%
		dplyr::mutate(logvalue = log(value)) %>%
		dplyr::mutate(logM = log(M))

## B/Bmsy averaged across last 5 years
bdata <- data %>% filter(variable == 'B/Bmsy')
bstocks <- unique(bdata$stockid)
	bdata_sum <- lapply(1:length(bstocks), function(x){
		sub <- bdata %>% filter(stockid == bstocks[x])
		yr <- sub$year
		if(length(yr)>=5) sub2 <- sub %>% filter(year %in% yr[(length(yr)-4):length(yr)])
		if(length(yr)<5) sub2 <- sub
		avg <- mean(sub2$value)
		sub3 <- sub2 %>% mutate(value = avg) %>%
				mutate(logvalue = log(avg)) %>%
				select(-year)
		return(unique(sub3))
	})
sbdata <- do.call(rbind, bdata_sum)

## U/Umsy averaged across last 5 years
udata <- data %>% filter(variable == 'U/Umsy')
ustocks <- unique(udata$stockid)
	udata_sum <- lapply(1:length(ustocks), function(x){
		sub <- udata %>% filter(stockid == ustocks[x])
		yr <- sub$year
		if(length(yr)>=5) sub2 <- sub %>% filter(year %in% yr[(length(yr)-4):length(yr)])
		if(length(yr)<5) sub2 <- sub
		avg <- mean(sub2$value)
		sub3 <- sub2 %>% mutate(value = avg) %>%
				mutate(logvalue = log(avg)) %>%
				select(-year)
		return(unique(sub3))
	})
sudata <- do.call(rbind, udata_sum)
sdata <- rbind(sbdata, sudata)
allstocks <- unique(sdata %>% select(stockid, M, logM, region, FisheryType, taxGroup))

#########################################
## plot the data against single factors
#########################################
## factors -- region of ocean, type of fish, natural mortality rate
	## scatterplots
	p <- ggplot(sdata) +
		geom_point(aes(x=M, y=value)) +
		stat_smooth(aes(x=M, y=value), method="lm") +
		facet_grid(variable ~ ., scales="free") +
		geom_hline(yintercept=1, color="black", lty=2) +
		theme_lsd()
	ggsave(file.path(fig_dir, "scatterplots_M.png"), p)

	p <- ggplot(sdata) +
		geom_point(aes(x=logM, y=logvalue)) +
		stat_smooth(aes(x=logM, y=logvalue), method="lm") +
		facet_grid(variable ~ ., scales="free") +
		geom_hline(yintercept=0, color="black", lty=2) +
		theme_lsd()
	ggsave(file.path(fig_dir, "scatterplots_logM.png"), p)

	## boxplots
	p <- ggplot(sdata) +
		geom_boxplot(aes(x=variable, y=value, color=variable, fill=variable)) +
		geom_hline(yintercept=1, color="black", lty=2) +
		guides(fill=FALSE, color=FALSE) +
		theme_lsd()
	ggsave(file.path(fig_dir, "violin_all.png"), p)

	p <- ggplot(sdata) +
		geom_boxplot(aes(x=variable, y=logvalue, color=variable, fill=variable)) +
		geom_hline(yintercept=0, color="black", lty=2) +
		guides(fill=FALSE, color=FALSE) +
		theme_lsd()
	ggsave(file.path(fig_dir, "violin_logall.png"), p)

	p <- ggplot(sdata) +
		geom_boxplot(aes(x=region, y=value, color=region, fill=region)) +
		facet_grid(variable ~ ., scales="free") + 
		geom_hline(yintercept=1, color="black", lty=2) +
		theme_lsd() +
		theme(axis.text.x=element_blank())
	ggsave(file.path(fig_dir, "violin_region.png"), p, width=10)

	
	p <- ggplot(sdata) +
		geom_boxplot(aes(x=region, y=logvalue, color=region, fill=region)) +
		facet_grid(variable ~ ., scales="free") + 
		geom_hline(yintercept=0, color="black", lty=2) +
		theme_lsd() +
		theme(axis.text.x=element_blank())
	ggsave(file.path(fig_dir, "violin_logregion.png"), p, width=10)	

	p <- ggplot(sdata) +
		geom_boxplot(aes(x=FisheryType, y=value, color=FisheryType, fill=FisheryType)) +
		facet_grid(variable ~ ., scales="free") + 
		geom_hline(yintercept=1, color="black", lty=2) +
		theme_lsd() +
		theme(axis.text.x=element_blank())
	ggsave(file.path(fig_dir, "violin_fishtype.png"), p, width=10)

	p <- ggplot(sdata) +
		geom_boxplot(aes(x=FisheryType, y=logvalue, color=FisheryType, fill=FisheryType)) +
		facet_grid(variable ~ ., scales="free") + 
		geom_hline(yintercept=0, color="black", lty=2) +
		theme_lsd() +
		theme(axis.text.x=element_blank())
	ggsave(file.path(fig_dir, "violin_logfishtype.png"), p, width=10)

	## correlation between factors
	p <- ggplot(allstocks) + 
		geom_boxplot(aes(x=region, y=M, color=region, fill=region)) +
		theme_lsd() +
		theme(axis.text.x=element_blank())
	ggsave(file.path(fig_dir, "region_by_M.png"), p, width=10)

	p <- ggplot(allstocks) + 
		geom_boxplot(aes(x=region, y=logM, color=region, fill=region)) +
		theme_lsd() +
		theme(axis.text.x=element_blank())
	ggsave(file.path(fig_dir, "region_by_M.png"), p, width=10)

	p <- ggplot(allstocks) + 
		geom_boxplot(aes(x=FisheryType, y=M, color=FisheryType, fill=FisheryType)) +
		theme_lsd() +
		theme(axis.text.x=element_blank())
	ggsave(file.path(fig_dir, "fishtype_by_M.png"), p, width=10)

	p <- ggplot(allstocks) + 
		geom_boxplot(aes(x=FisheryType, y=logM, color=FisheryType, fill=FisheryType)) +
		theme_lsd() +
		theme(axis.text.x=element_blank())
	ggsave(file.path(fig_dir, "fishtype_by_logM.png"), p, width=10)

	png(file.path(fig_dir, "factor_correlation_logM.png"), height=8, width=10, units="in", res=200)
	pairs(data.frame("logM"=allstocks$logM, "FishType"=allstocks$FisheryType, "Region"=allstocks$region))
	dev.off()

	png(file.path(fig_dir, "factor_correlation.png"), height=8, width=10, units="in", res=200)
	pairs(data.frame("logM"=allstocks$M, "FishType"=allstocks$FisheryType, "Region"=allstocks$region))
	dev.off()

#########################################
## regression
#########################################

## U/UMSY
u1 <- lm(logvalue ~ logM, data=sudata)
summary(u1)
plot(sudata$value, col="gray", lwd=2)
lines(fitted(u1))
boxplot(resid(u1))

u2 <- lm(logvalue ~ region2, data=sudata)
summary(u2)
plot(sudata$logvalue, col="gray", lwd=2)
lines(fitted(u2))

u3 <- lm(logvalue ~ FisheryType2, data=sudata)
summary(u3)
plot(sudata$logvalue, col="gray", lwd=2)
lines(fitted(u3))
boxplot(resid(u3))

u4 <- lm(logvalue ~ taxGroup2, data=sudata)
summary(u4)

u5 <- lm(logvalue ~ region2 + M, data=sudata)
summary(u5)
plot(sudata$logvalue, col="gray", lwd=2)
lines(fitted(u5))
boxplot(resid(u5))

u6 <- lm(logvalue ~ FisheryType2 + M, data=sudata)
summary(u6)
plot(sudata$logvalue, col="gray", lwd=2)
lines(fitted(u6))
boxplot(resid(u6))

u7 <- lm(logvalue ~ taxGroup2 + M, data=sudata)
summary(u7)

u8 <- lm(logvalue ~ region2 + FisheryType2, data=sudata)
summary(u8)
plot(sudata$logvalue, col="gray", lwd=2)
lines(fitted(u8))

u9 <- lm(logvalue ~ region2 + taxGroup2, data=sudata)
summary(u9)

u10 <- lm(logvalue ~ M + region2 + FisheryType2, data=sudata)
summary(u10)
plot(sudata$logvalue, col="gray", lwd=2)
lines(fitted(u10))

u11 <- lm(logvalue ~ M + region2 + taxGroup2, data=sudata)
summary(u11)

umods <- list("M"=u1, "region"=u2, "fishtype"=u3, "taxagroup"=u4,
				"region_M"=u5, "fishtype_M"=u6, "taxagroup_M"=u7, 
				"region_fishtype"=u8, "region_taxagroup"=u9, 
				"all_fishtype"=u10, "all_taxagroup"=u11)
ur2 <- sapply(1:length(umods), function(x) summary(umods[[x]])$adj.r.squared)

aicu <- AIC(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11) %>%
		dplyr::mutate(deltaAIC = AIC - min(AIC))
aicu1 <- AIC(u1, u2, u3, u4) %>%
		dplyr::mutate(deltaAIC = AIC - min(AIC))

aicu$R2 <- ur2
aicu$ModelName <- names(umods)
aicu1$ModelName <- names(umods)[1:4]
write.csv(aicu, file.path(fig_dir, "UUmsy_AIC.csv"))

### B/Bmsy
b1 <- lm(logvalue ~ logM, data=sbdata)
summary(b1)
plot(sbdata$value, col="gray", lwd=2)
lines(fitted(b1))
boxplot(resid(b1))

b2 <- lm(logvalue ~ region2, data=sbdata)
summary(b2)
plot(sbdata$logvalue, col="gray", lwd=2)
lines(fitted(b2))

b3 <- lm(logvalue ~ FisheryType2, data=sbdata)
summary(b3)
plot(sbdata$logvalue, col="gray", lwd=2)
lines(fitted(b3))
boxplot(resid(b3))

b4 <- lm(logvalue ~ taxGroup2, data=sbdata)
summary(b4)

b5 <- lm(logvalue ~ region2 + M, data=sbdata)
summary(b5)
plot(sbdata$logvalue, col="gray", lwd=2)
lines(fitted(b5))
boxplot(resid(b5))

b6 <- lm(logvalue ~ FisheryType2 + M, data=sbdata)
summary(b6)
plot(sbdata$logvalue, col="gray", lwd=2)
lines(fitted(b6))
boxplot(resid(b6))

b7 <- lm(logvalue ~ taxGroup2 + M, data=sbdata)
summary(b7)

b8 <- lm(logvalue ~ region2 + FisheryType2, data=sbdata)
summary(b8)
plot(sbdata$logvalue, col="gray", lwd=2)
lines(fitted(b8))

b9 <- lm(logvalue ~ region2 + taxGroup2, data=sbdata)
summary(b9)

b10 <- lm(logvalue ~ M + region2 + FisheryType2, data=sbdata)
summary(b10)
plot(sbdata$logvalue, col="gray", lwd=2)
lines(fitted(b10))

b11 <- lm(logvalue ~ M + region2 + taxGroup2, data=sbdata)
summary(b11)

bmods <- list("M"=b1, "region"=b2, "fishtype"=b3, "taxagroup"=b4,
				"region_M"=b5, "fishtype_M"=b6, "taxagroup_M"=b7, 
				"region_fishtype"=b8, "region_taxagroup"=b9, 
				"all_fishtype"=b10, "all_taxagroup"=b11)
br2 <- sapply(1:length(bmods), function(x) summary(bmods[[x]])$adj.r.squared)

aicb <- AIC(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11) %>%
		dplyr::mutate(deltaAIC = AIC - min(AIC))
aicb1 <- AIC(b1, b2, b3, b4) %>%
		dplyr::mutate(deltaAIC = AIC - min(AIC))

aicb$R2 <- br2
aicb$ModelName <- names(bmods)
aicb1$ModelName <- names(bmods)[1:4]
write.csv(aicb, file.path(fig_dir, "UUmsy_AIC.csv"))


############################
## diagnostic plots
############################
png(file.path(fig_dir, "diagnostic_BM.png"), height=6, width=8, units="in", res=200)
par(mfrow=c(2,2))
plot(b1)
dev.off()

png(file.path(fig_dir, "diagnostic_Bregion.png"), height=6, width=8, units="in", res=200)
par(mfrow=c(2,2))
plot(b2)
dev.off()

png(file.path(fig_dir, "diagnostic_Bfishtype.png"), height=6, width=8, units="in", res=200)
par(mfrow=c(2,2))
plot(b3)
dev.off()

png(file.path(fig_dir, "diagnostic_B_all_fishtype.png"), height=6, width=8, units="in", res=200)
par(mfrow=c(2,2))
plot(b10)
dev.off()

png(file.path(fig_dir, "diagnostic_UM.png"), height=6, width=8, units="in", res=200)
par(mfrow=c(2,2))
plot(u1)
dev.off()

png(file.path(fig_dir, "diagnostic_Uregion.png"), height=6, width=8, units="in", res=200)
par(mfrow=c(2,2))
plot(u2)
dev.off()

png(file.path(fig_dir, "diagnostic_Ufishtype.png"), height=6, width=8, units="in", res=200)
par(mfrow=c(2,2))
plot(u3)
dev.off()

png(file.path(fig_dir, "diagnostic_U_fishtype_M.png"), height=6, width=8, units="in", res=200)
par(mfrow=c(2,2))
plot(u6)
dev.off()

png(file.path(fig_dir, "diagnostic_U_all_fishtype.png"), height=6, width=8, units="in", res=200)
par(mfrow=c(2,2))
plot(u10)
dev.off()

p <- ggplot(sdata, aes(x=logM, y=logvalue , colour=factor(region))) +
	geom_point() +
	stat_smooth(method=lm, fullrange=FALSE) + 
	facet_grid(variable ~ .)
ggsave(file.path(fig_dir, "logM_logvalue_byregion_lm.png"), p, width=10)
############################
## residual plots
############################
### natural mortality
bres1 <- data.frame(variable="B/BMSY", value=b1$residuals)
ures1 <- data.frame(variable="U/UMSY", value=u1$residuals)
mres1 <- rbind(bres1, ures1)

p <- ggplot(mres1) +
	geom_violin(aes(x=variable, y=value, color=variable, fill=variable)) +
	theme_lsd() +
	guides(color=FALSE, fill=FALSE) + 
	geom_hline(yintercept=0, color="black", lty=2) +
	xlab("") + ylab("Residuals")
ggsave(file.path(fig_dir, "M_residuals.png"), p, width=6)

### region
b2 <- fortify(b2)
u2 <- fortify(u2)
bres2 <- data.frame(variable="B/BMSY", value=b2$.resid, region=b2$region)
ures2 <- data.frame(variable="U/UMSY", value=u2$.resid, region=u2$region)
mres2 <- rbind(bres2, ures2)

p <- ggplot(mres2) +
	geom_violin(aes(x=region, y=value, color=region, fill=region)) +
	facet_grid(variable ~ .) +
	theme_lsd() +
	geom_hline(yintercept=0, color="black", lty=2) +
	ylab("Residuals") +
	theme(axis.text.x=element_blank())
ggsave(file.path(fig_dir, "region_residuals.png"), p, width=10)

### fish type
b3 <- fortify(b3)
u3 <- fortify(u3)
bres3 <- data.frame(variable="B/BMSY", value=b3$.resid, FisheryType=b3$FisheryType)
ures3 <- data.frame(variable="U/UMSY", value=u3$.resid, FisheryType=u3$FisheryType)
mres3 <- rbind(bres3, ures3)

p <- ggplot(mres3) +
	geom_violin(aes(x=FisheryType, y=value, color=FisheryType, fill=FisheryType)) +
	facet_grid(variable ~ .) +
	theme_lsd() +
	geom_hline(yintercept=0, color="black", lty=2) +
	ylab("Residuals") +
	theme(axis.text.x=element_blank())
ggsave(file.path(fig_dir, "fishtype_residuals.png"), p, width=10)



############################
## cross-validation
############################

# ## randomly select half the stocks
bchoose_raw <- sample(1:length(bstocks), size=length(bstocks)/2)
bchoose <- bchoose_raw[order(bchoose_raw)]

bdata_train <- sbdata %>% dplyr::filter(stockid %in% bstocks[bchoose])
	bregions_choose <- unique(bdata_train$region)
	bfishtype_choose <- unique(bdata_train$FisheryType)
	btaxatype_choose <- unique(bdata_train$taxGroup)
bdata_test <- sbdata %>% dplyr::filter(stockid %in% bstocks[bchoose] == FALSE) %>%
			dplyr::filter(region %in% bregions_choose) %>%
			dplyr::filter(FisheryType %in% bfishtype_choose) %>%
			dplyr::filter(taxGroup %in% btaxatype_choose)

### run models with training data
cb1 <- lm(logvalue ~ logM, data=bdata_train)
summary(cb1)

cb2 <- lm(logvalue ~ region2, data=bdata_train)
summary(cb2)

cb3 <- lm(logvalue ~ FisheryType2, data=bdata_train)
summary(cb3)

cb4 <- lm(logvalue ~ taxGroup2, data=bdata_train)
summary(cb4)

cb5 <- lm(logvalue ~ region2 + logM, data=bdata_train)
summary(cb5)

cb6 <- lm(logvalue ~ FisheryType2 + logM, data=bdata_train)
summary(cb6)

cb7 <- lm(logvalue ~ taxGroup2 + logM, data=bdata_train)
summary(cb7)

cb8 <- lm(logvalue ~ region2 + FisheryType2, data=bdata_train)
summary(cb8)

cb9 <- lm(logvalue ~ region2 + taxGroup2, data=bdata_train)
summary(cb9)

cb10 <- lm(logvalue ~ logM + region2 + FisheryType2, data=bdata_train)
summary(cb10)

cb11 <- lm(logvalue ~ logM + region2 + taxGroup2, data=bdata_train)
summary(cb11)


cbmods <- list("M"=cb1, "region"=cb2, "fishtype"=cb3, "taxagroup"=cb4,
				"region_M"=cb5, "fishtype_M"=cb6, "taxagroup_M"=cb7, 
				"region_fishtype"=cb8, "region_taxagroup"=cb9,
				"all_fishtype"=cb10, "all_taxagroup"=cb11)
cbr2 <- sapply(1:length(cbmods), function(x) summary(cbmods[[x]])$adj.r.squared)

aicbc <- AIC(cb1, cb2, cb3, cb4, cb5, cb6, cb7, cb8, cb9, cb10, cb11) %>%
		dplyr::mutate(deltaAIC = AIC - min(AIC))
aicbc1 <- AIC(cb1, cb2, cb3, cb4) %>%
		dplyr::mutate(deltaAIC = AIC - min(AIC))
aicbc$R2 <- cbr2
aicbc$ModelName <- names(cbmods)
write.csv(aicb, file.path(fig_dir, "cross_BBmsy_AIC.csv"))


# ## randomly select half the stocks
uchoose_raw <- sample(1:length(ustocks), size=length(ustocks)/2)
uchoose <- uchoose_raw[order(uchoose_raw)]

udata_train <- sudata %>% dplyr::filter(stockid %in% ustocks[uchoose])
	uregions_choose <- unique(udata_train$region)
	ufishtype_choose <- unique(udata_train$FisheryType)
	utaxatype_choose <- unique(udata_train$taxGroup)
udata_test <- sudata %>% dplyr::filter(stockid %in% ustocks[uchoose] == FALSE) %>%
			dplyr::filter(region %in% uregions_choose) %>%
			dplyr::filter(FisheryType %in% ufishtype_choose) %>%
			dplyr::filter(taxGroup %in% utaxatype_choose)

### run models with training data
cu1 <- lm(logvalue ~ logM, data=bdata_train)
summary(cu1)

cu2 <- lm(logvalue ~ region2, data=bdata_train)
summary(cu2)

cu3 <- lm(logvalue ~ FisheryType2, data=bdata_train)
summary(cu3)

cu4 <- lm(logvalue ~ taxGroup2, data=bdata_train)
summary(cu4)

cu5 <- lm(logvalue ~ region2 + logM, data=bdata_train)
summary(cu5)

cu6 <- lm(logvalue ~ FisheryType2 + logM, data=bdata_train)
summary(cu6)

cu7 <- lm(logvalue ~ taxGroup2 + logM, data=bdata_train)
summary(cu7)

cu8 <- lm(logvalue ~ region2 + FisheryType2, data=bdata_train)
summary(cu8)

cu9 <- lm(logvalue ~ region2 + taxGroup2, data=bdata_train)
summary(cu9)

cu10 <- lm(logvalue ~ logM + region2 + FisheryType2, data=bdata_train)
summary(cu10)

cu11 <- lm(logvalue ~ logM + region2 + taxGroup2, data=bdata_train)
summary(cu11)


cumods <- list("M"=cu1, "region"=cu2, "fishtype"=cu3, "taxagroup"=cu4,
				"region_M"=cu5, "fishtype_M"=cu6, "taxagroup_M"=cu7, 
				"region_fishtype"=cu8, "region_taxagroup"=cu9,
				"all_fishtype"=cu10, "all_taxagroup"=cu11)
cur2 <- sapply(1:length(cumods), function(x) summary(cumods[[x]])$adj.r.squared)

aicuc <- AIC(cu1, cu2, cu3, cu4, cu5, cu6, cu7, cu8, cu9, cu10, cu11) %>%
		dplyr::mutate(deltaAIC = AIC - min(AIC))
aicuc1 <- AIC(cu1, cu2, cu3, cu4) %>%
		dplyr::mutate(deltaAIC = AIC - min(AIC))
aicuc$R2 <- cur2
aicuc$ModelName <- names(cumods)
write.csv(aicu, file.path(fig_dir, "cross_BBmsy_AIC.csv"))

png(file.path(fig_dir, "cross_val_B.png"), height=8, width=10, res=200, units="in")
pred <- predict(cb10, bdata_test, se.fit=TRUE)
val <- exp(pred$fit)
low <- exp(pred$fit - 1.96 * pred$se.fit)
up <- exp(pred$fit + 1.96 * pred$se.fit)
plot(val, col="black", pch=19, ylim=c(0, max(up)), xlab="Stock", ylab="B/Bmsy")
segments(x0=1:length(val), x1=1:length(val), y0=low, y1=up, col="black")
obs <- bdata_test$value
within <- sapply(1:length(obs), function(x){
	ifelse(obs[x] >= low[x] & obs[x] <= up[x], 19, 1)
})
col <- sapply(1:length(obs), function(x){
	ifelse(obs[x] >= 1 & val[x] >=1, "blue", ifelse(obs[x] < 1 & val[x] < 1, "blue", "red"))
})
points(bdata_test$val, pch=within, col=col)
length(which(within==19))/length(within)
length(which(col=="red"))/length(col)
length(which(within==1 & col=="red"))/length(within)
dev.off()

png(file.path(fig_dir, "cross_val_U.png"), height=8, width=10, res=200, units="in")
pred <- predict(cu10, udata_test, se.fit=TRUE)
val <- exp(pred$fit)
low <- exp(pred$fit - 1.96 * pred$se.fit)
up <- exp(pred$fit + 1.96 * pred$se.fit)
plot(val, col="black", pch=19, ylim=c(0, max(up)), xlab="Stock", ylab="B/Bmsy")
segments(x0=1:length(val), x1=1:length(val), y0=low, y1=up, col="black")
obs <- udata_test$value
within <- sapply(1:length(obs), function(x){
	ifelse(obs[x] >= low[x] & obs[x] <= up[x], 19, 1)
})
col <- sapply(1:length(obs), function(x){
	ifelse(obs[x] >= 1 & val[x] >=1, "blue", ifelse(obs[x] < 1 & val[x] < 1, "blue", "red"))
})
points(bdata_test$val, pch=within, col=col)
length(which(within==19))/length(within)
length(which(col=="red"))/length(col)
length(which(within==1 & col=="red"))/length(within)
dev.off()
