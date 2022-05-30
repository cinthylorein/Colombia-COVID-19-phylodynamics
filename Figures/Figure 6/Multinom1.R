# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN Colombia (GISAID GENOMIC EPIDEMIOLOGY METADATA)
# Ricardo Rivero (Based in T. Wenseelers analysis)


library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE)
library(emmeans)
library(readr)
library(ggplot2)
library(ggthemes)
library(scales)
library(dplyr)
library(gam)
library(extrafont)

today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2020-08-10")
today_num = as.numeric(today)
plotdir = "Colombia_GISAID"
suppressWarnings(dir.create(paste0("~/Documents/R scripts/Phylodynamics",plotdir)))

# import GISAID genomic epidemiology metadata
GISAID = read.csv("~/Documents/R scripts/Phylodynamics/Multinomial fit databases/1st/GISAID1_upd.csv", sep = ",")
GISAID = as.data.frame(GISAID)

GISAID$date = as.Date(GISAID$date)
GISAID = GISAID[!is.na(GISAID$date),]
GISAID[GISAID$host!="Human","strain"]
GISAID = GISAID[GISAID$host=="Human",]
GISAID = GISAID[GISAID$date>=as.Date("2020-01-01"),]
range(GISAID$date) # "2020-03-06" "2020-08-26"
nrow(GISAID) # 386
GISAID$Week = lubridate::week(GISAID$date)
GISAID$Year = lubridate::year(GISAID$date)
GISAID$Year_Week = interaction(GISAID$Year,GISAID$Week)
library(lubridate)
GISAID$floor_date = as.Date(as.character(cut(GISAID$date, "week")))+3.5 # week midpoint date
GISAID$DATE_NUM = as.numeric(GISAID$date)
colnames(GISAID)

GISAID$pango_lineage[grepl("B.1", GISAID$pango_lineage)]
GISAID$pango_lineage[grepl("B.1.111",GISAID$pango_lineage)]
GISAID$pango_lineage[grepl("B.1.420",GISAID$pango_lineage)] 
GISAID$pango_lineage[grepl("B.1.1",GISAID$pango_lineage)]


sel_target_VOC = "B.1.420"
GISAID$LINEAGE2 = GISAID$pango_lineage


GISAID$country <- with(GISAID, ifelse(GISAID$host=="Human", "Colombia"))
table_country_lineage = as.data.frame(table(GISAID$country, GISAID$LINEAGE2))
colnames(table_country_lineage) = c("Country","Lineage","Count")
tblB1420 = table_country_lineage[grepl(sel_target_VOC, table_country_lineage$Lineage, fixed=T)&table_country_lineage$Count>10,]
tblB1420

sel_ref_lineage = "B.1"






# ANALYSIS OF VOCs IN Colombia ####

sel_countries = "Colombia"

tblB1 = table_country_lineage[table_country_lineage$Lineage==sel_ref_lineage&table_country_lineage$Count>10&table_country_lineage$Country %in% sel_countries,]
tblB1

GISAID_sel = GISAID[GISAID$country %in% sel_countries,]

# use data from March 06 2020 onwards
GISAID_sel = GISAID_sel[GISAID_sel$date>=as.Date("2020-03-06"),]
nrow(GISAID_sel) # 386
range(GISAID_sel$date) # "2020-03-06" "2020-08-26"

lin_count <- rowSums(table(GISAID_sel$LINEAGE2,GISAID_sel$country))

#Colombia metadata has no country_exposure column so me omit this step
GISAID_sel = GISAID_sel[!is.na(GISAID_sel$LINEAGE2),]
nrow(GISAID_sel) # 386


sum(GISAID_sel$LINEAGE2==sel_target_VOC) # 37
sum(GISAID_sel$LINEAGE2=="B.1.111") # 81
sum(GISAID_sel$LINEAGE2=="B.1.1") # 5

table(GISAID_sel$LINEAGE2)

main_lineages = names(table(GISAID_sel$LINEAGE2))[100*table(GISAID_sel$LINEAGE2)/sum(table(GISAID_sel$LINEAGE2)) > 5]
main_lineages
# "B.1.351"                  "Delta (B.1.617.2 & AY.X)"
VOCs = c("B.1", sel_target_VOC,"B.1.111","B.1.1")
main_lineages = union(main_lineages, VOCs)
GISAID_sel$LINEAGE2[!(GISAID_sel$LINEAGE2 %in% main_lineages)] = "other" # minority lineages & non-VOCs
remove2 = names(table(GISAID_sel$LINEAGE2))[table(GISAID_sel$LINEAGE2)/sum(table(GISAID_sel$LINEAGE2)) < 0.01]
remove2 = remove2[!(remove2 %in% c("B.1", sel_target_VOC,"B.1.111","B.1.1"))]
GISAID_sel$LINEAGE2[(GISAID_sel$LINEAGE2 %in% remove2)] = "other" # minority VOCs
GISAID_sel$LINEAGE2 = factor(GISAID_sel$LINEAGE2)
GISAID_sel$LINEAGE2 = relevel(GISAID_sel$LINEAGE2, ref="B.1") # we code UK strain as the reference level
levels(GISAID_sel$LINEAGE2)
# [1] "B.1.1.7"                  "B.1.1.519"                "B.1.351"                  "B.1.617.1"                "C.1.2"                   
# [6] "Delta (B.1.617.2 & AY.X)" "other" 
levels_LINEAGE2 = c("B.1","B.1.111",
                    "B.1.1", "other", sel_target_VOC)
GISAID_sel$LINEAGE2 = factor(GISAID_sel$LINEAGE2, levels=levels_LINEAGE2)

# GISAID_sel = GISAID_sel[GISAID_sel$division!="India",]
table(GISAID_sel$LINEAGE2)

range(GISAID_sel$date) # "2020-03-06" "2022-01-21"

# aggregated data to make Muller plots of raw data
# aggregate by day to identify days on which INSA performed (days with a lot of sequences)
# we subset the data to just those days to avoid sampling biases (delta infection clusters etc)
data_agbyday2 = as.data.frame(table(GISAID_sel$date, GISAID_sel$LINEAGE2))
colnames(data_agbyday2) = c("date", "LINEAGE2", "count")
data_agbyday2_sum = aggregate(count ~ date, data=data_agbyday2, sum)
data_agbyday2$total = data_agbyday2_sum$count[match(data_agbyday2$date, data_agbyday2_sum$date)]
sum(data_agbyday2[data_agbyday2$LINEAGE2=="B.1","total"]) == nrow(GISAID_sel) # correct
data_agbyday2$date = as.Date(as.character(data_agbyday2$date))
data_agbyday2$LINEAGE2 = factor(data_agbyday2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyday2$date_num = as.numeric(data_agbyday2$date)
data_agbyday2$prop = data_agbyday2$count/data_agbyday2$total
data_agbyday2$floor_date = NULL
qplot(data=data_agbyday2, x=date, y=total, colour=total>20, fill=total>20, geom="col")
GISAID_sel$total_sequenced_on_that_day = data_agbyday2$total[match(GISAID_sel$date, data_agbyday2$date)]
# GISAID_sel = GISAID_sel[GISAID_sel$total_sequenced_on_that_day>20,] # dates on which was performed
# nrow(GISAID_sel) # 5710

# aggregated by week
data_agbyweek2 = as.data.frame(table(GISAID_sel$floor_date, GISAID_sel$LINEAGE2))
colnames(data_agbyweek2) = c("floor_date", "LINEAGE2", "count")
data_agbyweek2_sum = aggregate(count ~ floor_date, data=data_agbyweek2, sum)
data_agbyweek2$total = data_agbyweek2_sum$count[match(data_agbyweek2$floor_date, data_agbyweek2_sum$floor_date)]
sum(data_agbyweek2[data_agbyweek2$LINEAGE2=="B.1","total"]) == nrow(GISAID_sel) # correct
data_agbyweek2$collection_date = as.Date(as.character(data_agbyweek2$floor_date))
data_agbyweek2$LINEAGE2 = factor(data_agbyweek2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweek2$collection_date_num = as.numeric(data_agbyweek2$collection_date)
data_agbyweek2$prop = data_agbyweek2$count/data_agbyweek2$total
data_agbyweek2$floor_date = NULL


# MULLER PLOT (RAW DATA)
library(scales)
n2 = length(levels(GISAID_sel$LINEAGE2))
lineage_cols2 = hcl(h = seq(15, 320, length = n2), l = 65, c = 200)
lineage_cols2[which(levels(GISAID_sel$LINEAGE2)=="B.1")] = "#91D1C2FF"
lineage_cols2[which(levels(GISAID_sel$LINEAGE2)=="B.1.111")] = "#64cc80"
lineage_cols2[which(levels(GISAID_sel$LINEAGE2)=="B.1.1")] = "#314f39"
lineage_cols2[which(levels(GISAID_sel$LINEAGE2)==sel_target_VOC)] = "magenta"
lineage_cols2[which(levels(GISAID_sel$LINEAGE2)=="Other")] = "grey75"

# overall
muller_colombia_raw = ggplot(data=data_agbyweek2, aes(x=collection_date, y=count, group=LINEAGE2)) +
  # facet_wrap(~ division) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01"))),1,1),
                     limits=as.Date(c("2020-03-01",NA)), 
                     expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_minimal() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN COLOMBIA FROM MARCH TO SEPTEMBER 2020") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_sel$Week)))
muller_colombia_raw

ggsave(filename = "Muller_prev1.jpg", plot = muller_colombia_raw, height = 20, width = 25, units = "cm", dpi = 900)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\south africa_muller plots_raw data.pdf"), width=8, height=6)



# multinomial fits
data_agbyweek2$LINEAGE2 = relevel(data_agbyweek2$LINEAGE2, ref=sel_target_VOC)
data_agbyweek2$DATE_NUM = as.numeric(data_agbyweek2$collection_date)


# model without province/division included as an extra factor
library(nnet)
library(splines)
set.seed(1)
fit1_colombia_multi = nnet::multinom(LINEAGE2 ~ scale(DATE_NUM), weights=count, data=data_agbyweek2, maxit=5000)
fit2_colombia_multi = nnet::multinom(LINEAGE2 ~ ns(DATE_NUM, df=2), weights=count, data=data_agbyweek2, maxit=5000)
fit3_colombia_multi = nnet::multinom(LINEAGE2 ~ ns(DATE_NUM, df=3), weights=count, data=data_agbyweek2, maxit=5000)
fit4_colombia_multi = nnet::multinom(LINEAGE2 ~ ns(DATE_NUM, df=4), weights=count, data=data_agbyweek2, maxit=5000)
BIC(fit1_colombia_multi, fit2_colombia_multi, fit3_colombia_multi, fit4_colombia_multi) 
# df      BIC
# fit1_colombia_multi  8 1072.250
# fit2_colombia_multi 12 1081.742
# fit3_colombia_multi 16 1077.680
# fit2_colombia_multi 20 1095.907


# growth rate advantage compared to Delta in overall model (difference in growth rate per day) 
emtrcolombia = emtrends(fit1_colombia_multi, trt.vs.ctrl ~ LINEAGE2,  
                        var="DATE_NUM",  mode="latent",
                        at=list(DATE_NUM=max(GISAID_sel$DATE_NUM)),
                        adjust="none", df=NA)
delta_r_colombia = data.frame(confint(emtrcolombia, 
                                      adjust="none", df=NA)$contrasts, 
                              p.value=as.data.frame(emtrcolombia$contrasts, adjust="none", df=NA)$p.value)
r_table <- delta_r_colombia
#                               contrast     estimate          SE df   asymp.LCL   asymp.UCL      p.value
# 1   B.1.1.7 - Delta (B.1.617.2 & AY.X) -0.083053146 0.007547032 NA -0.09784506 -0.06826124 3.625488e-28
# 2 B.1.1.519 - Delta (B.1.617.2 & AY.X) -0.585953574 0.646813186 NA -1.85368412  0.68177697 3.649843e-01
# 3   B.1.351 - Delta (B.1.617.2 & AY.X) -0.101157067 0.005956704 NA -0.11283199 -0.08948214 1.115155e-64
# 4 B.1.617.1 - Delta (B.1.617.2 & AY.X) -0.108863423 0.041498692 NA -0.19019937 -0.02752748 8.708315e-03
# 5     C.1.2 - Delta (B.1.617.2 & AY.X) -0.008334439 0.011326303 NA -0.03053358  0.01386471 4.618231e-01
# 6     other - Delta (B.1.617.2 & AY.X) -0.053854857 0.005526241 NA -0.06468609 -0.04302362 1.932153e-22


# fitted prop of different LINEAGES in Colombiafor 2020-09-01
multinom_preds_today_avg = data.frame(emmeans(fit1_colombia_multi, ~ LINEAGE2|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_pred_table <- multinom_preds_today_avg
# LINEAGE2       prob         SE df    asymp.LCL  asymp.UCL
# 1  B.1.420 0.40741167 0.14761295 NA  0.118095603 0.69672774
# 2      B.1 0.24463829 0.10290329 NA  0.042951552 0.44632504
# 3  B.1.111 0.07763974 0.04748541 NA -0.015429947 0.17070943
# 4    B.1.1 0.02710419 0.03186609 NA -0.035352208 0.08956058
# 5    other 0.24320611 0.12238535 NA  0.003335238 0.48307698

# % non-B.1.421
colSums(multinom_preds_today_avg[-1, c("prob","asymp.LCL","asymp.UCL")])
# prob    asymp.LCL    asymp.UCL 
# 0.592588331 -0.004495366  1.189672027 


# PLOT MULTINOMIAL FIT

# extrapolate = 30
date.from = as.numeric(as.Date("2020-03-01"))
date.to = as.numeric(as.Date("2020-10-01")) # max(GISAID_sel$DATE_NUM)+extrapolate

#multinomial model predictions (fastest, but no confidence intervals)
predgrid = expand.grid(list(DATE_NUM=seq(date.from, date.to)))
fit_colombia_multi_preds = data.frame(predgrid, as.data.frame(predict(fit1_colombia_multi, newdata=predgrid, type="prob")),check.names=F)
library(tidyr)
library(tidyselect)
fit_colombia_multi_preds = gather(fit_colombia_multi_preds, LINEAGE2, prob, all_of(levels_LINEAGE2), factor_key=TRUE)
fit_colombia_multi_preds$collection_date = as.Date(fit_colombia_multi_preds$DATE_NUM, origin="1970-01-01")
fit_colombia_multi_preds$LINEAGE2 = factor(fit_colombia_multi_preds$LINEAGE2, levels=levels_LINEAGE2) 
# 
# 
# 
muller_colombia_mfit = ggplot(data=fit_colombia_multi_preds,
                              aes(x=collection_date, y=prob, group=LINEAGE2)) +
  #   # facet_wrap(~ STATE) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(GISAID_sel$DATE_NUM)+1,
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01"))),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_minimal() + theme(legend.position="right",
                          axis.title.x=element_blank()) +
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN COLOMBIA FROM MARCH TO SEPTEMBER 2020\n(GISAID data, multinomial fit)")
muller_colombia_mfit

ggsave("muller_plot_pred1.jpg", plot = muller_colombia_mfit, height = 20, width = 25, units = "cm", dpi = 900)

#
#
library(ggpubr)
ggarrange(muller_colombia_raw + coord_cartesian(xlim=c(as.Date("2020-03-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))),
          muller_colombia_mfit+ggtitle("Multinomial fit"), ncol=1)

# ggsave(file=paste0(".\\plots\\",plotdir,"\\south africa_muller plots multipanel_multinom fit.png"), width=10, height=10)
# # ggsave(file=paste0(".\\plots\\",plotdir,"\\south africa_muller plots multipanel_multinom fit.pdf"), width=10, height=10)





# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS

# multinomial model predictions by division/province with confidence intervals (but slower)


# fit_southafrica_multi_preds2 = fit_southafrica_multi_preds_withCI

fit_colombia_multi_preds_withCI_overall = data.frame(emmeans(fit1_colombia_multi,
                                                             ~ LINEAGE2,
                                                             by=c("DATE_NUM"),
                                                             at=list(DATE_NUM=seq(date.from, date.to, by=1)),  # by=XX to speed up things a bit
                                                             mode="prob", df=NA))
fit_colombia_multi_preds_withCI_overall$collection_date = as.Date(fit_colombia_multi_preds_withCI_overall$DATE_NUM, origin="1970-01-01")
fit_colombia_multi_preds_withCI_overall$LINEAGE2 = factor(fit_colombia_multi_preds_withCI_overall$LINEAGE2, levels=levels_LINEAGE2)
fit_colombia_multi_preds2 = fit_colombia_multi_preds_withCI_overall

muller_colombia_mfit2 = ggplot(data=fit_colombia_multi_preds_withCI_overall, 
                               aes(x=collection_date, y=prob, group=LINEAGE2)) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(GISAID_sel$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01"))),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_minimal() + theme(legend.position="right", 
                          axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN COLOMBIA (GISAID DATA), multinomial fit)")
muller_colombia_mfit2

ggsave("muller_pred_smoothed1.jpg", plot = muller_colombia_mfit2, width = 25, height = 20, units = "cm", dpi = 900)


library(ggpubr)
ggarrange(muller_colombia_raw + coord_cartesian(xlim=c(as.Date("2020-03-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_colombia_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\colombia_muller plots multipanel_multinom fit_by province.png"), width=10, height=10)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\south africa_muller plots multipanel_multinom fit.pdf"), width=10, height=10)


muller_colombia_mfit_overall = ggplot(data=fit_colombia_multi_preds2, 
                                      aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  # facet_wrap(~ division) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(GISAID_sel$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01"))),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_minimal() + theme(legend.position="right", 
                          axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN COLOMBIA (GISAID data), multinomial fit")
muller_colombia_mfit_overall

ggsave(file=paste0(".\\plots\\",plotdir,"\\colombia_muller plots_multinom fit_overall.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\south africa_muller plots_multinom fit_overall.pdf"), width=10, height=6)


library(ggpubr)
ggarrange(muller_colombia_raw + coord_cartesian(xlim=c(as.Date("2020-03-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_colombia_mfit_overall+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\colombia_muller plots multipanel_multinom fit_overall.png"), width=10, height=10)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\south africa_muller plots multipanel_multinom fit.pdf"), width=10, height=10)






# on logit scale:

ymin = 0.001
ymax = 0.999
fit_colombia_multi_preds2$asymp.LCL[fit_colombia_multi_preds2$asymp.LCL<ymin] = ymin
fit_colombia_multi_preds2$asymp.UCL[fit_colombia_multi_preds2$asymp.UCL<ymin] = ymin
fit_colombia_multi_preds2$asymp.UCL[fit_colombia_multi_preds2$asymp.UCL>ymax] = ymax
fit_colombia_multi_preds2$prob[fit_colombia_multi_preds2$prob<ymin] = ymin

plot_colombia_mfit_logit = qplot(data=fit_colombia_multi_preds2, x=collection_date, y=prob, geom="blank") +
  # facet_wrap(~ STATE) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_minimal() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN COLOMBIA (GISAID data), multinomial fit") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01"))),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  geom_point(data=data_agbyweek2,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE2
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.5, 4), limits=c(1,max(data_agbyweek2$total)), breaks=c(10,100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2020-03-01"),as.Date(date.to, origin="1970-01-01")), ylim=c(0.001, 0.9901), expand=c(0,0))
plot_colombia_mfit_logit

ggsave("mnomial_logit1.jpg", plot = plot_colombia_mfit_logit, height = 20, width = 25, units = "cm", dpi = 900)

# on response scale:
plot_colombia_mfit1 = qplot(data=fit_colombia_multi_preds2, x=collection_date, y=100*prob, geom="blank") +
  # facet_wrap(~ STATE) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_minimal() + xlab("") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01"))),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +
  scale_x_date(date_labels = "%b-%Y") +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2020-03-01",NA)),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  geom_point(data=data_agbyweek2,
             aes(x=collection_date, y=100*prop, size=total,
                 colour=LINEAGE2
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.5, 5), limits=c(1,max(data_agbyweek2$total)), breaks=c(10,100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right", 
        axis.title.x=element_text(vjust=-1, size=14, family="Times New Roman"),
        axis.title.y=element_text(vjust=-0.25, size=14, family="Times New Roman"), 
        legend.text=element_text(size=14, family="Times New Roman"), 
        legend.title=element_blank(),
        legend.key.height=unit(1, "line"), 
        legend.key.size=unit(0.8, "cm"), 
        legend.key=element_rect(fill=NA), 
        legend.background=element_blank(),
        axis.text.x = element_text(family = "Times New Roman"),
        plot.margin=unit(c(1,1,1,1), "cm")) +
  xlab("Collection date")
plot_colombia_mfit1

ggsave("mnomial_response1.svg", plot = plot_colombia_mfit1, height = 20, width = 25, units = "cm", dpi = 1200)
ggsave("mnomial_response1.jpg", plot = plot_colombia_mfit1, height = 20, width = 25, units = "cm", dpi = 1200)

# ggsave(file=paste0(".\\plots\\",plotdir,"\\south africa_multinom fit_response scale.pdf"), width=10, height=6)




# PLOTS OF NEW CASES PER DAY BY VARIANT & EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####

# load case data
library(covidregionaldata)
library(dplyr)
library(ggplot2)
library(scales)
library(gam)
library(mgcv)

cases_tot = as.data.frame(get_national_data(countries = "Colombia"))
cases_tot$date <- as.Date(cases_tot$date) 
cases_tot <- filter(cases_tot, between(date, as.Date("2020-03-01"), as.Date("2020-08-10")))
cases_tot = cases_tot[cases_tot$date>=as.Date("2020-03-06"),]
cases_tot$DATE_NUM = as.numeric(cases_tot$date)
# cases_tot$BANKHOLIDAY = bankholiday(cases_tot$date)
cases_tot$WEEKDAY = weekdays(cases_tot$date)
# cases_tot = cases_tot[cases_tot$date<=(max(cases_tot$date)-3),] # cut off data from last 3 days (incomplete)

# smooth out weekday effects in case nrs using GAM (if testing data is available one could correct for testing intensity as well)
fit_cases = gam(cases_new ~ s(DATE_NUM, bs="cs", k=25, fx= FALSE) + 
                  WEEKDAY, family=poisson(log), data=cases_tot) 
BIC(fit_cases) #198403.9

# STACKED AREA CHART OF NEW CASES BY VARIANT (MULTINOMIAL FIT MAPPED ONTO CASE DATA) ####

fit_colombia_multi_preds_withCI_overall$totcases = cases_tot$cases_new[match(round(fit_colombia_multi_preds_withCI_overall$DATE_NUM),
                                                                             cases_tot$DATE_NUM)]
fit_colombia_multi_preds_withCI_overall$cases = fit_colombia_multi_preds_withCI_overall$totcases * fit_colombia_multi_preds_withCI_overall$prob
fit_colombia_multi_preds_withCI_overall$cases[fit_colombia_multi_preds_withCI_overall$cases<=0.001] = NA
cases_emmeans = as.data.frame(emmeans(fit_cases, ~ DATE_NUM, at=list(DATE_NUM=seq(date.from, date.to, by=0.5), BANHOLIDAY="no"), type="response"))
fit_colombia_multi_preds_withCI_overall$smoothed_totcases = 
  cases_emmeans$rate[match(fit_colombia_multi_preds_withCI_overall$DATE_NUM,
                           cases_emmeans$DATE_NUM)]
fit_colombia_multi_preds_withCI_overall$smoothed_cases = fit_colombia_multi_preds_withCI_overall$smoothed_totcases * fit_colombia_multi_preds_withCI_overall$prob
fit_colombia_multi_preds_withCI_overall$smoothed_cases[fit_colombia_multi_preds_withCI_overall$smoothed_cases<=0.001] = NA

fitted_cases_raw <- ggplot(data=fit_colombia_multi_preds_withCI_overall[fit_colombia_multi_preds_withCI_overall$collection_date>=as.Date("2020-03-01"),], 
                           aes(x=collection_date, y=cases, group=LINEAGE2)) + 
  # facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01"))),1,1),
                     # limits=c(as.Date("2021-03-01"),max(cases_tot$date)), 
                     expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_minimal() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN COLOMBIA (case data & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  coord_cartesian(xlim=c(as.Date("2020-03-01"),NA))

ggsave("fitted_cases_raw1.jpg", plot = fitted_cases_raw, height = 20, width = 25, units = "cm", dpi = 900)

fitted_cases_smoothed <- ggplot(data=fit_colombia_multi_preds_withCI_overall[fit_colombia_multi_preds_withCI_overall$collection_date>=as.Date("2020-01-01")&
                                                                               fit_colombia_multi_preds_withCI_overall$collection_date<=max(cases_tot$date),], 
                                aes(x=collection_date-7, y=smoothed_cases, group=LINEAGE2)) + 
  # facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01"))),1,1),
                     # limits=c(as.Date("2021-03-01"),today), 
                     expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_minimal() + theme(legend.position="right") + 
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of infection") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN COLOMBIA (case data & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  coord_cartesian(xlim=c(as.Date("2020-03-01"),max(cases_tot$date)))

ggsave("fitted_cases_smoothed1.svg", plot = fitted_cases_smoothed, height = 20, width = 25, units = "cm", dpi = 900)



# EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####

# Function to calculate Re values from intrinsic growth rate
# cf. https://github.com/epiforecasts/EpiNow2/blob/5015e75f7048c2580b2ebe83e46d63124d014861/R/utilities.R#L109
# https://royalsocietypublishing.org/doi/10.1098/rsif.2020.0144
# (assuming gamma distributed gen time)
Re.from.r <- function(r, gamma_mean=4.7, gamma_sd=2.9) { # Nishiura et al. 2020, or use values from Ferretti et al. 2020 (gamma_mean=5.5, gamma_sd=1.8)
  k <- (gamma_sd / gamma_mean)^2
  R <- (1 + k * r * gamma_mean)^(1 / k)
  return(R)
}


# calculate average instantaneous growth rates & 95% CLs using emtrends ####
# based on the slope of the GAM fit on a log link scale
avg_r_cases = as.data.frame(emtrends(fit_cases, ~ DATE_NUM, var="DATE_NUM", 
                                     at=list(DATE_NUM=seq(date.from,
                                                          date.to)#,
                                             # BANKHOLIDAY="no"
                                     ), # weekday="Wednesday",
                                     type="link"))
colnames(avg_r_cases)[2] = "r"
colnames(avg_r_cases)[5] = "r_LOWER"
colnames(avg_r_cases)[6] = "r_UPPER"
avg_r_cases$DATE = as.Date(avg_r_cases$DATE_NUM, origin="1970-01-01") # -7 TO CALCULATE BACK TO INFECTION DATE
avg_r_cases$Re = Re.from.r(avg_r_cases$r)
avg_r_cases$Re_LOWER = Re.from.r(avg_r_cases$r_LOWER)
avg_r_cases$Re_UPPER = Re.from.r(avg_r_cases$r_UPPER)
avg_r_cases = avg_r_cases[complete.cases(avg_r_cases),]
qplot(data=avg_r_cases, x=DATE-7, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", alpha=I(0.5), fill=I("steelblue")) +
  # facet_wrap(~ REGION) +
  geom_line() + theme_hc() + xlab("Date of infection") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01", "2020-10-01")),
                     labels=c("Mar 2020","Apr 2020","May 2020","Jun 2020","Jul 2020","Aug 2020","Sep 2020", "Oct 2020")) +
  # scale_y_continuous(limits=c(1/2, 2), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re IN COLOMBIA AT MOMENT OF INFECTION BASED ON NEW CASES") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
# coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

# calculate above-average intrinsic growth rates per day of each variant over time based on multinomial fit using emtrends weighted effect contrasts ####
# for best model fit3_sanger_multi
above_avg_r_variants0 = do.call(rbind, lapply(seq(date.from,
                                                  date.to), 
                                              function (d) { 
                                                wt = as.data.frame(emmeans(fit1_colombia_multi, ~ LINEAGE2 , at=list(DATE_NUM=d), type="response"))$prob   # important: these should sum to 1
                                                # wt = rep(1/length(levels_VARIANTS), length(levels_VARIANTS)) # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_LINEAGE2)
                                                cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
                                                names(cons) = seq_along(cons)
                                                EMT = emtrends(fit3_colombia_multi,  ~ LINEAGE2 , by=c("DATE_NUM"),
                                                               var="DATE_NUM", mode="latent",
                                                               at=list(DATE_NUM=d))
                                                out = as.data.frame(confint(contrast(EMT, cons), adjust="none", df=NA))
                                                # sum(out$estimate*wt) # should sum to zero
                                                return(out) } ))
above_avg_r_variants = above_avg_r_variants0
above_avg_r_variants$contrast = factor(above_avg_r_variants$contrast, 
                                       levels=1:length(levels_LINEAGE2), 
                                       labels=levels_LINEAGE2)
above_avg_r_variants$variant = above_avg_r_variants$contrast # gsub(" effect|\\(|\\)","",above_avg_r_variants$contrast)
above_avg_r_variants$collection_date = as.Date(above_avg_r_variants$DATE_NUM, origin="1970-01-01")
range(above_avg_r_variants$collection_date) # "2021-01-04" "2021-07-30"
above_avg_r_variants$avg_r = avg_r_cases$r[match(above_avg_r_variants$collection_date,
                                                 avg_r_cases$DATE)]  # average growth rate of all lineages calculated from case nrs
above_avg_r_variants$r = above_avg_r_variants$avg_r+above_avg_r_variants$estimate
above_avg_r_variants$r_LOWER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.LCL
above_avg_r_variants$r_UPPER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.UCL
above_avg_r_variants$Re = Re.from.r(above_avg_r_variants$r)
above_avg_r_variants$Re_LOWER = Re.from.r(above_avg_r_variants$r_LOWER)
above_avg_r_variants$Re_UPPER = Re.from.r(above_avg_r_variants$r_UPPER)
df = data.frame(contrast=NA,
                DATE_NUM=avg_r_cases$DATE_NUM, # -7 to calculate back to time of infection
                # REGION=avg_r_cases$REGION,
                estimate=NA,
                SE=NA,
                df=NA,
                asymp.LCL=NA,
                asymp.UCL=NA,
                # p.value=NA,
                collection_date=avg_r_cases$DATE,
                variant="avg",
                avg_r=avg_r_cases$r,
                r=avg_r_cases$r,
                r_LOWER=avg_r_cases$r_LOWER,
                r_UPPER=avg_r_cases$r_UPPER,
                Re=avg_r_cases$Re,
                Re_LOWER=avg_r_cases$Re_LOWER,
                Re_UPPER=avg_r_cases$Re_UPPER)
# df = df[df$DATE_NUM<=max(above_avg_r_variants$DATE_NUM)&df$DATE_NUM>=(min(above_avg_r_variants$DATE_NUM)+7),]
above_avg_r_variants = rbind(above_avg_r_variants, df)
above_avg_r_variants$variant = factor(above_avg_r_variants$variant, levels=c(levels_LINEAGE2,"avg"))
above_avg_r_variants$prob = fit_colombia_multi_preds_withCI_overall$prob[match(interaction(above_avg_r_variants$DATE_NUM,
                                                                                           above_avg_r_variants$variant),
                                                                               interaction(fit_colombia_multi_preds_withCI_overall$DATE_NUM,
                                                                                           fit_colombia_multi_preds_withCI_overall$LINEAGE2))]
above_avg_r_variants2 = above_avg_r_variants
ymax = 4
ymin = 1/3
above_avg_r_variants2$Re[above_avg_r_variants2$Re>=ymax] = NA
above_avg_r_variants2$Re[above_avg_r_variants2$Re<=ymin] = NA
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$Re_LOWER>=ymax] = ymax
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$Re_LOWER<=ymin] = ymin
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$Re_UPPER>=ymax] = ymax
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$Re_UPPER<=ymin] = ymin
above_avg_r_variants2$Re[above_avg_r_variants2$prob<0.01] = NA
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$prob<0.01] = NA
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$prob<0.01] = NA
r_variants1 <- qplot(data=above_avg_r_variants2[!((above_avg_r_variants2$variant %in% c("other"))|above_avg_r_variants2$collection_date>max(cases_tot$DATE)),], 
                    x=collection_date-7, # -7 to calculate back to date of infection
                    y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", colour=variant, fill=variant, alpha=I(0.5),
                    group=variant, linetype=I(0)) +
  # geom_ribbon(aes(fill=variant, colour=variant), alpha=I(0.5))
  geom_line(aes(colour=variant), lwd=I(0.72)) + theme_minimal() + xlab("Date of infection") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01")),
                     labels=c("Mar 2020","Apr 2020","May 2020","Jun 2020","Jul 2020","Aug 2020")) +
  scale_y_continuous(limits=c(1/ymax,ymax), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  # labs(tag = tag) +
  theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8),
        axis.title.x=element_text(vjust=-1, size=14, family="Times New Roman"),
        axis.title.y=element_text(vjust=-0.25, size=14, family="Times New Roman"),
        axis.text.x = element_text(family = "Times New Roman"),
        legend.text=element_text(size=14, family="Times New Roman"), legend.title=element_blank(),
        legend.key.height=unit(1, "line"), 
        legend.key.size=unit(0.8, "cm"), 
        legend.key=element_rect(fill=NA), 
        legend.background=element_blank()) +
  coord_cartesian(xlim=c(as.Date("2020-03-01"),max(cases_tot$date))) +
  scale_fill_manual("variant", values=c(head(lineage_cols2,-1),"black")) +
  scale_colour_manual("variant", values=c(head(lineage_cols2,-1),"black")) +
  theme(legend.position="right") 

r_variants1

ggsave("R_variants1.svg", plot = r_variants1, width = 25, height = 20, units = "cm", dpi = 1200)
ggsave("R_variants1.jpg", plot = r_variants1, width = 25, height = 20, units = "cm", dpi = 1200)

#save tables with gt
library(gt)

r_gt <- gt(r_table) %>%
  tab_header(title = md("**Relative growth advantage per variant**"))

multinom_gt <- gt(multinom_pred_table) %>%
  tab_header(title = md("**Multinomial Model prediction on Growth Advantage**"))

rt_gt <- gt(above_avg_r_variants2) %>%
  tab_header(title = md("**Relative Rt of variants**"))


#pairwise comparison

emtrcol_pairw = emtrends(fit1_colombia_multi, pairwise ~ LINEAGE2,  
                         var="DATE_NUM",  mode="latent",
                         at=list(DATE_NUM=max(GISAID$date)))
delta_r_colombia_pairw = data.frame(confint(emtrcol_pairw, 
                                            adjust="none", df=NA)$contrasts, 
                                    p.value=as.data.frame(emtrcol_pairw$contrasts)$p.value)

delta_r_colombia_pairw

pairwise_r <- gt(delta_r_colombia_pairw) %>%
  tab_header(title = md("**Pairwise relative R between variants**"))

pairwise_r

gtsave(r_gt, filename = "r_table1.png", expand = 10, path = "~/Documents/R scripts/Phylodynamics/Multinomial fit databases")
gtsave(multinom_gt, filename = "mnom_table1.pdf", expand = 10, path = "~/Documents/R scripts/Phylodynamics/Multinomial fit databases")
gtsave(pairwise_r, filename = "pwise_r_table1.pdf", expand = 10, path = "~/Documents/R scripts/Phylodynamics/Multinomial fit databases")
gtsave(rt_gt, filename = "rt_table1.pdf", expand = 10, path = "~/Documents/R scripts/Phylodynamics/Multinomial fit databases")

