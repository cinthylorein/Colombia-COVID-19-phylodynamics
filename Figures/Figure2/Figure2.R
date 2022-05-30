library(ggplot2)
library(data.table)
library(ggpubr)
library(dplyr)
library(maps)
library(viridis)
library(ggmap)
library(raster)
library(sp)
library(ggsci)
library(readxl)
library(lubridate)
library(plotly)
library(htmlwidgets)


setwd("C:/Users/Nikla/Dropbox/SARS-CoV-2_Colombia_Phylodinamics-/Database/rawdataset/Colombia")


colombia <- as.data.frame(fread("Colombia_allseq.tsv"))
#To know the total of lineages
#length(unique(colombia$lineage))

#colombia$date <- as.numeric(colombia$date)
colombia$date <- as.Date(colombia$date)


colombia$date <- format(colombia$date, "%Y-%m")

colombia <- colombia %>% distinct(strain, .keep_all = TRUE)

totalgenomes <- as.data.frame(table(colombia$division))
#cundinamarca = 261+ 3217
totalgenomes[totalgenomes == 261] <- 3478

totalgenomes <- totalgenomes %>% filter(Var1 != "Bogota")
totalgenomes <- totalgenomes %>% filter(Var1 != "Na")
id <- c(1,12,23,
        27,28,29,
        30,31,32,
        2,3,4,
        5,6,7,
        10,8,9,
        11,13,14,
        15,16,17,
        18,19,20,
        21,22,24,
        25,26)
totalgenomes <- cbind(totalgenomes,id)

names(totalgenomes) <- c("states","Genomes","id")

##############################################################
#plot all lineages vs VOC
ggplot(colombia, aes(lineage_b, fill=VOC)) +
  geom_bar() +
  coord_flip()+
  theme_minimal()

#############################################
#plot VOC through time 

lin_share <- as.data.frame(table(colombia$date,colombia$VOC))
names(lin_share) <- c("date","VOC","Count")


data <- lin_share %>%
  group_by(date, VOC) %>%
  summarise(n = sum(Count)) %>%
  mutate(Percentage = n/sum(n))


data_num_date <- data %>%
  mutate_at(vars(date), ym)
  

time <- ggplot(data = data_num_date, mapping = aes(x=date, y=Percentage, fill = VOC)) +
  stat_smooth(se=FALSE, geom="area",
              method = 'loess',
              span = 0.1,aes(fill=VOC), 
              position = "stack")+
  #geom_area(alpha=0.6, size=0.01)+
  #geom_histogram(position = "fill") +
  theme_minimal()+
  #scale_fill_brewer(palette="Spectral")+
  scale_color_npg()+
  scale_fill_npg()+
  theme(text = element_text(size=20))+
  xlab("Time") +
  ylab ("Relative Prevalence (%)")+ 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_labels = "%Y")+
  theme(legend.position = "none") 


###################################################3
#plot number of lineages per state
states <- ggplot(colombia, aes(division, fill = VOC)) +
  geom_bar() +
  coord_flip()+
  theme_minimal()+
  #scale_fill_brewer(palette="Spectral")+
  scale_color_npg()+
  scale_fill_npg()+
  theme(text = element_text(size=20))+
  xlab("Colombian States") +
  ylab ("No. genomes") 


############ map of colombian's states with number of available genomes


co <- getData("GADM", country = "CO", level = 1, download = TRUE)
col_depto <- fortify(co)  # make compatible to ggplot2

locat = as.vector(bbox(co))

ncmap = get_map(location = locat, source = "stamen", maptype = "toner", zoom = 6)
# ggmap(ncmap) not nice

class(col_depto$id)
totalgenomes$id <- as.character(totalgenomes$id)
col_depto <- left_join(col_depto, totalgenomes, by = "id")

col_depto2 <- col_depto %>% filter(id == "2")

mapbase <- ggplot(col_depto, aes(long, lat, group = group)) + 
  geom_polygon(aes(fill = Genomes), color = "blue") + 
  coord_equal() + 
 # geom_text(aes(label = id)) 
  scale_fill_gradient(high = "black", low = "white", guide = "colorbar") +
  geom_path(color = "grey")+
  theme(text = element_text(size=20))+
  labs(title = "Genomes per Colombian state")
  

### to see per state, you could use (id == "32") in the fill 

p_all<-ggarrange(time, states, mapbase, ncol = 3, nrow = 1)

setwd("C:/Users/Nikla/Dropbox/SARS-CoV-2_Colombia_Phylodinamics-/Figures")
ggsave("Figure2.jpg",plot=p_all,dpi=500,width=30,height=10)




##############333
table(colombia$division,colombia$VOC)

lin_share <- as.data.frame(table(colombia$division,colombia$lineage_c))
names(lin_share) <- c("state","lineage","Count")

colombia$date <- format(colombia$date, "%Y")


data <- lin_share %>%
  group_by(state, lineage) %>%
  summarise(n = sum(Count)) %>%
  mutate(Percentage = n/sum(n))


data_num_date <- data %>%
  mutate_at(vars(date), ym)


#####################################3

library(dplyr)
arrests <- USArrests 
arrests$region <- tolower(rownames(USArrests))
head(arrests)

# Retrieve the states map data and merge with crime data
states_map <- map_data("state")
arrests_map <- left_join(states_map, arrests, by = "region")

# Create the map
ggplot(arrests_map, aes(long, lat, group = group))+
  geom_polygon(aes(fill = Assault), color = "white")+
  scale_fill_viridis_c(option = "C")












##########################################################
world_map <- map_data("world")
ggplot(world_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill="lightgray", colour = "white")

colombia <- map_data("world", region = "Colombia")

colombia.lab.data <- colombia %>%
  group_by(region) %>%
  summarise(long = mean(long), lat = mean(lat))

ggplot(colombia, aes(x = long, y = lat)) +
  geom_polygon(aes( group = group, fill = region))+
  geom_text(aes(label = region), data = colombia.lab.data,  size = 3, hjust = 0.5)+
  scale_fill_viridis_d()+
  theme_void()+
  theme(legend.position = "none")



############################################


