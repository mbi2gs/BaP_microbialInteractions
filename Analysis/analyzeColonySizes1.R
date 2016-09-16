#-------------------------------------------------------
# Read in data extracted from images by CellProfiler
#-------------------------------------------------------

library(ggplot2)
library(reshape2)
library(dplyr)


#-------------------------------------------------------
# Read in data
#-------------------------------------------------------
colonyInfo = read.table('EditedObjects.csv',sep=',',header=T,row.names=NULL)
ci = select(colonyInfo,ImageNumber,ObjectNumber,Metadata_Media,Metadata_Left,Metadata_Right,AreaShape_Area,AreaShape_Center_X)
ci

#-------------------------------------------------------
# Helper functions
#-------------------------------------------------------

# Within an image, find the left (or right) colony 
keepLeftOrRight <- function(areadf,LoR)
{
  onlyLeftAreas = areadf[1,]
  for(i in unique(areadf$ImageNumber))
  {
    tmp = areadf[areadf$ImageNumber == i,]
    if(nrow(tmp) > 1)
    {
      tmp = arrange(tmp,AreaShape_Center_X)
      if(LoR == 'L')
      {
        onlyLeftAreas = rbind(onlyLeftAreas,tmp[1,])
      }
      else
      {
        onlyLeftAreas = rbind(onlyLeftAreas,tmp[2,])
      }
      
    }
    else
    {
      onlyLeftAreas = rbind(onlyLeftAreas,tmp)
    }    
  }
  onlyLeftAreas = onlyLeftAreas[-1,]
  return(onlyLeftAreas)
}

retrieveAreasForSpecies <- function(speciesName)
{
  areas_left = ci[ci$Metadata_Left == speciesName,]
  areas_left = keepLeftOrRight(areas_left,'L')
  areas_left = mutate(areas_left, description = paste0(Metadata_Media,'_',Metadata_Right))
  areas_left = mutate(areas_left, paired_species = Metadata_Right)
    
  areas_right = ci[ci$Metadata_Right == speciesName,]
  areas_right = keepLeftOrRight(areas_right,'R')
  areas_right = mutate(areas_right,description = paste0(Metadata_Media,'_',Metadata_Left))
  areas_right = mutate(areas_right, paired_species = Metadata_Left)
  
  areas = rbind(areas_left,areas_right)
  return(areas)
}

#-------------------------------------------------------
# Bar plots
#-------------------------------------------------------
## Get areas for standards
std_areas = retrieveAreasForSpecies('standard')
std_areas = select(std_areas,description,AreaShape_Area,Metadata_Media)
std_area = mean(std_areas$AreaShape_Area)
sqcm = .0625 # the standard area is 0.254cm * 0.254cm = 6.4516 sq. mm


species = c("p14","p1","s","hi","hp")
species_labels = c("P. aeruginosa PA14","P. aeruginosa PA01","S. aureus","H. influenzae","H. parainfluenzae")
yL = c(0.1,0.1,0,0,0)
yU = c(0.5,0.6,0.225,0.04,0.075)

for( i in 1:length(species))
{
  ## Get areas for species i and plot
  areas = retrieveAreasForSpecies(species[i])
  areas = select(areas,description,AreaShape_Area,Metadata_Media)
  areas$AreaShape_Area = areas$AreaShape_Area * sqcm / std_area
  
  levelsList = c(ordered=TRUE)
  for(j in 1:length(species))
  {
    if(i != j)
    {
      levelsList = c(levelsList,paste0("cntrl_",species[j]),paste0("benzo_",species[j]))
    }
    else
    {
      levelsList = c(levelsList,"cntrl_none","benzo_none")
    }
  }
  areas$description = factor(areas$description,levels=levelsList)
  
  # Plot
  p = ggplot(areas, aes(x=factor(description), y=AreaShape_Area)) + 
    geom_boxplot(outlier.size=-1,size=1,aes(color=Metadata_Media)) + 
    scale_colour_brewer(palette = "Set1") +
    geom_jitter(color='grey',alpha=0.7,size=2, width = 0.5) +
    ylim(yL[i],yU[i]) +
    ylab(expression(paste("Area (mm"^bold("2"),")"))) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1,size=14), axis.title.x = element_blank(),legend.position="none")
  print(p)
  
  ggsave(paste0("../Figures/",species_labels[i],".jpg"),width=8,height=2.5)
}


#-------------------------------------------------------
# Heat maps
#-------------------------------------------------------
library(pheatmap)
library(RColorBrewer)

# Fold change for heatmap for colony area fold change in control media 
# (fold change of column species grown next to row species relative to growth alone)
fold_change_control = matrix(data=NA,nrow=length(species),ncol=length(species))
colnames(fold_change_control) = species
row.names(fold_change_control) = species

p_values_control = matrix(data=NA,nrow=length(species),ncol=length(species))

for(i in species)
{
  speciesi_areas = retrieveAreasForSpecies(i)
  speciesi_areas = select(speciesi_areas,description,AreaShape_Area,Metadata_Media,paired_species)
  speciesi_areas$AreaShape_Area = speciesi_areas$AreaShape_Area * sqcm / std_area
  
  for(j in species)
  {
    if(i != j)
    {
      speciesi_none_cntrls = speciesi_areas[speciesi_areas$paired_species=='none' & speciesi_areas$Metadata_Media=='cntrl',]$AreaShape_Area
      speciesi_none_cntrl = mean(speciesi_none_cntrls)
      speciesi_paired_cntrls = speciesi_areas[speciesi_areas$paired_species==j & speciesi_areas$Metadata_Media=='cntrl',]$AreaShape_Area
      speciesi_paired_cntrl = mean(speciesi_paired_cntrls)
      p = wilcox.test(speciesi_none_cntrls,speciesi_paired_cntrls,alternative="two.sided")

      fold_change_control[row.names(fold_change_control) == j, colnames(fold_change_control) == i] = speciesi_paired_cntrl/speciesi_none_cntrl
      p_values_control[row.names(fold_change_control) == j, colnames(fold_change_control) == i] = p$p.value
    }
  }
}

adjustedPs = p.adjust(p_values_control, method = "fdr")
p_values_control_ajstd = matrix(adjustedPs,nrow=length(species),ncol=length(species))
colnames(p_values_control_ajstd) = species
row.names(p_values_control_ajstd) = species
write.table(p_values_control_ajstd, file="p_values_control_ajstd.txt", quote = FALSE, sep="\t")

#create the breaks
bk2 = unique(c(seq(min(fold_change_control, na.rm = TRUE), 0.99, length=9), 1, seq(1.01, max(fold_change_control, na.rm = TRUE), length=10)))

#set different color vectors for each interval
col1 = colorRampPalette(c("blue", 'white'))(9) #set the order of greys
col2 <- rep("white", 1)
col3 = colorRampPalette(c("white", "red"))(9)
colors2 <- c(col1, col2, col3)


#draw heatmap
fontsize = 12

hm.parameters <- list(fold_change_control, 
                      breaks=bk2,
                      color = colors2,
                      cellwidth = 20, cellheight = 20, scale = "none",
                      treeheight_row = 0,
                      treeheight_col = 0,
                      fontsize = fontsize, fontsize_row = fontsize,
                      fontsize_col = fontsize,
                      kmeans_k = NA,
                      Rowv=FALSE,na.rm=F, na.color="black",  Colv=FALSE, 
                      show_rownames = T, show_colnames = T,
                      main = NA,
                      cluster_rows = FALSE, cluster_cols = FALSE,
                      filename = "../Figures/fold_change_control_heat_map.png")

# To draw the heat map on screen 
do.call("pheatmap", hm.parameters)



# Fold change for heatmap for colony area fold change in benzo media 
# (fold change of column species grown next to row species relative to growth alone)
fold_change_benzo = matrix(data=NA,nrow=length(species),ncol=length(species))
colnames(fold_change_benzo) = species
row.names(fold_change_benzo) = species

p_values_benzo = matrix(data=NA,nrow=length(species),ncol=length(species))

for(i in species)
{
  speciesi_areas = retrieveAreasForSpecies(i)
  speciesi_areas = select(speciesi_areas,description,AreaShape_Area,Metadata_Media,paired_species)
  speciesi_areas$AreaShape_Area = speciesi_areas$AreaShape_Area * sqcm / std_area
  
  for(j in species)
  {
    if(i != j)
    {
      speciesi_none_benzos = speciesi_areas[speciesi_areas$paired_species=='none' & speciesi_areas$Metadata_Media=='benzo',]$AreaShape_Area
      speciesi_none_benzo = mean(speciesi_none_benzos)
      speciesi_paired_benzos = speciesi_areas[speciesi_areas$paired_species==j & speciesi_areas$Metadata_Media=='benzo',]$AreaShape_Area
      speciesi_paired_benzo = mean(speciesi_paired_benzos)
      p = wilcox.test(speciesi_none_benzos,speciesi_paired_benzos,alternative="two.sided")
      
      fold_change_benzo[row.names(fold_change_benzo) == j, colnames(fold_change_benzo) == i] = speciesi_paired_benzo/speciesi_none_benzo
      p_values_benzo[row.names(fold_change_benzo) == j, colnames(fold_change_benzo) == i] = p$p.value
    }
  }
}

adjustedPs = p.adjust(p_values_benzo, method = "fdr")
p_values_benzo_ajstd = matrix(adjustedPs,nrow=length(species),ncol=length(species))
colnames(p_values_benzo_ajstd) = species
row.names(p_values_benzo_ajstd) = species
write.table(p_values_benzo_ajstd, file="p_values_benzo_ajstd.txt", quote = FALSE, sep="\t")

#create the breaks
bk2 = unique(c(seq(min(fold_change_benzo, na.rm = TRUE), 0.99, length=9), 1, seq(1.01, max(fold_change_benzo, na.rm = TRUE), length=10)))

#draw heatmap
fontsize = 12

hm.parameters <- list(fold_change_benzo, 
                      breaks = bk2,
                      color = colors2,
                      cellwidth = 20, cellheight = 20, scale = "none",
                      fontsize = fontsize, fontsize_row = fontsize,
                      fontsize_col = fontsize,
                      kmeans_k = NA, 
                      show_rownames = T, show_colnames = T,
                      cluster_rows = FALSE, cluster_cols = FALSE,
                      filename = "../Figures/fold_change_benzo_heat_map.png")

# To draw the heat map on screen 
do.call("pheatmap", hm.parameters)


#heatmap of control vs benzo
fold_change = matrix(data=NA,nrow=5,ncol=5)
colnames(fold_change) = species
row.names(fold_change) = species

fold_change = fold_change_control-fold_change_benzo

#create the breaks
bk2 = unique(c(seq(min(fold_change, na.rm = TRUE), -0.001, length=9), 0, seq(0.001, max(fold_change, na.rm = TRUE), length=9)))

#draw heatmap
fontsize = 12

hm.parameters <- list(fold_change, 
                      breaks = bk2,
                      color = colors2,
                      cellwidth = 20, cellheight = 20, scale = "none",
                      fontsize = fontsize, fontsize_row = fontsize,
                      fontsize_col = fontsize,
                      kmeans_k = NA,
                      show_rownames = T, show_colnames = T,
                      clustering_method = "average",
                      cluster_rows = FALSE, cluster_cols = FALSE,
                      filename = "../Figures/fold_change_difference_heat_map.png")

# To draw the heat map on screen 
do.call("pheatmap", hm.parameters)


