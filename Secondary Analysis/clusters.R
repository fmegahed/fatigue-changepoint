# This file is written to generate the profile of walking per cluster
# by averaging the values for each feature across the participants within
# the same cluster.
# The file is sourced by the RMD, we did not include this code directly
# in the RMD since it is long and repetitive 


# [1] Creating the Clusters
# C1 Membership based on the Dendogram
stride.len.C1 <- cbind(subject1_medianf_norm[,2],
                       subject11_medianf_norm[,2],
                       subject15_medianf_norm[,2],
                       subject10_medianf_norm[,2]) %>% rowMeans()
stride.height.C1 <- cbind(subject1_medianf_norm[,3],
                       subject11_medianf_norm[,3],
                       subject15_medianf_norm[,3],
                       subject10_medianf_norm[,3]) %>% rowMeans()
stride.dur.C1 <- cbind(subject1_medianf_norm[,4],
                       subject11_medianf_norm[,4],
                       subject15_medianf_norm[,4],
                       subject10_medianf_norm[,4]) %>% rowMeans()
C1 <- cbind(subject1_medianf_norm[,1],
            stride.len.C1, stride.height.C1, stride.dur.C1)
colnames(C1)[1] <- "percent.time.from.start"


# C2 based on the dendogram
stride.len.C2 <- cbind(subject2_medianf_norm[,2],
                       subject8_medianf_norm[,2],
                       subject12_medianf_norm[,2]) %>% rowMeans()
stride.height.C2 <- cbind(subject2_medianf_norm[,3],
                          subject8_medianf_norm[,3],
                          subject12_medianf_norm[,3]) %>% rowMeans()
stride.dur.C2 <- cbind(subject2_medianf_norm[,4],
                       subject8_medianf_norm[,4],
                       subject12_medianf_norm[,4]) %>% rowMeans()
C2 <- cbind(subject1_medianf_norm[,1],
            stride.len.C2, stride.height.C2, stride.dur.C2)
colnames(C2)[1] <- "percent.time.from.start"


# C3 membership based on the dendogram
stride.len.C3 <- cbind(subject4_medianf_norm[,2],
                       subject14_medianf_norm[,2]) %>% rowMeans()
stride.height.C3 <- cbind(subject4_medianf_norm[,3],
                          subject14_medianf_norm[,3]) %>% rowMeans()
stride.dur.C3 <- cbind(subject4_medianf_norm[,4],
                       subject15_medianf_norm[,4]) %>% rowMeans()
C3 <- cbind(subject1_medianf_norm[,1],
            stride.len.C3, stride.height.C3, stride.dur.C3)
colnames(C3)[1] <- "percent.time.from.start"

# C4
stride.len.C4 <- cbind(subject5_medianf_norm[,2],
                       subject6_medianf_norm[,2],
                       subject7_medianf_norm[,2],
                       subject9_medianf_norm[,2],
                       subject13_medianf_norm[,2]) %>% rowMeans()
stride.height.C4 <- cbind(subject5_medianf_norm[,3],
                          subject6_medianf_norm[,3],
                          subject7_medianf_norm[,3],
                          subject9_medianf_norm[,3],
                          subject13_medianf_norm[,3]) %>% rowMeans()
stride.dur.C4 <- cbind(subject5_medianf_norm[,4],
                       subject6_medianf_norm[,4],
                       subject7_medianf_norm[,4],
                       subject9_medianf_norm[,4],
                       subject13_medianf_norm[,4]) %>% rowMeans()
C4 <- cbind(subject1_medianf_norm[,1],
            stride.len.C4, stride.height.C4, stride.dur.C4)
colnames(C4)[1] <- "percent.time.from.start"


# [2] Phase II: Plotting the data using ggplot2
# ---------------------------------------------------
df <- cbind(C1[,1:2],C2[,2],C3[,2], C4[,2]) %>% data.frame()
colnames(df) <- c("percent.time.from.start",paste0("stride.len.C",seq(1,4,1)))
meltdf.stride.length <- reshape2::melt(df, id.vars="percent.time.from.start",
                                       measure.vars=c(paste0("stride.len.C",seq(1,4,1)))) 
g1 <- ggplot(meltdf.stride.length,
       aes(x=percent.time.from.start,y=value, colour=variable)) + 
  geom_line() +
  theme_bw(base_size = 14) +
  scale_colour_brewer(palette = "Dark2") +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Average Scaled Stride Lengths per Cluster")


df <- cbind(C1[,c(1,3)],C2[,3],C3[,3], C4[,3]) %>% data.frame()
colnames(df) <- c("percent.time.from.start",paste0("stride.height.C",seq(1,4,1)))
meltdf.stride.height <- reshape2::melt(df, id.vars="percent.time.from.start",
                                       measure.vars=c(paste0("stride.height.C",seq(1,4,1)))) 

g2 <- ggplot(meltdf.stride.height,
             aes(x=percent.time.from.start,y=value, colour=variable)) + 
  geom_line() +
  theme_bw(base_size = 14) +
  scale_colour_brewer(palette = "Dark2") +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Average Scaled Stride Height per Cluster")



df <- cbind(C1[,c(1,4)],C2[,4],C3[,4], C4[,4]) %>% data.frame()
colnames(df) <- c("percent.time.from.start",paste0("stride.dur.C",seq(1,4,1)))
meltdf.stride.dur <- reshape2::melt(df, id.vars="percent.time.from.start",
                                       measure.vars=c(paste0("stride.dur.C",seq(1,4,1)))) 

g3 <- ggplot(meltdf.stride.dur,
             aes(x=percent.time.from.start,y=value, colour=variable)) + 
  geom_line() +
  theme_bw(base_size = 14) +
  scale_colour_brewer(palette = "Dark2") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5) ) +
  ggtitle("Average Stride Duration per Cluster")


# [3] Printing the plot with a common legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(g3)
gridExtra::grid.arrange((g1), 
                        (g2), 
                        (g3+ theme(legend.position="none")),
                        mylegend, nrow = 4,
                        heights=c(5,5,5,1))

# [4] Saving the Plot with a Resolution of choice
gridPlot <- arrangeGrob((g1), 
                        (g2), 
                        (g3+ theme(legend.position="none",
                                   legend.title=element_text(size=14), 
                                   legend.text=element_text(size=14))),
                        mylegend, nrow = 4,
                        heights=c(5,5,5,1))
ggsave('dtwClusters.png', gridPlot, height = 11, width= 8.5,
       units = "in", dpi=600)
