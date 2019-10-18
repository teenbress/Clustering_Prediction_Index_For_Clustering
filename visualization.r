##visualization plots and tables


library(scatterplot3d)
library(ggplot2)

##plot scatter3D;
d <- sim_a()
d <- sim_b()$data
scatterplot3d(brain100, type = "h", highlight.3d = FALSE,color = as.factor(brain.y), angle = 45, scale.y = 0.7, pch = 16, main = "glioblastoma")

##add 3D plant;
lms <- lm(h[,3] ~ h[,1] + h[,2])

sh$plane3d(lms)

data <- read.csv("sim100.csv")

colnames(data) <- c("Index","model","K1","K2","K3","K4","K5","K6","K7","K8","K9","K10")

m1 <- subset(data, model == "1")
m2 <- subset(data, model == "2")
m3 <- subset(data, model == "3")
m4 <- subset(data, model == "4")
m5 <- subset(data, model == "5")
m6 <- subset(data, model == "6")
m7 <- subset(data, model == "7")
m8 <- subset(data, model == "8")

m1[,3:8] <- m1[,3:8]*2
m2[,3:8] <- m2[,3:8]*2
m3[,3:8] <- m3[,3:8]*2
m4[,3:8] <- m4[,3:8]*2
m5[,3:8] <- m5[,3:8]*2
m6[,3:8] <- m6[,3:8]*2
m7[,3:8] <- m7[,3:8]*2
m8[,3:8] <- m8[,3:8]*2

# M4 Group with Joules;
#ggplot(m4,aes(x=Day,y=L.olf,color=as.factor(ID)))+geom_point()+geom_line()+facet_grid(.~as.factor(Joules))# An outlier,quit;

ggplot(data,aes(x = Day, y = R.olf, color = as.factor(ID))) + geom_point() + geom_line() + facet_grid(.~as.factor(Joules)) #1_Jou_line;

ggplot(m4, aes(x = as.factor(Day), y = R.olf)) + geom_point() + geom_boxplot(aes(group = Day, fill = as.factor(Day))) + facet_grid(.~ as.factor(Joules)) + xlab("Day") + labs(fill = "Day")#1_Jou_boxplot;
ggplot(m4, aes(x = as.factor(Day), y = R.olf, color = as.factor(Joules))) + geom_point() + geom_line(aes(group = ID)) + xlab("Day") + labs(fill = 'Joules') #1_Jou_spaghetti;
ggplot(m4, aes(x = as.factor(Day), y = R.olf)) + geom_point() + geom_boxplot(aes(group = Day,fill = as.factor(Day))) + facet_grid(. ~ as.factor(Joules)) + xlab("Day") + labs(fill = "Day") + stat_summary(fun.y = mean, geom = "line", aes(group = 1)) + stat_summary(fun.y = mean, geom = "point")
#M4 Group with Daytime;
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot(p5,p6,p7,p8,cols = 2)

p1 <- ggplot(m1,aes(x=Index,  y=K1,col = Index)) + geom_bar(stat = "identity",aes(fill = Index)) + ylab("Accuraty") + labs(title = "Model1,k = 1")
p2 <- ggplot(m2,aes(x=Index,  y=K3,col = Index)) + geom_bar(stat = "identity",aes(fill = Index)) + ylab("Accuraty") + labs(title = "Model2,k = 3")
p3 <- ggplot(m3,aes(x=Index,  y=K4,col = Index)) + geom_bar(stat = "identity",aes(fill = Index)) + ylab("Accuraty") + labs(title = "Model3,k = 4")
p4 <- ggplot(m4,aes(x=Index,  y=K4,col = Index)) + geom_bar(stat = "identity",aes(fill = Index)) + ylab("Accuraty") + labs(title = "Model4,k = 4")
p5 <- ggplot(m5,aes(x=Index,  y=K2,col = Index)) + geom_bar(stat = "identity",aes(fill = Index)) + ylab("Accuraty") + labs(title = "Model5,k = 2")
p6 <- ggplot(m6,aes(x=Index,  y=K2,col = Index)) + geom_bar(stat = "identity",aes(fill = Index)) + ylab("Accuraty") + labs(title = "Model6,k = 2")
p7 <- ggplot(m7,aes(x=Index,  y=K2,col = Index)) + geom_bar(stat = "identity",aes(fill = Index)) + ylab("Accuraty") + labs(title = "Model7,k = 2")
p8 <- ggplot(m8,aes(x=Index,  y=K3,col = Index)) + geom_bar(stat = "identity",aes(fill = Index)) + ylab("Accuraty") + labs(title = "Model8,k = 3")
