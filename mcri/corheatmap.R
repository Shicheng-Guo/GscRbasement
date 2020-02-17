

library("corrplot")
nba <- as.matrix(read.csv("https://raw.githubusercontent.com/Shicheng-Guo/Shicheng-Guo.Github.io/master/data/ppg2008.csv")[-1])
res1 <- cor.mtest(nba, conf.level = .95)
par(mfrow=c(2,2))

# correlation and P-value
corrplot(cor(nba), p.mat = res1$p, insig = "label_sig",sig.level = c(.001, .01, .05), pch.cex = 0.8, pch.col = "white",tl.cex=0.8)

# correlation and hclust
corrplot(cor(nba), method = "shade", outline = T, addgrid.col = "darkgray", order="hclust", 
         mar = c(4,0,4,0), addrect = 4, rect.col = "black", rect.lwd = 5,cl.pos = "b", tl.col = "indianred4", 
         tl.cex = 0.8, cl.cex = 0.8)

