
library(bayesmove)
library(ggplot2)

data(tracks.list)

#subset only first track
tracks.list<- dplyr::bind_rows(tracks.list)

#only retain id and discretized step length (SL) and turning angle (TA) columns
tracks<- subset(tracks.list, select = c(id, SL, TA))


set.seed(1)

# Define model params
alpha=0.1
ngibbs=1000
nburn=ngibbs/2
nmaxclust=7

model1=cluster_obs(dat=tracks, alpha=alpha, ngibbs=ngibbs, nmaxclust=nmaxclust, nburn=nburn)


plot(model1$loglikel, type = "l")
plot(model1$gamma1, type = "l")

MAP.iter<- bayesmove:::get_MAP_internal(dat = model1$loglikel, nburn = nburn)

theta<- model1$theta[MAP.iter,]
names(theta)<- 1:length(theta)
theta<- sort(theta, decreasing = TRUE)
theta %>% cumsum()  #first 5 states optimal (i.e. > 90%)

ord<- as.numeric(names(theta))


behav.res<- get_behav_hist(dat = model1, nburn = nburn, ngibbs = ngibbs, nmaxclust = nmaxclust,
                            var.names = c("Step Length","Turning Angle"), ord = ord,
                            MAP.iter = MAP.iter)



# Plot histograms of proportion data

ggplot(behav.res, aes(x = bin, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = c(viridis::viridis(5),
                               "grey35","grey35"), guide = FALSE) +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  scale_x_continuous(breaks = 1:8) +
  facet_grid(behav ~ var, scales = "free_x")







