# data manipulation & visualization
library(tidyverse)
# path management
library(here)
# combining plots
library(grid)
library(gridExtra)
library(ggpubr)

results_dir <- here("results","simulation_lorenz")

df_calibration <- read_csv(paste0(results_dir,"/calibration.csv"))

num <- 2^(0:8)

fig_dpi <- 300
tick_label_size <-20
label_size <- 20
theme_text_size <- theme(axis.text.x = element_text(size=tick_label_size),
                         axis.text.y = element_text(size=tick_label_size),
                         axis.ticks.length=unit(0.25,'cm'),
                         text = element_text(size=label_size))

lapply(num,function(x){
  colMeans(df_calibration[1:x,])
}) %>% 
  do.call(rbind,.) %>% 
  cbind(num) %>% 
  as.data.frame() %>% 
  pivot_longer(cols=1:3,values_to="value",names_to="method") -> gg.calibration

ggplot(data=gg.calibration %>% 
         mutate(method= if_else(method=="ABC","Rej",method)) %>% 
         mutate(method = factor(method,levels=c("Rej","MCMC","SMC"))) %>% 
         filter(num>9),aes(x=log(num,2),y=value,col=method)) +
  geom_line(alpha=0.9,size=2) +
  ylim(0,1) +
  geom_hline(yintercept = 0.9,col='purple',size=2,alpha=0.9) +
  theme_classic() +
  xlab("log_2(Number of experiments)") +
  ylab("Calibration value") +
  theme_text_size +
  
  theme(legend.position = 'top',
        legend.title = element_blank()) -> gg.cal

ggsave(paste0(results_dir,"/calibration.png"),plot = gg.cal,dpi = fig_dpi)

df_ess <- read_csv(paste0(results_dir,"/ess.csv"))
# filter(ABC<=5000)
ess <- round(apply(df_ess,2,summary),0)
df_cpu <- read_csv(paste0(results_dir,"/cpu.csv"))[-which.max(df_ess$ABC),]
cpu <- round(apply(df_cpu,2,summary),2)
df_ess_cpu <- read_csv(paste0(results_dir,"/ess_cpu.csv"))[-which.max(df_ess$ABC),]
ess_cpu <- round(apply(df_ess_cpu,2,summary),0)

ess %>% 
  cbind(cpu) %>% 
  cbind(ess_cpu) -> tab
openxlsx::write.xlsx(tab %>% as.data.frame(),here(paste0(results_dir),"/tab.xlsx"),row.names=T)

