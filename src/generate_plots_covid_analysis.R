# data manipulation & visualization
library(tidyverse)
# path management
library(here)
# functional box plot
library(fda)
# combining plots
library(grid)
library(gridExtra)
library(ggpubr)

source(here("src","fbplot.R"))

select <- dplyr::select

data_dir <- "data/weekly"
results_dir <- "results/covid/covid_rds/"
output_dir <- "results/covid/covid_fitted_values/"
data_names <- list.files(here(data_dir))
results_names <- list.files(here(results_dir))

I_col <- "black"
R_col <- "purple"
I_est_col <- "red"
I_fb_col <- "pink"
R_est_col <- '#006400'
R_fb_col <- 'yellowgreen'
 

fig_dpi <- 300
tick_label_size <-15
label_size <- 20
theme_text_size <- theme(axis.text.x = element_text(size=tick_label_size),
                         axis.text.y = element_text(size=tick_label_size),
                         axis.ticks.length=unit(0.25,'cm'),
                         text = element_text(size=label_size))
scale_factor <- 1
y_lab <- "Number of persons (10,000)"

line_size  <- 2 
fig_label <- c("Infected",
               "Est. Infected",
               "Recovered",
               "Est. Recovered")
line_alpha <- 0.9

generate_fbp_plot <- function(key_word){
  df <- read_csv(here(data_dir,paste0(data_names[str_detect(data_names,key_word)])))
  rds_list <- results_names[str_detect(results_names,key_word)]
  # determine max scale
  max_y <- ceiling(lapply(rds_list,function(dat){
    tmp <- read_rds(here(results_dir,dat))
    quantile(c(tmp[[2]],tmp[[3]]),.99)}) %>% 
    unlist() %>% 
    max(.)) + 1

  list_gg <- lapply(rds_list,function(rds_name){
    subtitle <- if_else(str_detect(rds_name,"smc"),"SMC",
                        if_else(str_detect(rds_name,"mcmc"),"MCMC",
                                "Rej"))
    
    df_result <- read_rds(here(results_dir,rds_name))
    I <- df_result[[2]]
    R <- df_result[[3]]
  
    df_I <- custom_fbplot(I,color='gray70',barcol='gray70',outliercol=0,plot=F)
    df_I <- df_I * scale_factor
    df_I <- df_I %>% 
      mutate(type = "I_est")
    
    df_R <- custom_fbplot(R,color='gray70',barcol='gray70',outliercol=0,plot=F)
    df_R <- df_R * scale_factor
    df_R <- df_R %>% 
      mutate(type = "R_est")
    
    df_gg <- df %>% 
      mutate(I = I*scale_factor,
             R = R*scale_factor) %>% 
      dplyr::select(-S) %>% 
      # cbind(I_est = df_I$y) %>% 
      # cbind(R_est = df_R$y) %>% 
      pivot_longer(cols=-1,names_to="type")
    
    gg <- ggplot() +
      geom_line(data=df_gg,aes(x=t,y=value,color=type),size=line_size,alpha=line_alpha) +
      geom_line(data=df_I,
                aes(x=df$t,y=y,color=type),size=line_size,alpha=line_alpha) +
      geom_ribbon(data=df_I,
                  aes(x=df$t,ymin=lower,ymax=upper,color=NA,fill=type),alpha=0.5)+
      geom_line(data=df_R,
                aes(x=df$t,y=y,color=type),size=line_size,alpha=line_alpha) +
      geom_ribbon(data=df_R,
                  aes(x=df$t,ymin=lower,ymax=upper,color=NA,fill=type),alpha=0.5) +
      ylim(0,max_y)+
      theme_classic() + 
      ylab(y_lab) +
      xlab("") +
      labs(subtitle=subtitle)+
      scale_color_manual(labels=fig_label,
                         values=c("I"=I_col,
                                  "R"=R_col,
                                  "I_est"=I_est_col,
                                  "R_est"=R_est_col)) +
      scale_fill_manual(values=c("I_est"=I_fb_col,
                                 "R_est"=R_fb_col)) +
      scale_x_date(date_labels="%m/%y") +
      guides(fill=FALSE) +
      guides(color=guide_legend(override.aes = list(fill=NA)))+
      theme(legend.title=element_blank(),
            legend.position='bottom',
            plot.subtitle = element_text(hjust=0.5))+
      theme_text_size
    gg
  })
}

# combine the three
# take y axis & legend

combine_plots <- function(gg_list,title){
  lg <- ggpubr::get_legend(gg_list[[1]],position='top')
  
  gg_list <- lapply(gg_list,function(gg){
    gg +
      ylab(element_blank())
      # theme(legend.position='none')
  })
  
  fig <- ggpubr::ggarrange(gg_list[[3]],gg_list[[1]],gg_list[[2]],ncol=3,
                    common.legend=T,
                    legend='top')
  fig <- annotate_figure(fig,
                  left = text_grob(y_lab,rot=90,size = 20),
                  top = text_grob(title,size=25,face = "bold"))
  fig
}

# nice fig size: 15.9 x 6.79 in

# pre vac
gg_pre <- generate_fbp_plot("pre")
ggsave(here("results","covid","covid_fitted_values","pre.png"),plot = combine_plots(gg_pre,"Pre"),dpi = fig_dpi)

# post vac
gg_post <- generate_fbp_plot("post")
ggsave(here("results","covid","covid_fitted_values","post.png"),plot = combine_plots(gg_post,"Post"),dpi = fig_dpi)
