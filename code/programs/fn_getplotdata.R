make_df_forest <- function(){
  if(!exists("ecz_table")){stop("Need to make the results tables first (ecz_table and pso_table)")}
  get_plot_data <- function(pretty_table) {
    numeric_tab <- pretty_table %>%
      dplyr::select(-starts_with("n")) %>%
      dplyr::select(-starts_with("events")) %>%
      dplyr::select(-characteristic,-outcome) 
    for(ii in 1:4){
      if(paste0("ci", ii) %in% colnames(numeric_tab)){
        numeric_tab <- numeric_tab %>% 
          separate(eval(parse(text = (paste0("ci", ii)))), paste0(c("ciL", "ciU"),ii), sep = "-") 
      }
    }
    numeric_tab <- numeric_tab %>% 
      drop_na() %>% 
      filter(hr2 != "-") %>% 
      mutate_if(is.character, as.numeric) 
    
    get_n <- pretty_table %>% 
      dplyr::select(starts_with("n")) %>% 
      filter(n2 != " ") %>% 
      mutate_all(~str_remove_all(.,",")) %>% 
      mutate_if(is.character, as.numeric) %>% 
      mutate(outcome = c(rep("Anxiety",2), rep("Depression",2))) %>% 
      group_by(outcome) %>% 
      summarise_all(~sum(.)) %>% 
      dplyr::select(-outcome)
    
    plot_tab1 <- pretty_table %>%
      dplyr::select(2) %>%
      filter(outcome != " ") %>%
      bind_cols(numeric_tab, get_n) 
  }
  
  # uses the little function loaded at the top of the script
  ecz_plot <- get_plot_data(ecz_table) %>%
    mutate(exposure = "Eczema")
  pso_plot <- get_plot_data(pso_table) %>%
    mutate(exposure = "Psoriasis")
  
  plot_df <- bind_rows(ecz_plot, pso_plot)
  plot_df <- plot_df %>% 
    pivot_longer(cols = starts_with(c("hr", "ci", "n")),
                 names_to = c("metric", "model"),
                 names_pattern = "(.*)([0-9])") %>% 
    pivot_wider(id_cols = c(outcome, exposure, model), names_from = metric)
  
  ##  models as factors 
  plot_df <- plot_df %>%
    mutate(model = factor(model, levels = 1:3, labels = c("Crude", "Confounder Adjusted", "Mediator adjusted")))
  
  ## add alpha parameter to grey out other models
  plot_df$a <- 0.5
  plot_df$a[plot_df$model == "Mediator adjusted"] <- 1
  
  # convert results to string and pad so they have same width for printing on plots
  plot_df$text_hr <- str_pad(round(plot_df$hr,2), 4, pad = "0", side = "right")
  plot_df$text_ciL <- str_pad(round(plot_df$ciL,2), 4, pad = "0", side = "right")
  plot_df$text_ciU <- str_pad(round(plot_df$ciU,2), 4, pad = "0", side = "right")
  plot_df$text_n <- prettyNum(plot_df$n, big.mark = ",")
  plot_df$text_to_plot <- str_pad(paste0("(",
                                 plot_df$text_n,
                                 ") ",
                                 plot_df$text_hr, 
                                 " [", 
                                 plot_df$text_ciL,
                                 ",", 
                                 plot_df$text_ciU,
                                 "]"),
                                 28, pad = " ", side = "left")
  
  plot_df
}
