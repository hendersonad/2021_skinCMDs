source(here::here("code/packages.R"))
source(here::here("code/file_paths.R"))

XX <- c("psoriasis", "eczema")

for (exposure in XX) {
  
  ABBRVexp <- str_sub(exposure, 1 , 3)
  
  # load data ---------------------------------------------------------------
  df_model_anx <-
    readRDS(paste0(
      datapath,
      "out/df_model",
      ABBRVexp,
      "_anxiety.rds"
    ))
  df_model_dep <-
    readRDS(paste0(
      datapath,
      "out/df_model",
      ABBRVexp,
      "_depression.rds"
    ))
  
  get_crude_RD <- function(data) {
    df_1row <- data %>%
      #slice(1:1000) %>%
      group_by(setid, patid, exposed) %>%
      summarise(pyears = sum(t / 365.25), out = sum(out))
    
    df_crude <- df_1row %>%
      ungroup() %>%
      group_by(exposed) %>%
      summarise(pyears = sum(pyears), events = sum(out)) %>%
      mutate(crude_rate = (events / pyears) * 1000)
    
    events_con <-
      df_crude %>% filter(exposed == "Unexposed") %>% select(events) %>% pull()
    events_exp <-
      df_crude %>% filter(exposed != "Unexposed") %>% select(events) %>% pull()
    pyears_con <-
      df_crude %>% filter(exposed == "Unexposed") %>% select(pyears) %>% pull()
    pyears_exp <-
      df_crude %>% filter(exposed != "Unexposed") %>% select(pyears) %>% pull()
    
    RD <-
      fmsb::ratedifference(
        events_exp,
        events_con,
        pyears_exp,
        pyears_con,
        CRC = T,
        conf.level = 0.99
      )
    RD
    RD_out <- c(NA, RD$estimate)
    RD_lci <- c(NA, RD$conf.int[1])
    RD_uci <- c(NA, RD$conf.int[2])
    out_df <- df_crude %>%
      bind_cols(rd = RD_out, lci = RD_lci, uci = RD_uci)
  }
  anx_rd <-
    get_crude_RD(data = df_model_anx) %>% mutate(out = "Anxiety")
  dep_rd <-
    get_crude_RD(data = df_model_dep) %>% mutate(out = "Depression")
  
  rd_output <- anx_rd %>%
    bind_rows(dep_rd) %>%
    select(out, exposed, everything()) %>%
    mutate_at(c("rd", "lci", "uci"), ~round(. * 1000,2)) %>% 
    mutate_at(c("pyears","crude_rate"), ~round(.,2))
  
  assign(paste0("rd_",ABBRVexp), rd_output)
}

rd_full <- rd_ecz %>% 
  bind_rows(rd_pso)

rd_tab <- rd_full %>% 
  mutate(rd2 = paste0(round(rd, 1), " (", round(lci, 1), "-", round(uci,1), ")")) %>% 
  select(-rd, -lci, -uci) 
rd_tab$rd2[grepl(pattern = "NA", rd_tab$rd2)] <- " "
rd_tab %>% 
  gt() %>%
  tab_row_group(
    label = md("**Atopic eczema**"), 
    rows = 1:4
  ) %>% 
  tab_row_group(
    label = md("**Psoriasis**"), 
    rows = 5:8
  ) %>% 
  row_group_order(c("**Atopic eczema**", "**Psoriasis**")) %>% 
  cols_align(columns = 1:2, align = "left") %>% 
  cols_align(columns = 3:dim(rd_tab)[2], align = "right") %>% 
  fmt_number(
    columns = 3:4,
    decimals = 0,
    use_seps = T
  ) %>% 
  fmt_number(
    columns = 5,
    decimals = 1,
    use_seps = T
  ) %>% 
  cols_label(
    out = md("**Outcome**"),
    exposed = md("**Exposure**"),
    pyears = md("**Person-years**"),
    events = md("**Events**"),
    crude_rate = md("**Crude Rate**"),
    rd2 = md("**Rate difference**")
  ) %>% 
  tab_footnote(
    footnote = "Per 1,000 person-years",
    locations = cells_column_labels(
      columns = c(crude_rate, 
                  rd2))) %>% 
  gt::gtsave(
    filename =  paste0("tab5_cruderates.html"),
    path = here::here("out/tables")
  )
