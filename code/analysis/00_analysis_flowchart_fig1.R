source(here::here("code/packages.R"))
source(here::here("code/file_paths.R"))

XX <- c("psoriasis", "eczema")

for(exposure in XX){
	#exposure <- XX[1]
	ABBRVexp <- str_sub(exposure, 1 , 3)
	
	# total with diagnosis code -----------------------------------------------
	patid_file <- haven::read_dta(paste0(datapath, "in/Patient_extract_", ABBRVexp, "_extract_1.dta"))
	total_with_code <- unique(patid_file$patid) %>% length()
	
	# total that meet algorithm -----------------------------------------------
	exposed_file <- haven::read_dta(paste0(datapath, "out/", exposure, "Exposed.dta"))
	total_meet_algorithm <- unique(exposed_file$patid) %>% length()
	
	# with eligible follow up  ------------------------------------------------
	fup_file <- haven::read_dta(paste0(datapath, "out/", exposure, "Exposed-eligible-mhealth.dta"))
	total_with_fup <- unique(fup_file$patid) %>% length()
	print("Got the full cohort")
	print(total_with_fup)
	
	# eligible match ----------------------------------------------------------
	matched_cohort <- haven::read_dta(paste0(datapath, "out/getmatchedcohort-", exposure, "-main-mhealth.dta"))
	
	matched_cohort$tempid <- paste0(matched_cohort$setid, matched_cohort$patid)
	total_matched_exp <- unique(matched_cohort$tempid[matched_cohort$exposed == 1]) %>% length()
	total_matched_unexp <- unique(matched_cohort$tempid[matched_cohort$exposed == 0]) %>% length()
	total_matched_cohort <- total_matched_exp + total_matched_unexp
	
	print("Got the full Matched cohort")
	print(total_matched_cohort)
	
	# excluded based on outcome  ----------------------------------------------
	var_depressionAll <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-depression-all.dta"))
	var_anxietyAll <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-anxiety-all.dta"))
	var_depressionCensor <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-depression-censor.dta"))
	var_anxietyCensor <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-anxiety-censor.dta"))
	
	## censor depression/anxiety based on censor codes (SMI)
	var_dep <- var_depressionAll %>%
		rename(depdate = eventdate) %>%
		left_join(var_depressionCensor, by = "patid") %>%
		rename(censorDepDate = eventdate) %>%
		mutate_at("censorDepDate", ~ifelse(is.na(.), as.Date("3000-01-01"), as.Date(.))) %>%
		mutate_at("censorDepDate", ~as.Date(., origin = "1970-01-01")) %>%
		mutate(dep = case_when(
			depdate < censorDepDate &
				!is.na(depdate) ~ 1,
			TRUE ~ 0
		)) %>%
		select(patid, case = dep, eventdate = depdate)
	
	var_anx <- var_anxietyAll %>%
		rename(anxdate = eventdate) %>%
		left_join(var_anxietyCensor, by = "patid") %>%
		rename(censorAnxDate = eventdate) %>%
		mutate_at("censorAnxDate", ~ifelse(is.na(.), as.Date("3000-01-01"), as.Date(.))) %>%
		mutate_at("censorAnxDate", ~as.Date(., origin = "1970-01-01")) %>%
		mutate(anx = case_when(
			anxdate < censorAnxDate &
				!is.na(anxdate) ~ 1,
			TRUE ~ 0
		)) %>%
		select(patid, case = anx, eventdate = anxdate)
	
	## ANXIETY
	select_vars <- c("enddate", "eventdate")
	
	an_anxiety <- matched_cohort %>% 
		ungroup() %>%
		left_join(var_anx, by = "patid") %>%
		mutate(newenddate = apply(select(., all_of(select_vars)), 1, FUN = min, na.rm = T)) %>%	
		filter(newenddate>indexdate) %>%
		select(-newenddate) ## have used this data enough now
	
	## check valid sets
	invalid_sets <- an_anxiety %>%
		group_by(setid) %>%
		summarise(valid=max(exposed)) %>%
		filter(valid==0) %>%
		pull(setid)
	an_anxiety_valid <- an_anxiety %>%
		filter(!setid %in% invalid_sets)
	
	an_anxiety$tempid <- paste0(an_anxiety$setid, an_anxiety$patid)
	total_anxiety_cohort_all <- unique(an_anxiety$tempid) %>% length()
	lost_to_censoring_anxiety <- total_matched_cohort - total_anxiety_cohort_all
	
	an_anxiety_valid$tempid <- paste0(an_anxiety_valid$setid, an_anxiety_valid$patid)
	total_anxiety_cohort <- unique(an_anxiety_valid$tempid) %>% length()
	lost_to_invalid_sets_anxiety <- total_anxiety_cohort_all - total_anxiety_cohort
	
	total_anxiety_cohort_exp <- unique(an_anxiety_valid$tempid[an_anxiety_valid$exposed==1]) %>% length()
	total_anxiety_cohort_unexp <- unique(an_anxiety_valid$tempid[an_anxiety_valid$exposed==0]) %>% length()
	
	## DEPRESSION
	an_depression <- matched_cohort %>% 
		ungroup() %>%
		left_join(var_dep, by = "patid") %>%
		mutate(newenddate = apply(select(., all_of(select_vars)), 1, FUN = min, na.rm = T)) %>%
		filter(newenddate>indexdate) %>%
		select(-newenddate) ## have used this data enough now
	
	## check valid sets
	invalid_sets <- an_depression %>%
		group_by(setid) %>%
		summarise(valid=max(exposed)) %>%
		filter(valid==0) %>%
		pull(setid)
	an_depression_valid <- an_depression %>%
		filter(!setid %in% invalid_sets)
	
	an_depression$tempid <- paste0(an_depression$setid, an_depression$patid)
	total_depression_cohort_all <- unique(an_depression$tempid) %>% length()
	lost_to_censoring_depression <- total_matched_cohort - total_depression_cohort_all
	
	an_depression_valid$tempid <- paste0(an_depression_valid$setid, an_depression_valid$patid)
	total_depression_cohort <- unique(an_depression_valid$tempid) %>% length()
	lost_to_invalid_sets_depression <- total_depression_cohort_all - total_depression_cohort
	
	total_depression_cohort_exp <- unique(an_depression_valid$tempid[an_depression_valid$exposed==1]) %>% length()
	total_depression_cohort_unexp <- unique(an_depression_valid$tempid[an_depression_valid$exposed==0]) %>% length()
	
	print("Got all the values. Making the graph...")
	# make the diagram --------------------------------------------------------
	paste_flow <- function(text, number = NA){
		paste0(text, ": ", prettyNum(number, big.mark = ","))
	}
	
	text1 <- paste_flow(paste0("Individuals in CPRD with a\n diagnosis code for ", exposure), total_with_code)
	text2 <- paste_flow("Excluded as do not\n meet exposure algorithm*", total_with_code - total_meet_algorithm)
	text3 <- paste_flow("Individuals remaining", total_meet_algorithm)
	text4 <- paste_flow("Excluded as do not contribute\n follow-up time during\n the study period", total_meet_algorithm - total_with_fup)
	text5 <- paste_flow("Exposed and eligible to be matched", total_with_fup)
	text6 <- paste_flow(paste0("Eligible individuals with ", exposure, "\n and at least 1 eligible matched\n individual without ", exposure), total_matched_exp)
	text7 <- paste_flow(paste0("Eligible individuals without ", exposure, "\n and at least 1 eligible matched\n individual with ", exposure), total_matched_unexp)
	text8 <- paste_flow("Anxiety cohort", total_matched_cohort)
	text9 <- paste_flow("Excluded as individual had pre-existing\n anxiety diagnosis", lost_to_censoring_anxiety)
	text10 <- paste_flow("Eligible individuals", total_anxiety_cohort_all)
	text11 <- paste_flow("Excluded as no eligible exposed individual\n in the matched set", lost_to_invalid_sets_anxiety)
	text12 <- paste0(paste_flow("Final anxiety cohort", total_anxiety_cohort),
									 "\n",
									 paste_flow("Exposed", total_anxiety_cohort_exp),
									 ", ", 
									 paste_flow("Unxposed", total_anxiety_cohort_unexp))
	text13 <- paste_flow("Depression cohort", total_matched_cohort)
	text14 <- paste_flow("Excluded as individual had pre-existing\n depression diagnosis", lost_to_censoring_depression)
	text15 <- paste_flow("Eligible individuals", total_depression_cohort_all)
	text16 <- paste_flow("Excluded as no eligible exposed individual\n in the matched set", lost_to_invalid_sets_depression)
	text17 <- paste0(paste_flow("Final depression cohort", total_depression_cohort),
									 "\n",
									 paste_flow("Exposed", total_depression_cohort_exp),
									 ", ", 
									 paste_flow("Unxposed", total_depression_cohort_unexp))
	text18 <- paste_flow("Excluded as no eligible match for\n the exposed individual", total_with_fup - total_matched_exp)
	
	graph <- "digraph {
				graph [layout = dot, rankdir = TB]
	
	      # node definitions with substituted label text
	      node [fontname = Helvetica, shape = rectangle]
	      nodesep=0.5; #hack
	      
	      tab1 [label = '@@1']
	      tab2 [label = '@@2']
	      tab3 [label = '@@3']
	      tab4 [label = '@@4']
	      tab5 [label = '@@5']
	      tab6 [label = '@@6']
	      tab7 [label = '@@7']
	      tab8 [label = '@@8']
	      tab9 [label = '@@9']
	      tab10 [label = '@@10']
	      tab11 [label = '@@11']
	      tab12 [label = '@@12']
	      tab13 [label = '@@13']
	      tab14 [label = '@@14']
	      tab15 [label = '@@15']
	      tab16 [label = '@@16']
	      tab17 [label = '@@17']
	      tab18 [label = '@@18']
	
	      # edge definitions with the node IDs
	    	edge[minlen = 2, color = grey, penwidth = 2]
	      
	      tab1 -> tab3 -> tab5 -> tab6 -> tab8;
				tab1 -> tab2;
	      tab3 -> tab4;
	      tab5 -> tab18;
	      tab8 -> tab9;
	      tab10 -> tab11;
	      tab13 -> tab14;
	      tab15 -> tab16;
	      tab6 -> tab13;
	      tab8 -> tab10 -> tab12;
	      tab13 -> tab15 -> tab17;
	      tab7 -> tab8;
	      tab7 -> tab13;
	      
	      subgraph main {
				  {rank = same; tab1; tab2}
					{rank = same; tab3; tab4}
					{rank = same; tab5; tab18}
					{rank = same; tab6; tab7}
			  }
	      subgraph cluster_1 {
	        label = 'Anxiety cohort';
	        color = Blue;
	        fontname = Helvetica
			    {rank = same; tab8; tab9}
			    {rank = same; tab10; tab11}
			    tab12;
			  }
	      subgraph cluster_2 {
	        label = 'Depression cohort'
	        color = Red
	        fontname = Helvetica
			    {rank = same; tab13; tab14}
			    {rank = same; tab15; tab16}
			    tab17
	      }
	
	      }
	
	      [1]: text1
	      [2]: text2
	      [3]: text3
	      [4]: text4
	      [5]: text5
	      [6]: text6
	      [7]: text7
	      [8]: text8
	      [9]: text9
	      [10]: text10
	      [11]: text11
	      [12]: text12
	      [13]: text13
	      [14]: text14
	      [15]: text15
	      [16]: text16
	      [17]: text17
	      [18]: text18
	"
	grViz(graph)
	grViz(graph) %>% 
		export_svg %>% 
		charToRaw %>% 
		rsvg_pdf(here("out", "diagram", paste0("flowchart_", exposure, ".pdf")))

}
