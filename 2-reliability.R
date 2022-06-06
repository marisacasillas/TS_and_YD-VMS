library(tidyverse)
library(phonfieldwork)
vowels <- "[iɪeɛəoɔaɑu ]"


metadata <- read_delim("ACLEW_list_of_corpora-CasillasOnly-UTF8.csv",
                       delim = ";",
                       locale = locale(encoding = "UTF-8"))

reliability.clips <- read_csv("reliability-clips.csv")

files <- list.files("eafs/", "\\d{4}.eaf$", recursive = TRUE)
all.vocs.tbl <- tibble()
all.agmt.data <- tibble()

# IPA substitutions
# tʃ' goes to $
# tʃ goes to S
# dʒ goes to Z
# p' goes to P
# t' goes to T
# k' goes to K
# multi-character YD symbols aren't otherwise used afaik

for (file in files) {
  print(file)
  
  # read in a file and note which language it comes from 
  filename <- file
  language <- ifelse(grepl("tseltal", filename), "Tseltal", "Yélî Dnye")
  data <- eaf_to_df(paste0("eafs/", filename))
  names(data) <- c("tier_num", "id", "value", "tier", "tier_type",
                   "start", "stop", "file_source")
  data <- mutate(data, length = stop - start) %>%
    select(c("tier", "start", "stop", "length", "value"))
    
  # make sure we're only including vocs from the annotated random clips
  clips <- filter(data, tier == "code" & grepl("random", value))
  annot.clips <- filter(reliability.clips,
                  file == str_extract(filename, "\\d{4}"))$clip
  
  random.clip.vocs <- tibble()
  for (i in 1:nrow(clips)) {
    if (clips$value[i] %in% annot.clips) {
      vocs <- which((data$tier == "pho@CHI" | data$tier == "ph2@CHI") &
                      ((data$start >= clips$start[i] & data$stop <= clips$stop[i])|
                         (data$start < clips$start[i] & data$stop <= clips$stop[i] &
                            data$stop > clips$start[i])|
                         (data$start >= clips$start[i] & data$start <= clips$stop[i] &
                            data$stop > clips$stop[i])))
      vocs.tbl <- bind_cols(data[vocs,],
                            tibble(clip = clips$value[i], language = language, filename = filename))
      random.clip.vocs <- bind_rows(random.clip.vocs, vocs.tbl)
    }
  }
  # substitute 2+char symbols with a single character
  random.clip.vocs$value <- gsub("tʃ'", "$", random.clip.vocs$value)
  random.clip.vocs$value <- gsub("tʃ", "S", random.clip.vocs$value)
  random.clip.vocs$value <- gsub("dʒ", "Z", random.clip.vocs$value)
  random.clip.vocs$value <- gsub("p'", "P", random.clip.vocs$value)
  random.clip.vocs$value <- gsub("t'", "T", random.clip.vocs$value)
  random.clip.vocs$value <- gsub("k'", "K", random.clip.vocs$value)
  # bind individual rec data to table of all rec data
  all.vocs.tbl <- bind_rows(all.vocs.tbl, random.clip.vocs)
  
  bp.annots <- filter(all.vocs.tbl, tier == "pho@CHI")
  names(bp.annots)[which(names(bp.annots) == "value")] <- "BP"
  bp.annots$tier <- NULL

  mc.annots <- filter(all.vocs.tbl, tier == "ph2@CHI")
  names(mc.annots)[which(names(mc.annots) == "value")] <- "MC"
  mc.annots$tier <- NULL
  
  all.vocs.wider <- left_join(mc.annots, bp.annots)
  
  # agreement calculations
  all.vocs.agmt <- all.vocs.wider %>%
    mutate(
      ## canonical babbles vs. other vocs
      ## (note that MC coded all non-canonical babble as N, unlike BP
      ## who annotated L and Y as well--see paper)
      ncb.agr.levels = case_when(
        grepl("[YLN]", BP) & MC == "N" ~ "N",
        !grepl("[YLN]", BP) & MC != "N" ~ "C",
        grepl("[YLN]", BP) & MC != "N" ~ "BN",
        !grepl("[YLN]", BP) & MC == "N" ~ "MN",
        TRUE ~ "oops"
      ),
      ncb.agr.bin = case_when(
        ncb.agr.levels == "N" | ncb.agr.levels == "C" ~ 1,
        TRUE ~ 0
      ),
      pb.delta = str_count(BP, "[pb]") - str_count(MC, "[pb]"),
      pb.detected = str_count(BP, "[pb]") + str_count(MC, "[pb]") > 0,
      td.delta = str_count(BP, "[td]") - str_count(MC, "[td]"),
      td.detected = str_count(BP, "[td]") + str_count(MC, "[td]") > 0,
      kg.delta = str_count(BP, "[kg]") - str_count(MC, "[kg]"),
      kg.detected = str_count(BP, "[kg]") + str_count(MC, "[kg]") > 0,
      m.delta = str_count(BP, "[m]") - str_count(MC, "[m]"),
      m.detected = str_count(BP, "[m]") + str_count(MC, "[m]") > 0,
      n.delta = str_count(BP, "[n]") - str_count(MC, "[n]"),
      n.detected = str_count(BP, "[n]") + str_count(MC, "[n]") > 0,
      l.delta = str_count(BP, "[l]") - str_count(MC, "[l]"),
      l.detected = str_count(BP, "[l]") + str_count(MC, "[l]") > 0
    )
  all.agmt.data <- bind_rows(all.agmt.data, all.vocs.agmt) 
}
all.agmt.data$aclew_id <- as.numeric(str_extract(all.agmt.data$filename, "\\d{4}"))
all.agmt.data <- left_join(all.agmt.data, metadata)

write_csv(all.agmt.data, "all.reliability.data.csv")



# NCB agreement
ncb.by.site.age <- all.agmt.data %>%
  mutate(
    MC.N = case_when(
      ncb.agr.levels == "MN" ~ 1,
      TRUE ~ 0),
    BP.N = case_when(
      ncb.agr.levels == "BN" ~ 1,
      TRUE ~ 0)
  ) %>%
  group_by(aclew_id, age_mo_round, language) %>%
  summarize(
    mean.agmt = mean(ncb.agr.bin),
    n.mcn = sum(MC.N),
    n.bpn = sum(BP.N),
    n.vocs = n()
  ) %>%
  mutate(
    ncb.bias = n.bpn - n.mcn,
    ncb.bias.rel = ncb.bias/n.vocs
  )
round(mean(ncb.by.site.age$mean.agmt),2)
round(mean(subset(ncb.by.site.age, language == "Tseltal")$mean.agmt),2)
round(mean(subset(ncb.by.site.age, language == "Yélî Dnye")$mean.agmt),2)
ggplot(ncb.by.site.age, aes(mean.agmt)) +
  geom_histogram(stat = "bin") +
  facet_grid(~language) +
  labs(
    title = "canonical babble vs. other",
    x = "mean agreement",
    y = "# recordings")

# pb agreement
pb.by.site.age <- all.agmt.data %>%
  filter(pb.detected == TRUE) %>%
  group_by(aclew_id, age_mo_round, language) %>%
  summarize(
    mean.delta = mean(pb.delta),
    median.delta = median(pb.delta),
    sd.delta = sd(pb.delta),
    min.delta = min(pb.delta),
    max.delta = max(pb.delta),
    n.vocs = n()
  )

# td agreement
td.by.site.age <- all.agmt.data %>%
  filter(td.detected == TRUE) %>%
  group_by(aclew_id, age_mo_round, language) %>%
  summarize(
    mean.delta = mean(td.delta),
    median.delta = median(td.delta),
    sd.delta = sd(td.delta),
    min.delta = min(td.delta),
    max.delta = max(td.delta),
    n.vocs = n()
  )

# kg agreement
kg.by.site.age <- all.agmt.data %>%
  filter(kg.detected == TRUE) %>%
  group_by(aclew_id, age_mo_round, language) %>%
  summarize(
    mean.delta = mean(kg.delta),
    median.delta = median(kg.delta),
    sd.delta = sd(kg.delta),
    min.delta = min(kg.delta),
    max.delta = max(kg.delta),
    n.vocs = n()
  )

# m agreement
m.by.site.age <- all.agmt.data %>%
  filter(m.detected == TRUE) %>%
  group_by(aclew_id, age_mo_round, language) %>%
  summarize(
    mean.delta = mean(m.delta),
    median.delta = median(m.delta),
    sd.delta = sd(m.delta),
    min.delta = min(m.delta),
    max.delta = max(m.delta),
    n.vocs = n()
  )

# n agreement
n.by.site.age <- all.agmt.data %>%
  filter(n.detected == TRUE) %>%
  group_by(aclew_id, age_mo_round, language) %>%
  summarize(
    mean.delta = mean(n.delta),
    median.delta = median(n.delta),
    sd.delta = sd(n.delta),
    min.delta = min(n.delta),
    max.delta = max(n.delta),
    n.vocs = n()
  )

# l agreement
l.by.site.age <- all.agmt.data %>%
  filter(l.detected == TRUE) %>%
  group_by(aclew_id, age_mo_round, language) %>%
  summarize(
    mean.delta = mean(l.delta),
    median.delta = median(l.delta),
    sd.delta = sd(l.delta),
    min.delta = min(l.delta),
    max.delta = max(l.delta),
    n.vocs = n()
  )

# deltas figure
deltas.pb <- tibble(
  phone = c("p/b", "p/b"),
  language = c("Tseltal", "Yélî Dnye"),
  mean_d = c(
    mean(subset(pb.by.site.age, language == "Tseltal")$mean.delta),
    mean(subset(pb.by.site.age, language == "Yélî Dnye")$mean.delta))
)

deltas.td <- tibble(
  phone = c("t/d", "t/d"),
  language = c("Tseltal", "Yélî Dnye"),
  mean_d = c(
    mean(subset(td.by.site.age, language == "Tseltal")$mean.delta),
    mean(subset(td.by.site.age, language == "Yélî Dnye")$mean.delta))
)

deltas.kg <- tibble(
  phone = c("k/g", "k/g"),
  language = c("Tseltal", "Yélî Dnye"),
  mean_d = c(
    mean(subset(kg.by.site.age, language == "Tseltal")$mean.delta),
    mean(subset(kg.by.site.age, language == "Yélî Dnye")$mean.delta))
)

deltas.m <- tibble(
  phone = c("m", "m"),
  language = c("Tseltal", "Yélî Dnye"),
  mean_d = c(
    mean(subset(m.by.site.age, language == "Tseltal")$mean.delta),
    mean(subset(m.by.site.age, language == "Yélî Dnye")$mean.delta))
)

deltas.n <- tibble(
  phone = c("n", "n"),
  language = c("Tseltal", "Yélî Dnye"),
  mean_d = c(
    mean(subset(n.by.site.age, language == "Tseltal")$mean.delta),
    mean(subset(n.by.site.age, language == "Yélî Dnye")$mean.delta))
)

deltas.l <- tibble(
  phone = c("l", "l"),
  language = c("Tseltal", "Yélî Dnye"),
  mean_d = c(
    mean(subset(l.by.site.age, language == "Tseltal")$mean.delta),
    mean(subset(l.by.site.age, language == "Yélî Dnye")$mean.delta))
)

all.deltas <- bind_rows(
  deltas.pb, deltas.td, deltas.kg,
  deltas.m, deltas.n, deltas.l) %>%
  mutate(mean_d_round = round(mean_d, 2))

round(mean(pb.by.site.age$mean.delta),2)
round(mean(td.by.site.age$mean.delta),2)
round(mean(kg.by.site.age$mean.delta),2)
round(mean(m.by.site.age$mean.delta),2)
round(mean(n.by.site.age$mean.delta),2)
round(mean(l.by.site.age$mean.delta),2)
sum(pb.by.site.age$n.vocs)
sum(td.by.site.age$n.vocs)
sum(kg.by.site.age$n.vocs)
sum(m.by.site.age$n.vocs)
sum(n.by.site.age$n.vocs)
sum(l.by.site.age$n.vocs)

ggplot(all.deltas, aes(y = mean_d, x = phone,
                       color = language, shape = language)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  labs(title = "Mean difference in token count within a vocalization",
       x = "Phone type",
       y = "Mean difference")
