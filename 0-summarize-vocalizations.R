# library(tidyverse)
vowels <- "[iɪeɛəoɔaɑu ]"

metadata <- read_delim("ACLEW_list_of_corpora-CasillasOnly-UTF8.csv",
                       delim = ";",
                       locale = locale(encoding = "UTF-8"))

files <- list.files(".", "\\d{4}.txt", recursive = TRUE)
all.vocs.tbl <- tibble()
all.consonant.freq <- tibble()
all.consonant.freq.clip <- tibble()

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
  data <- read_delim(filename, delim = "\t", col_names = FALSE)
  names(data) <- c("tier", "deleteme", "start", "stop", "length", "value")
  data$deleteme <- NULL
  
  # make sure we're only including vocs from random clips
  clips <- filter(data, tier == "code" & grepl("random", value))
  random.clip.vocs <- tibble()
  for (i in 1:nrow(clips)) {
    vocs <- which(data$tier == "pho@CHI" &
        ((data$start >= clips$start[i] & data$stop <= clips$stop[i])|
        (data$start < clips$start[i] & data$stop <= clips$stop[i] &
            data$stop > clips$start[i])|
        (data$start >= clips$start[i] & data$start <= clips$stop[i] &
            data$stop > clips$stop[i])))
    vocs.tbl <- bind_cols(data[vocs,],
      tibble(clip = clips$value[i], language = language, filename = filename))
    random.clip.vocs <- bind_rows(random.clip.vocs, vocs.tbl)
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

  # calculate canonical proportion (C/(N+C))
  # includes ALL consonants (i.e., no exclusions)
  noncanon.idx <- which(random.clip.vocs$value == "N")
  canon.idx <- which(random.clip.vocs$value != "N" &
      random.clip.vocs$value != "Y" &
      random.clip.vocs$value != "L")
  canon.prop <- length(canon.idx)/
    (length(noncanon.idx) + length(canon.idx))
  
  # extract consonant frequencies
  all.babble.string <- paste(
    random.clip.vocs$value[canon.idx], collapse = "")
  all.babble.string <- gsub(vowels, "", all.babble.string)
  
  consonants <- unlist(strsplit(all.babble.string, ''))
  consonant.freq <- tibble(
    character = consonants) %>%
    group_by(character) %>%
    summarize(raw_count = n()) %>%
    arrange(-raw_count)
  
  
  # percentages of consonant use
  total.consonants <- sum(consonant.freq$raw_count)
  consonant.freq  <- consonant.freq %>%
    mutate(percent.use = (raw_count/total.consonants)*100)
  consonant.freq$canon.prop <- canon.prop
  consonant.freq$aclew_id <- as.numeric(str_extract(filename, "\\d{4}"))
  consonant.freq <- consonant.freq %>%
    left_join(metadata)
  all.consonant.freq <- bind_rows(all.consonant.freq, consonant.freq)
  
  # extract consonant frequencies again, but this time by clip
  # and don't bother with percentages, canon.prop, etc.
  for (i in 1:nrow(clips)) {
    all.babble.string <- paste(
      subset(random.clip.vocs[canon.idx,], clip == clips$value[i])$value,
      collapse = "")
    all.babble.string <- gsub(vowels, "", all.babble.string)
    consonants <- unlist(strsplit(all.babble.string, ''))
    consonant.freq <- tibble(
      character = consonants) %>%
      group_by(character) %>%
      summarize(raw_count = n()) %>%
      mutate(clip = clips$value[i]) %>%
      arrange(-raw_count)
    consonant.freq$aclew_id <- as.numeric(str_extract(filename, "\\d{4}"))
    consonant.freq <- consonant.freq %>%
      left_join(metadata)
    all.consonant.freq.clip <- bind_rows(all.consonant.freq.clip,
                                         consonant.freq)
  }
}

write_csv(all.vocs.tbl, "all.vocalizations.csv")
write_csv(all.consonant.freq, "all.consonant.freq.csv")
write_csv(all.consonant.freq.clip, "all.consonant.freq.clip.csv")

