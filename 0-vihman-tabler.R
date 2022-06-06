create.vihman.table <- function(dataset, lg, min.tokens,
                                inc.aff, inc.gli, inc.glo) {
  
  # To create these tables we want to know how many tokens
  # each child produced as well as their age, so we'll first
  # make a table with that information so that we can integrate
  # the information back in later on
  table.ptcp.info <- dataset %>%
    # This says for each child and age...
    group_by(aclew_id, age_mo_round, language) %>%
    # give me a column that simply sums the
    # raw counts
    summarize(total.tokens = sum(raw_count))

  # First create a new column with the collapsed places
  dataset <- dataset %>%
    mutate(
      # case_when() is a useful function that in this instance
      # lets us give different outputs for a buch of different
      # scenarios in mapping 'character' to our new column values
      collapsed.segments = case_when(
        character == "p" | character == "b" ~ "p/b",
        character == "t" | character == "d" ~ "t/d",
        character == "S" | character == "Z" ~ "tʃ/dʒ",
        character == "ʃ" | character == "ʒ" ~ "ʃ/ʒ",
        character == "k" | character == "ɡ" ~ "k/ɡ",
        character == "s" | character == "z" ~ "s/Z",
        character == "x" | character == "ɣ" ~ "x/ɣ",
        TRUE ~ character))
  # And a vector with the consonants ordered similarly to
  # the Vihman display (approx front-back)
  consonant.order <- c(
    # plosives, labial to velar
    "p/b", "t/d", "k/ɡ",
    # other non-glide labials
    "m", "f",
    # other non-glide coronals
    "n", "l", "r", "s/z",
    # post-alveolar
    "tʃ/dʒ", "ʃ/ʒ",
    # other non-glide velar
    "ŋ", "x/ɣ",
    # glottals and glides
    "w", "j", "ʔ", "h"
  )
  
  # so we first make sure that the factor levels match the
  # order we want to use
  dataset$collapsed.segments <- factor(
    dataset$collapsed.segments, levels = consonant.order
  )
  
  collapsed.segment.freq.table <- dataset %>%
    filter(language == lg) %>%
    # We group by child and segment
    group_by(aclew_id, collapsed.segments) %>%
    # Then we get the sum of raw counts in each of those groups
    # which is equivalent to counting the number of instances
    # of each collapsed.segment type for each recording
    summarize(n.instances.segment = sum(raw_count)) %>%
    # We then only keep cases where children produced 
    # at least the 'min.num.tokens.vms' number of tokens
    # Note you can set 'min.num.tokens.vms' at the top of
    # this script
    filter(n.instances.segment >= min.tokens) %>%
    # Since we want to make a table with all possible cases and simply
    # put NAs for cases where children didn't reach vms for that segment
    # we join with a tibble that has all possible id * segment combinations
    full_join(filter(dataset, language == lg) %>%
                tidyr::expand(aclew_id, collapsed.segments)) %>%
    # we flip the table around so that the segments are the column names
    pivot_wider(
      names_from = collapsed.segments,
      values_from = n.instances.segment) %>%
    ungroup() %>%
    # we add the participant information back in
    left_join(table.ptcp.info, by = "aclew_id") %>%
    # we sort by age
    arrange(age_mo_round)
  
  # Now in order to make the vms table the way she has it we split into
  # affricates
  collapsed.segment.freq.table.affricates <- 
    collapsed.segment.freq.table %>%
    select("tʃ/dʒ")
  # glides
  collapsed.segment.freq.table.glides <- 
    collapsed.segment.freq.table %>%
    select("w", "j")
  # glottals
  collapsed.segment.freq.table.glottals <- 
    collapsed.segment.freq.table %>%
    select("ʔ", "h")
  # and everything else
  collapsed.segment.freq.table.nonglottalsglidesaffricates <- 
    collapsed.segment.freq.table %>%
    select(-"w", -"j", -"ʔ", -"h", -"tʃ/dʒ") %>%
    mutate(
      # we count all the non-NA cases in the "everything else" segments
      vms.count = rowSums(!is.na(.)) - 4
    )
  
  # Then we create versions with and without the excluded segments
  collapsed.segment.freq.table.all.allsegs <- bind_cols(
    collapsed.segment.freq.table.nonglottalsglidesaffricates,
    collapsed.segment.freq.table.affricates) %>%
    bind_cols(collapsed.segment.freq.table.glides) %>%
    bind_cols(collapsed.segment.freq.table.glottals) %>%
    left_join(n_vocs_rec) %>%
    rename("total.cons.tokens" = total.tokens) %>%
    # and we put the columns in the desired order
    select(aclew_id, age_mo_round,
           "p/b", "t/d", "k/ɡ", "m", "f",
           "n", "l", "r", "s/z", "ʃ/ʒ", "ŋ", "x/ɣ",
           vms.count, total.cons.tokens, total.n.vocs,
           "tʃ/dʒ", "w", "j", "ʔ", "h")
  
  collapsed.segment.freq.table.all.withexcls <- 
    collapsed.segment.freq.table.all.allsegs
  if (inc.aff == FALSE) {
    collapsed.segment.freq.table.all.withexcls <- 
      collapsed.segment.freq.table.all.withexcls %>%
      select(-"tʃ/dʒ")
  }
  if (inc.gli == FALSE) {
    collapsed.segment.freq.table.all.withexcls <- 
      collapsed.segment.freq.table.all.withexcls %>%
      select(-"w", -"j")
  }
  if (inc.glo == FALSE) {
    collapsed.segment.freq.table.all.withexcls <- 
      collapsed.segment.freq.table.all.withexcls %>%
      select(-"ʔ", -"h")
  }
  
  # Lastly, we write the table to a csv and a tab-delimited text file
  # with an automated naming scheme
  # All segments
  write_delim(collapsed.segment.freq.table.all.allsegs,
              paste0("vihman-style-vms-table-", lg, "-all_segments.txt"),
              delim = "\t")
  write_csv(collapsed.segment.freq.table.all.allsegs,
              paste0("vihman-style-vms-table-", lg, "-all_segments.csv"))
  # Minus excluded segments
  write_delim(collapsed.segment.freq.table.all.withexcls,
              paste0("vihman-style-vms-table-", lg,
                     "-post_segment_exclusions.txt"),
              delim = "\t")
  write_csv(collapsed.segment.freq.table.all.withexcls,
            paste0("vihman-style-vms-table-", lg,
                   "-post_segment_exclusions.csv"))
  
}