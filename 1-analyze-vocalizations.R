library(tidyverse)
library(jtools)
library(lme4)

# segmental restrictions
include.affricates <- FALSE
include.glides <- FALSE
include.glottals <- FALSE

# IPA substitutions
# tʃ' goes to $
# tʃ goes to S
# dʒ goes to Z
# p' goes to P
# t' goes to T
# k' goes to K
# multi-character YD symbols aren't otherwise used afaik

excl.segments <- c()
if (include.affricates == FALSE) {
  excl.segments <- c(excl.segments, c("$", "S", "Z"))
}
if (include.glides == FALSE) {
  excl.segments <- c(excl.segments, c("j", "w"))
}
if (include.glottals == FALSE) {
  excl.segments <- c(excl.segments, c("h", "ʔ"))
}

# minimum threshold for vms
min.num.tokens.vms.ts <- 10
min.num.tokens.vms.yd <- 5

source("0-summarize-vocalizations.R")
source("0-vihman-tabler.R")

#---- MAIN

all_data <- read_csv("all.consonant.freq.csv")
n_vocs_rec <- read_csv("all.vocalizations.csv") %>%
  group_by(filename) %>%
  summarize(total.n.vocs = n()) %>%
  mutate(aclew_id = as.numeric(str_extract(filename, "\\d{4}"))) %>%
  select(-filename)
consonants <- read_delim("consonant_table-UTF8.txt", delim = "\t")

all_data <- all_data %>%
  left_join(consonants, by = c("character" = "consonant")) %>%
  filter(language == "Tseltal" | language == "Yélî Dnye")
# sanity check: make sure everything has a hyper_place
# which(is.na(all_data$hyper_place))

# Counts of consonant uses (Laing and Bergelson consonants only)
LB_Cs <- c("p", "b", "t", "d", "k", "g", "m", "n", "s", "l", "w", "r", "j")
LB_C_subset <- all_data %>%
  filter(character %in% LB_Cs) %>%
  group_by(aclew_id, language) %>%
  summarize(tot.C.count = sum(raw_count))


### SUMMARY TABLES
# Counts how many instances of each segment there is for each child and
# only shows counts for segments that pass the VMS threshold

# Write out the Tseltal segmental summary tables
create.vihman.table(all_data, "Tseltal", min.num.tokens.vms.ts,
                    include.affricates, include.glides, include.glottals)

# Write out the Yélî Dnye segmental summary tables
create.vihman.table(all_data, "Yélî Dnye", min.num.tokens.vms.yd,
                    include.affricates, include.glides, include.glottals)


### CANONICAL PROPORTION
cp.plot.data <- all_data %>%
  select(aclew_id, age_mo_round, canon.prop, language) %>%
  distinct() %>%
  left_join(n_vocs_rec) %>%
  mutate(
    threshold = case_when(
      total.n.vocs >= 100 ~ "Yes",
      total.n.vocs < 100 ~ "No",
      TRUE ~ "EEP"
    ))

ggplot(data = cp.plot.data,
       aes(x = age_mo_round, y = canon.prop, size = total.n.vocs,
           color = language, fill = language)) +
  geom_smooth(method = "lm")


canon.prop_plot <- ggplot(data = cp.plot.data,
  aes(x = age_mo_round, y = canon.prop, size = total.n.vocs,
      color = language, fill = language)) +
  geom_smooth(method = "lm") +
  # facet_grid(~ language) +
  geom_point(aes(shape = threshold), alpha = 0.6) +
  labs(x = "Age (months)" , y = "Canonical proportion")+
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(5,20)) +
  scale_size("Vocalizations", range = c(0, 6)) +
  scale_shape_manual(values = c(3, 19)) +
  scale_fill_manual("Language", values = c("cyan3", "darkgoldenrod1")) +
  scale_colour_manual("Language", values = c("cyan3", "darkgoldenrod1")) +
  # scale_colour_manual("Threshold", values = c("red", "gray20")) +
  # geom_segment(aes(y = 0.15, yend = 0.15, x = 10, xend = 20),
  #              linetype = "dashed") +
  annotate("rect", xmin = 10, xmax = 20, ymin = 0, ymax = 0.15,
           color = "black", fill = "white", alpha = 0.3) +
  theme_apa(legend.use.title = TRUE,
            legend.font.size = 10) +
  guides(shape = "none",
         size = "none")
ggsave("canon.prop_plot.png", device = png(),
       width = 20, height = 10, units = "cm")
print(canon.prop_plot)
dev.off()

### VMS: INDIVIDUAL SEGMENTS

ts.vms <- read_csv(
  "vihman-style-vms-table-Tseltal-post_segment_exclusions.csv")
yd.vms <- read_csv(
  "vihman-style-vms-table-Yélî Dnye-post_segment_exclusions.csv")
vms.by.age <- bind_rows(ts.vms,yd.vms) %>%
  left_join(select(metadata, c(aclew_id, language)))

vms.plot.data <- vms.by.age %>%
  mutate(
    threshold = case_when(
      total.n.vocs >= 100 ~ "Yes",
      total.n.vocs < 100 ~ "No",
      TRUE ~ "EEP"
    ))

vms.by.age_plot <- ggplot(data = vms.plot.data,
  aes(x = age_mo_round, y = vms.count, size = total.n.vocs,
      fill = language, color = language)) +
  # facet_grid(~ language) +
  # geom_segment(aes(x = 9, xend = 20, y = 2, yend = 2), linetype = "dashed") +
  annotate("rect", xmin = 12, xmax = 20, ymin = 0, ymax = 2,
           color = "black", fill = "white", alpha = 0.3) +
  # facet_grid(~ language) +
  geom_smooth(method = "lm") +
  geom_point(aes(shape = threshold), alpha = 0.6) +
  labs(x = "Age (months)" , y = "# vocal motor schemes") +
  scale_y_continuous(limits = c(0,8)) +
  scale_x_continuous(limits = c(5,20)) +
  scale_size("Vocalizations", range = c(0, 6)) +
  scale_shape_manual(values = c(3, 19)) +
  scale_fill_manual("Language", values = c("cyan3", "darkgoldenrod1")) +
  scale_colour_manual("Language", values = c("cyan3", "darkgoldenrod1")) +
  # geom_point(aes(color = threshold), alpha = 0.4, position = position_dodge()) +
  # geom_jitter(aes(color = threshold), width = 0.2, height = 0, alpha = 0.4) +
  theme_apa(legend.use.title = TRUE,
            legend.font.size = 10) +
  guides(shape = "none",
         size = "none")
ggsave("vms.by.age_plot.png", device = png(),
       width = 20, height = 10, units = "cm")
print(vms.by.age_plot)
dev.off()

all.consonant.counts <- all_data %>%
  group_by(language, character) %>%
  summarize(total_count = sum(raw_count))
all.consonant.counts$character2 <- factor(
  all.consonant.counts$character)
all.consonant.counts$character2 <- factor(
  all.consonant.counts$character2,
  labels = c(
    "b", "d", "f", "ɡ", "ɣ", "h", "j", "k", "l", "m", "n", "ŋ",
    "p", "r", "s", "tʃ", "ʃ", "t", "w", "x", "z", "dʒ", "ʒ", "ʔ"))
all.consonant.counts$character2 <- factor(
  all.consonant.counts$character2,
  levels = c(
   "p", "b", "t", "d", "k", "ɡ", "m", "n", "ŋ",
   "w", "j", "l", "r",
   "tʃ", "dʒ", "f", "s", "z", "ʃ", "ʒ", "x", "ɣ",
   "ʔ", "h"))
all.consonants.used <- ggplot(data = all.consonant.counts,
                            aes(x = character2, y = total_count)) +
  geom_col() +
  facet_wrap(~ language, ncol = 1) +
  labs(x = "Phone", y = "Total # instances") +
  theme_apa()
ggsave("all.phone.counts.png", device = png(),
       width = 24, height = 10, units = "cm")
print(all.consonants.used)
dev.off()

### VMS: HYPER PLACE

# all_places_all_ptcps <- expand_grid(aclew_id = unique(all_data$aclew_id),
#   hyper_place = unique(all_data$hyper_place)) %>%
#   left_join(metadata)
# 
# all_data_hyper_place <- all_data %>%
#   group_by(aclew_id, age_mo_round, hyper_place, language) %>%
#   summarize(percent.hp.use = sum(percent.use)) %>%
#   full_join(all_places_all_ptcps) %>%
#   replace_na(list(percent.hp.use = 0))
# all_data_hyper_place$hyper_place <- factor(all_data_hyper_place$hyper_place,
#   levels = c("labial", "coronal", "dorsal", "laryngeal"))
# all_data_hyper_place$hyper_place <- factor(all_data_hyper_place$hyper_place,
#   labels = c("Labial", "Coronal", "Dorsal", "Laryngeal"))
# 
# hyper_place_plot <- ggplot(
#   data = all_data_hyper_place ,
#   aes(x = age_mo_round, y = percent.hp.use, color = hyper_place)) +
#   geom_smooth(method = "lm") +
#   geom_jitter() +
#   labs(x = "Age (months)", y = "% Consonant productions") +
#   coord_cartesian(ylim=c(0,100), xlim=c(0,20)) +
#   facet_grid(vars(language), vars(hyper_place)) +
#   theme_apa() +
#   theme(legend.position = "none")
# ggsave("hyper_place_plot.png", device = png(),
#        width = 20, height = 10, units = "cm")
# print(hyper_place_plot)
# dev.off()
# 
# # Plot of vms by age
# collapsed.segment.freq.all <- bind_rows(collapsed.segment.freq.table.all.ts,
#   collapsed.segment.freq.table.all.yd) %>%
#   left_join(select(metadata, c(aclew_id, language)))
# vms.by.age_plot <- ggplot(data = collapsed.segment.freq.all,
#   aes(x = age_mo_round, y = vms.count)) +
#   facet_grid(~ language) +
#   geom_hline(yintercept = 2, lty = "dashed") +
#   geom_smooth(method = "lm", color = "black") +
#   geom_jitter(color = "red") +
#   labs(x = "Age (months)" , y = "# VMS by hyperplace")+
#   scale_y_continuous(limits = c(0,8))+
#   scale_x_continuous(limits = c(5,20))+
#   theme_apa()
# ggsave("vms.by.age_plot.png", device = png(),
#        width = 20, height = 10, units = "cm")
# print(vms.by.age_plot)
# dev.off()
# 
# # Plot of vms-hyperplace
# hyperplace.freq.all <- bind_rows(hyperplace.freq.table.all.ts,
#   hyperplace.freq.table.all.yd) %>%
#   left_join(select(metadata, c(aclew_id, language)))
# vms.by.age.hyperplace_plot <- ggplot(data = hyperplace.freq.all,
#   aes(x = age_mo_round, y = vms.count)) +
#   facet_grid(~ language) +
#   geom_smooth(method = "lm", color = "black") +
#   geom_jitter(color = "red") +
#   labs(x = "Age (months)" , y = "# VMS by hyperplace")+
#   scale_y_continuous(limits = c(0,5))+
#   scale_x_continuous(limits = c(5,20))+
#   theme_apa()
# ggsave("vms.by.age.hyperplace_plot.png", device = png(),
#        width = 20, height = 10, units = "cm")
# print(vms.by.age.hyperplace_plot)
# dev.off()

# statistical models

## CP
model.cp <- lm(canon.prop ~ age_mo_round * language,
                data = cp.plot.data)

model.cp.tbl <- tidy(model.cp) %>%
  mutate(estimate = round(estimate, 3),
         std.error = round(std.error, 3),
         statistic = round(statistic, 3),
         p.value = case_when(
           p.value < .001 ~ "<.001",
           p.value < .01 ~ "<.01",
           p.value < .05 ~ "<.05",
           TRUE ~ as.character(round(p.value, 3))
         ))
write_delim(model.cp.tbl, "cp.regression.estimates.txt", delim = "\t")

## VMS
model.vms <- lm(vms.count ~ age_mo_round * language,
                data = vms.by.age)
# model.vms.hyperplace <- lm(vms.count ~ age_mo_round * language,
#                            data = hyperplace.freq.all)

model.vms.tbl <- tidy(model.vms) %>%
  mutate(estimate = round(estimate, 3),
         std.error = round(std.error, 3),
         statistic = round(statistic, 3),
         p.value = case_when(
           p.value < .001 ~ "<.001",
           p.value < .01 ~ "<.01",
           p.value < .05 ~ "<.05",
           TRUE ~ as.character(round(p.value, 3))
         ))
write_delim(model.vms.tbl, "vms.regression.estimates.txt", delim = "\t")


# Alternative phone stability analysis
crossclip.data <- read_csv("all.consonant.freq.clip.csv")

clips.per.rec.per.phone <- crossclip.data %>%
  mutate(phone = case_when(
    grepl("[pb]", character) ~ "p/b",
    grepl("[td]", character) ~ "t/d",
    grepl("[kɡ]", character) ~ "k/g",
    TRUE ~ character
  )) %>%
  group_by(aclew_id, language, phone) %>%
  summarize(n_clips = n())

usable.phones <- c("p/b", "t/d", "k/g", "m", "n", "l", "s")
multiclip.phones <- clips.per.rec.per.phone %>%
  filter(phone %in% usable.phones) %>%
  filter(n_clips >= 4)

nmulticlipphones.perkid <- multiclip.phones %>%
  group_by(aclew_id, language) %>%
  summarize(n_multiclip4 = n())

all.measures <- vms.by.age %>%
  left_join(nmulticlipphones.perkid) %>%
  replace_na(list(vms.count = 0, n_multiclip4 = 0))

write_csv(all.measures, "all.measures.table.csv")

vms.vs.mclip <- ggplot(all.measures, aes(x = vms.count, y = n_multiclip4,
                                         color = language, fill = language,
                                         size = total.n.vocs)) +
  geom_smooth(method = "lm") +
  geom_jitter(width = 0.1, height = 0.1) +
  labs(x = "VMS count (low threshold)", y = "# 4+ cross-clip consonants") +
  # coord_cartesian(xlim = c(-1, 8), ylim = c(-1, 8)) +
  # scale_y_continuous(limits = c(0,7)) +
  # scale_x_continuous(limits = c(0,7)) +
  scale_size("Vocalizations", range = c(0, 6)) +
  scale_shape_manual(values = c(3, 19)) +
  scale_fill_manual("Language", values = c("cyan3", "darkgoldenrod1")) +
  scale_colour_manual("Language", values = c("cyan3", "darkgoldenrod1")) +
  # geom_point(aes(color = threshold), alpha = 0.4, position = position_dodge()) +
  # geom_jitter(aes(color = threshold), width = 0.2, height = 0, alpha = 0.4) +
  theme_apa(legend.use.title = TRUE,
            legend.font.size = 10) +
  guides(shape = "none",
         size = "none")
ggsave("vms.vs.multiclip.png", device = png(),
       width = 20, height = 10, units = "cm")
print(vms.vs.mclip)
dev.off()

multiclip.by.age_plot <- ggplot(data = all.measures,
                          aes(x = age_mo_round, y = n_multiclip4,
                              size = total.n.vocs,
                              color = language, fill = language)) +
  # facet_grid(~ language) +
  # geom_segment(aes(x = 9, xend = 20, y = 2, yend = 2), linetype = "dashed") +
  annotate("rect", xmin = 12, xmax = 20, ymin = 0, ymax = 2,
           color = "black", fill = "white", alpha = 0.3) +
  geom_smooth(method = "lm") +
  geom_point(alpha = 0.6) +
  labs(x = "Age (months)" , y = "# 4+ cross-clip consonants") +
  scale_y_continuous(limits = c(0,8)) +
  scale_x_continuous(limits = c(5,20)) +
  scale_size("Vocalizations", range = c(0, 6)) +
  scale_shape_manual(values = c(3, 19)) +
  scale_fill_manual("Language", values = c("cyan3", "darkgoldenrod1")) +
  scale_colour_manual("Language", values = c("cyan3", "darkgoldenrod1")) +
  # geom_point(aes(color = threshold), alpha = 0.4, position = position_dodge()) +
  # geom_jitter(aes(color = threshold), width = 0.2, height = 0, alpha = 0.4) +
  theme_apa(legend.use.title = TRUE,
            legend.font.size = 10) +
  guides(shape = "none",
         size = "none")
ggsave("multiclip.by.age_plot.png", device = png(),
       width = 20, height = 10, units = "cm")
print(multiclip.by.age_plot)
dev.off()

model.mclp <- lm(n_multiclip4 ~ age_mo_round * language,
                data = all.measures)
model.mclp.tbl <- tidy(model.mclp) %>%
  mutate(estimate = round(estimate, 3),
         std.error = round(std.error, 3),
         statistic = round(statistic, 3),
         p.value = case_when(
           p.value < .001 ~ "<.001",
           p.value < .01 ~ "<.01",
           p.value < .05 ~ "<.05",
           TRUE ~ as.character(round(p.value, 3))
         ))
write_delim(model.mclp.tbl, "mclp.regression.estimates.txt", delim = "\t")

cor.test(all.measures$vms.count, all.measures$n_multiclip4)
