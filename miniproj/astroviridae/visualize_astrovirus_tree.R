setwd("C:/MALICHAVA/Virome/picorna-virome/miniproj/astroviridae")

library(treeio)
library(ggtree)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)

# ---------- 1. Читаем дерево ----------
tree <- read.newick("centroids_all_bats_and_our_contigs_complete_1B_aln_removed_4_bad_seqs.faa.treefile")

# ---------- 2. Набор родов рукокрылых ----------
BAT_GENERA <- c(
  "Miniopterus", "Nyctalus", "Eidolon", "Rousettus",
  "Myotis", "Hipposideros", "Rhinolophus", "Scotophilus",
  "Pipistrellus", "Tylonycteris", "Mops", "Coleura",
  "Pteropus", "Eptesicus", "Tadarida", "Molossus", "Carollia",
  "Artibeus", "Desmodus", "Glossophaga", "Anoura", "Megaderma",
  "Rhinopoma", "Mormopterus", "Chaerephon", "Nycteris",
  "Emballonura", "Taphozous", "Saccolaimus", "Chalinolobus",
  "Lasiurus", "Nycticeius", "Plecotus", "Barbastella", "Murina",
  "Kerivoula", "Phoniscus", "Cynopterus", "Macroglossus",
  "Eonycteris", "Epomophorus", "Epomops", "Hypsignathus",
  "Micropteropus", "Myonycteris", "Lissonycteris"
)

# ---------- 3. Парсер одного типла ----------
# Возвращает data.frame с полями: label, host_raw, host_genus, is_contig, category
parse_tip <- function(label) {
  # Контиг? Проверяем по наличию "|" — если нет, то это контиг
  is_contig <- !grepl("\\|", label)
  
  host_raw <- NA_character_
  
  if (is_contig) {
    # SRR28250970/k141_146080/Kenya/2015/Mops-condylurus
    parts <- str_split(label, "/", simplify = TRUE)
    if (length(parts) >= 5) host_raw <- parts[5]
  } else {
    # MN832787|DNK|Myotis-daubentonii|2014|NA
    parts <- str_split(label, "\\|", simplify = TRUE)
    if (length(parts) >= 3) host_raw <- parts[3]
  }
  
  # "nan", "NA", "", "viral-metagenome" и т.п. -> неизвестно
  if (is.na(host_raw) ||
      host_raw %in% c("NA", "nan", "", "viral-metagenome", "unclassified")) {
    return(data.frame(
      label      = label,
      host_raw   = host_raw,
      host_genus = NA_character_,
      is_contig  = is_contig,
      category   = "unknown",
      stringsAsFactors = FALSE
    ))
  }
  
  # Род = первое слово, разделённое дефисом (для bat-формата "Mops-condylurus")
  # или для вирусного формата "Myotis-daubentonii" тоже первое слово до "-"
  host_genus <- str_split(host_raw, "[-_ ]", simplify = TRUE)[1, 1]
  
  # Проверяем вхождение в список рукокрылых
  is_bat <- host_genus %in% BAT_GENERA
  
  category <- if (is_bat) "bat" else "other_mammal"
  
  data.frame(
    label      = label,
    host_raw   = host_raw,
    host_genus = host_genus,
    is_contig  = is_contig,
    category   = category,
    stringsAsFactors = FALSE
  )
}

# ---------- 4. Прогоняем по всем листьям ----------
tip_df <- data.frame(label = tree@phylo$tip.label, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(parsed = list(parse_tip(label))) %>%
  unnest(parsed) %>%
  select(-label, label = label)

# ---------- 5. Цвета ----------
tip_colors <- c(
  bat          = "#1f78b4",  # синий — летучие мыши
  other_mammal = "#33a02c",  # зелёный — прочие млекопитающие
  unknown      = "black"     # чёрный — хост не определён
)

# ---------- 6. Базовый plot ----------
p <- ggtree(tree, branch.length = "none") %<+% tip_df +
  geom_tiplab(aes(color = category), size = 2.2, align = TRUE,
              offset = 0.05, show.legend = FALSE) +
  scale_color_manual(values = tip_colors, na.value = "black")

# ---------- 7. Bootstrap-точки (>95) на узлах ----------
# В newick bootstrap-значения читаются treeio как node.label
# (или как label, зависит от того, чем вы записали — nhx/nwk)
p <- p +
  geom_nodepoint(
    aes(subset = !is.na(as.numeric(label)) & as.numeric(label) >= 95),
    size = 1.2, color = "black", shape = 16
  )

# ---------- 8. (опционально) подсветка контигов ----------
# Если хотите — можно чуть подсветить подписи контигов, например жирным
p <- p + geom_tiplab(
  aes(subset = is_contig, color = category),
  fontface = "bold", size = 2.2, align = TRUE,
  offset = 0.05, show.legend = FALSE
)

p + theme(legend.position = "none")