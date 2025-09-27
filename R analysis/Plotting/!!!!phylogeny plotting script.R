library(ggtree)
library(ggplot2)
library(ggnewscale)
library(dplyr)

TREE_FILE <- "FcC_supermatrix.coyote1.rooted.fixed.tre"
META_FILE <- "full_metadata.csv"

tree <- read.tree(TREE_FILE)
meta <- read.csv(META_FILE, stringsAsFactors = FALSE)

# Optional: enforce factor orders
meta$taxon  <- factor(meta$taxon,  levels = c("Dog","Wolf"))
meta$period <- factor(meta$period, levels = c("Fossil","Modern"))

bioGeo_order  <- c("Alaska","Canada","Central_America","Central_Asia","East_Asia",
                   "Europe","Greenland","Middle_East","Russia","Siberia",
                   "Southwest_Asia","USA")
bioGeo_colors <- c(
  # Arctic / Eurasia → Americas gradient
  "Siberia"        = "#245E90",  # darker denim blue
  "Russia"         = "#6BAED3",  # lighter steel/azure
  "Alaska"         = "#2FA699",  # deeper teal
  "Canada"         = "#4FA45D",  # cooler medium green
  "USA"            = "#839E2E",  # olive (slightly darker)
  
  # Arctic adjunct (much lighter for separation)
  "Greenland"      = "#9FD0F5",  # icy light blue
  
  # Asia (clear separation; darker vs lighter pair)
  "East_Asia"      = "#6A56A5",  # deep muted purple
  "Central_Asia"   = "#B789D0",  # lighter lilac
  
  # Europe → Middle East → SW Asia (warm earth tones)
  "Europe"         = "#B8852F",  # darker ochre
  "Middle_East"    = "#8F4A2A",  # deep brick/rust
  "Southwest_Asia" = "#B8684A",  # muted terracotta
  
  # Americas south (deeper coral, not neon)
  "Central_America"= "#B74747"
)
meta$macro_area <- factor(meta$macro_area, levels = bioGeo_order)
pdf("Tree.pdf", width = 40, height = 60)


# tree + meta already loaded; meta must have columns: label (or tip), taxon, period, macro_area
# If your join key is `tip`, mirror it to `label` so ggtree can use it for tiplab:
# meta <- meta %>% mutate(label = tip)

p <- ggtree(tree, size = 0.3) %<+% meta +
  theme_tree2() +
  
  # tip labels (from the tree's label column after %<+%)
  geom_tiplab(aes(label = label), size = 2.2, hjust = 0) +
  
  # --- TIP POINTS: colour = Dog/Wolf, shape = Fossil/Modern ---
  geom_tippoint(aes(color = taxon, shape = period), size = 1.8) +
  scale_color_manual(values = c(Dog = "#ff7f00", Wolf = "#377eb8"), name = "Taxon") +
  scale_shape_manual(values = c(Fossil = 16, Modern = 17), name = "Period") +
  
  # bootstrap (only ≥ 70)
  geom_text2(
    aes(subset = !isTip & as.numeric(label) >= 70, label = label),
    hjust = -0.3, size = 2.2, color = "grey30"
  )

geom_tippoint(
  data = subset(p$data, label == "coyote1"),
  shape = 8,
  size = 4,
  color = "gold"
) 

# Ensure meta has label matching tree$tip.label
meta <- meta %>% mutate(label = tip)  # or label = tip if needed

# Reorder to match the tree's tip order
meta_ord <- tibble(label = tree$tip.label) %>% left_join(meta, by = "label")

# Build heatmap df: rownames = tip labels; one column 'macro_area'
hm <- data.frame(macro_area = meta_ord$macro_area)
rownames(hm) <- meta_ord$label

# Add a new fill scale and the heatmap strip
bioGeo_colors <- c(
  # Arctic / Eurasia → Americas gradient
  "Siberia"        = "#245E90",  # darker denim blue
  "Russia"         = "#6BAED3",  # lighter steel/azure
  "Alaska"         = "#2FA699",  # deeper teal
  "Canada"         = "#4FA45D",  # cooler medium green
  "USA"            = "#839E2E",  # olive (slightly darker)
  
  # Arctic adjunct (much lighter for separation)
  "Greenland"      = "#9FD0F5",  # icy light blue
  
  # Asia (clear separation; darker vs lighter pair)
  "East_Asia"      = "#6A56A5",  # deep muted purple
  "Central_Asia"   = "#B789D0",  # lighter lilac
  
  # Europe → Middle East → SW Asia (warm earth tones)
  "Europe"         = "#B8852F",  # darker ochre
  "Middle_East"    = "#8F4A2A",  # deep brick/rust
  "Southwest_Asia" = "#B8684A",  # muted terracotta
  
  # Americas south (deeper coral, not neon)
  "Central_America"= "#B74747"
)

p2 <- p +
  geom_tippoint(
    data = subset(p$data, label == "coyote1"),
    shape = 8,               # star shape
    size = 4,
    color = "gold"
  )

p2 <- new_scale_fill()
p2 <- gheatmap(p, hm, width = 0.06, colnames = FALSE, color = NA, offset = 0.01) +
  scale_fill_manual(values = bioGeo_colors, name = "Macro-area", drop = FALSE)

# Yellow star for coyote



pdf("Tree.pdf", width = 20, height = 50, onefile = FALSE)
print(p2)

dev.off()

if (.Platform$OS.type == "windows") shell.exec("Tree.pdf")

# ===== 5. Save =====
ggsave("phylogeny_with_macroarea.pdf", plot = p, width = 10, height = 12, dpi = 600)
ggsave("phylogeny_with_macroarea.svg", plot = p2, width = 8.27, height = 12)



if (.Platform$OS.type == "windows") shell.exec("phylogeny_with_macroarea.svg")
