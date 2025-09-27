# ============================================================
# End-to-end pipeline: time-scaled dogs-only tree with MRCA dots
# Tip labels = accession codes only (no Biogeotip)
# Colour system = bioGeo_colors (your new palette)
# ============================================================

suppressPackageStartupMessages({
  library(ape)
  library(paleotree)
  library(strap)
  library(phangorn)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
})

# ---------------------------
# 0) Colour palette (new colour way)
# ---------------------------
bioGeo_colors <- c(
  # Arctic / Eurasia → Americas gradient
  "Northern_Siberia"      = "#245E90",  # deep navy denim
  "Eastern_Siberia"       = "#6BAED3",  # medium azure
  
  "Alaska"                = "#669f01",  # vivid green-teal
  "North_American_Arctic" = "#3D803D",  # forest green
  "Central_North_America" = "#A6B32D",  # olive/chartreuse
  
  # Arctic adjunct
  "Greenland"             = "#CDEAF9",  # very light icy blue
  
  # Asia
  "East_Asia"             = "#6A56A5",
  "Central_Asia"          = "#B789D0",
  
  # Europe → Middle East → SW Asia
  "Europe"                = "#B8852F",
  "Middle_East"           = "#8F4A2A",
  "Southwest_Asia"        = "#B8684A",
  
  # Americas south
  "Central_America"       = "#B74747"
)

# ============================================================
# 1) Load and prepare tree + fossil ages
# ============================================================
tree <- read.tree("FcC_supermatrix.coyote1.rooted.fixed.tre")
tree_clado <- compute.brlen(tree, 1)
tree_clado_lad <- ladderize(tree_clado)

timeData <- read.csv("paleotreeinput.csv", header = TRUE, row.names = 1)

# Ensure FAD > LAD (older = larger BP). Swap if needed, then convert to Ma.
timeData$FAD <- as.numeric(timeData$FAD)
timeData$LAD <- as.numeric(timeData$LAD)
timeData <- timeData[!is.na(timeData$FAD) & !is.na(timeData$LAD), ]
swap <- which(timeData$FAD < timeData$LAD)
if (length(swap)) timeData[swap, c("FAD", "LAD")] <- timeData[swap, c("LAD", "FAD")]

timeData_Ma <- transform(timeData, FAD = FAD/1e6, LAD = LAD/1e6)

# Drop tips
ids_to_drop <- c(
  "LR745183","CGG16","PA1","Bj1","ThA2","ThA1","ThA16","Fr1","Ms20","ThA5","Ms1","ThA6","Fr2","Ar1",
  "ThA10","HMNH_011","Zh4","Zh5","Bj3","Ms22","Ms2","ThA12","ThA7","ThA15","Ms5","ThA18","Bj6","Ms16",
  "Ms17","Ms10","ThA26","Ms18","ThA33","Ms11","ThA25","Ms12","Ms14","Ms19","ThA24","ThA19","Ms7","Ms8",
  "PRW89","ThA21","ThA38","ThA9","Ms9","ThA37","ThA29","MA7","ThA27","ThA31","ThA32","ThA22","Ms13",
  "ThA36","ThA20","Ms15","CGG20","CGG21","ThA23","ThA28","EU789658","TH14","TU14","TU11","TH10","SK1",
  "WolfHead","Bj5","ThA13","ThA14","MA3","MA6","MA4","MA5","MA9","MA2","MA1","ThA11","ThA8","TU5",
  "TU8","TU839","TH7","TU13","Ms25","Zh1","Zh2","Zh8","Zh6"
)
tree_clado_lad <- drop.tip(tree_clado_lad, intersect(ids_to_drop, tree_clado_lad$tip.label))

cat("Tree tips:", Ntip(tree_clado_lad), "\n")
cat("timeData rows:", nrow(timeData), "\n")
cat("Perfect match?", all(tree_clado_lad$tip.label %in% rownames(timeData)), "\n")

# Root safety
if (!is.rooted(tree)) tree <- midpoint(tree)

# ============================================================
# 2) Wolves vs dogs split + node minima → time-scaling
# ============================================================
# Helpers
safe_mrca <- function(tr, tips){
  tips <- intersect(tips, tr$tip.label)
  if (length(tips) < 2) return(NA_integer_)
  getMRCA(tr, tips)
}
lca_two <- function(tr, n1, n2){
  if (is.na(n1) || is.na(n2)) return(NA_integer_)
  a1 <- c(n1, Ancestors(tr, n1, type="all"))
  a2 <- c(n2, Ancestors(tr, n2, type="all"))
  com <- intersect(a1, a2)
  if (!length(com)) return(NA_integer_)
  is_lowest <- function(n){
    ch <- Descendants(tr, n, type="children")[[1]]
    !any(ch %in% com)
  }
  out <- com[sapply(com, is_lowest)]
  if (length(out)) out[1] else com[1]
}

wolf_ids <- c(
  "LR742739","CGG16","LR745069","AL3284","JK2181","PA1","CGG15","CGG22","CGG26","CGG25",
  "LOW006","VAL_050","AL2744","WolfHead","Bj5","ThA14","AL2657","CGG33","Bj1","ThA2",
  "ThA1","PON012","JK2179","CGG19","VAL_011","Yana1","AL3185","LR745053","TU4","ThA16",
  "Fr1","Ms20","ThA5","Ms1","ThA6","Fr2","Bj3","Ms22","LR742727","LR745186","LR745182",
  "Ms9","ThA37","LR742738","LR745098","LR745192","MA7","ThA22","ThA20","Ms15","Ms13",
  "ThA36","LR742757","ThA29","ThA27","ThA31","ThA32","LR742747","Ms5","ThA18","Ms2",
  "ThA7","ThA15","Ar1","ThA10","HMNH_011","Zh4","Zh5","Bj6","Ms16","Ms10","ThA26","ThA33",
  "Ms18","Ms17","Ms11","ThA25","Ms12","Ms14","Ms19","ThA24","ThA19","Ms7","Ms8","PRW89",
  "ThA21","ThA38","ThA9","CGG20","CGG21","LR742780","ThA23","ThA28","EU789658","TU14",
  "TU11","TH10","SK1","Taimyr1","VAL_005","MA3","MA6","MA4","MA5","MA9","MA2","MA1",
  "CGG23","AL2370","HMNH_007","ThA11","TU5","IRK","TU8","TU839","TH7","TU13","TH12",
  "AL2541","CGG27","LOW002","LOW007","TH8","CGG18","JK2175","LOW003","VAL_037","VAL_012",
  "VAL_033","AL2741","VAL_18A","CGG29","CGG32","CGG12","LOW008","AH574","TH1","TU7","TU6",
  "JK2174","TH3","TU10","Ms25","Zh1","Zh2","Zh8","Zh6","coyote1", "VAL_011"
)
wolf_ids_present <- intersect(wolf_ids, tree_clado_lad$tip.label)
dog_ids_present  <- setdiff(tree_clado_lad$tip.label, wolf_ids_present)

wolf_mrca <- safe_mrca(tree_clado_lad, wolf_ids_present)
dog_mrca  <- safe_mrca(tree_clado_lad, dog_ids_present)
split_node <- lca_two(tree_clado_lad, wolf_mrca, dog_mrca)

# Named clades for node minima (your constraints)
clade1_tips <- c("LR742791","LR742785","LR745078","LR742787","LR745070")  # ≥3 ka
clade2_tips <- c("LR745088","LR742752","LR745174","LR742778","LR745151")  # ≥7.5 ka
clade3_tips <- c("LR742790","LR742809","LR742858","LR742760","LR742793")  # ≥2 ka
pcd_tips    <- c("P59","CICVD","CAO1","AL3223","ISM090","May4","5MT520")  # ≥17 ka
bottom_clade<- c("LR742865","CGG1","LR745129","LR742851","LR742850")      # ≥9.5 ka

node_clade1 <- getMRCA(tree_clado_lad, clade1_tips)
node_clade2 <- getMRCA(tree_clado_lad, clade2_tips)
node_clade3 <- getMRCA(tree_clado_lad, clade3_tips)
node_pcd    <- getMRCA(tree_clado_lad, pcd_tips)
node_bottom <- getMRCA(tree_clado_lad, bottom_clade)

node_mins <- rep(NA_real_, tree_clado_lad$Nnode)
.i <- function(n) n - Ntip(tree_clado_lad)
.set_min <- function(node, val){
  if (is.na(node)) return(invisible(NULL))
  idx <- .i(node)
  if (!is.na(idx) && idx > 0 && idx <= length(node_mins)) node_mins[idx] <<- val
}
.set_min(split_node, 0.027)
.set_min(node_pcd,   0.017)
.set_min(node_clade1,0.003)
.set_min(node_clade2,0.0075)
.set_min(node_clade3,0.002)
.set_min(node_bottom,0.0095)

# Time-scale tree
tstree <- timePaleoPhy(
  tree = tree_clado_lad,
  timeData_Ma,
  type = "equal",
  add.term = TRUE,
  vartime = 2,
  node.mins = node_mins
)

# ============================================================
# 3) Build dogs-only tree and base plot (accession labels only)
# ============================================================
present_wolves <- intersect(wolf_ids, tstree$tip.label)
dogs <- setdiff(tstree$tip.label, present_wolves)
stopifnot(length(dogs) > 1)

dogtree <- keep.tip(tstree, dogs)
dogtree <- ladderize(dogtree, right = TRUE)

# Optional extra pruning
to_drop <- c(
  "LR742754","LR742789","LR745077","LR742755","LR745083","LR745087","LR745178",
  "LR742776","LR745076","LR745068","LR745081","LR745086","LR745080","LR742784","AB499817"
)
drop_now <- intersect(to_drop, dogtree$tip.label)
if (length(drop_now)) {
  dogtree <- drop.tip(dogtree, drop_now)
  dogtree <- ladderize(dogtree, right = TRUE)
  message("Dropped: ", paste(drop_now, collapse = ", "))
}

# Zoom window (e.g., last 25 ka)
pg <- ggtree(dogtree, size = 0.5, layout = "rectangular")
H  <- max(pg$data$x)
x0 <- H - 0.025           # 25 ka window
x1 <- H
label_pad <- 0.0025

# Alternating 2-ka bands
span_ka  <- round((x1 - x0) * 1000)
ka_ticks <- seq(0, span_ka, by = 2)
breaks_x <- x1 - ka_ticks/1000
edges <- sort(unique(breaks_x))
bands <- if (length(edges) >= 2) {
  idx <- seq(1, length(edges) - 1, by = 2)
  data.frame(xmin = edges[idx], xmax = edges[idx + 1], fill = "grey95")
} else NULL

# Base plot
p <- ggplot() +
  geom_rect(data = bands,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            inherit.aes = FALSE, alpha = 0.4) +
  scale_fill_identity() +
  geom_tree(data = dogtree, size = 0.5) +
  theme_tree2() +
  coord_cartesian(xlim = c(x0, x1 + label_pad), expand = FALSE, clip = "on") +
  scale_x_continuous(
    breaks = seq(x0, x1, by = 0.002),
    labels = function(x) round((x1 - x) * 1000)
  ) +
  labs(x = "Time (ka BP)", y = NULL) +
  theme(plot.margin = margin(10, 40, 10, 10),
        panel.grid = element_blank(),
        text = element_text(family = "sans"))

# Tip labels = accession codes only
tree_df <- pg$data
tip_col <- if ("isTip" %in% names(tree_df)) "isTip" else if ("isLeaf" %in% names(tree_df)) "isLeaf" else "leaf"
tip_df  <- tree_df %>% filter(.data[[tip_col]]) %>% select(label, x, y)
pad <- 0.00001
tip_df$x_lab <- tip_df$x + pad

p <- p +
  scale_y_reverse() +
  coord_cartesian(xlim = c(x0, x1 + pad), expand = FALSE, clip = "off") +
  theme(plot.margin = margin(10, 80, 24, 10)) +
  geom_point(data = tip_df, aes(x = x, y = y), size = 0.7, inherit.aes = FALSE) +
  geom_text (data = tip_df, aes(x = x_lab, y = y, label = label),
             size = 2.6, hjust = 0, inherit.aes = FALSE, check_overlap = FALSE) +
  scale_x_continuous(limits = c(x0, x1 + pad), breaks = breaks_x, labels = ka_ticks) +
  theme(axis.text.x  = element_text(margin = margin(t = 2)),
        axis.title.x = element_text(margin = margin(t = 6)))

# ============================================================
# 4) Clade MRCA dots (incl. 10–14, 21–23 + appended 22b, 31–57)
#    Region tagging → bioGeo_colors
# ============================================================
.keep_tips <- function(v) unique(intersect(v, dogtree$tip.label))

# New clades (10–14, 21–23)
clade10_tips <- c("LR745062","LR742867","LR742748","LR742844")
clade11_tips <- c("LR742863","LR745149","LR742807","LR745183")
clade12_tips <- c("LR742832","LR742811","LR742830","LR745137","LR742814")
clade13_tips <- c("LR742855","LR742851","LR742763","LR745118","LR742862")
clade14_tips <- c("LR742855","LR742819","LR742828","LR742815")

clade21_tips <- c("AM310B","OSU611","AL2772","P59","P35")
clade22_tips <- c("5MT501","CINH7","5MT316","ISM256")
clade23_tips <- c("May2","May3","May4","May10")

# Base clades from earlier
clades <- list(
  "Inuit ≥3 ka"    = .keep_tips(clade1_tips),
  "Arctic 7.5 ka"  = .keep_tips(clade2_tips),
  "Inuit ≥2 ka"    = .keep_tips(clade3_tips),
  "PCD ≥17 ka"     = .keep_tips(pcd_tips),
  "Bottom clade"   = .keep_tips(bottom_clade),
  
  "N. American Arctic (10)" = .keep_tips(clade10_tips),
  "N. American Arctic (11)" = .keep_tips(clade11_tips),
  "Greenland (12)"          = .keep_tips(clade12_tips),
  "N. American Arctic (13)" = .keep_tips(clade13_tips),
  "Greenland (14)"          = .keep_tips(clade14_tips),
  
  "Central North America (21)" = .keep_tips(clade21_tips),
  "Central North America (22)" = .keep_tips(clade22_tips),
  "Central America (23)"       = .keep_tips(clade23_tips)
)

# Extra appended clades
clade22b_tips <- c("5MT501","5MT316","ISM256","CAO1")
clade31_tips  <- c("LR742790","LR742809")
clade32_tips  <- c("LR742735","LR742760")
clade33_tips  <- c("LR742743","LR742736","LR74273")
clade34_tips  <- c("LR742746","LR745120","LR742849","LR742848")
clade35_tips  <- c("LR742734","LR742739","PRD1","LR742761")
clade36_tips  <- c("CINHA","LR742736","AL3198")

clade41_tips <- c("KU290873","LR742752","LR742778")
clade42_tips <- c("LR745166","LR742749","LR745151")
clade43_tips <- c("LR742838","LR745191","LR745166")
clade44_tips <- c("LR742749","LR745151")
clade45_tips <- c("LR742841","LR745195","LR742834")
clade46_tips <- c("LR742753","LR742826")
clade47_tips <- c("LR742751","LR742837","LR745179")
clade48_tips <- c("SRR11193499","SRR11193486","SRR11193492")

clade51_tips <- c("LR742785","LR742791")
clade52_tips <- c("LR745078","LR742787")
clade53_tips <- c("LR742745","LR742864")
clade54_tips <- c("LR745069","LR742873","CGG11")
clade55_tips <- c("LR745181","DQ480499","LR742802","LR745187")
clade56_tips <- c("LR742779","LR742750","LR745181")
clade57_tips <- c("LR742779","CGG10","LR745055")

new_clades <- list(
  "Central North America (22b)" = .keep_tips(clade22b_tips),
  
  "Siberia (31)"                = .keep_tips(clade31_tips),
  "North American Arctic (32)"  = .keep_tips(clade32_tips),
  "Siberia (33)"                = .keep_tips(clade33_tips),
  "North American Arctic (34)"  = .keep_tips(clade34_tips),
  "North American Arctic (35)"  = .keep_tips(clade35_tips),
  "Central North America (36)"  = .keep_tips(clade36_tips),
  
  "Siberia (41)"                = .keep_tips(clade41_tips),
  "Siberia (42)"                = .keep_tips(clade42_tips),
  "Greenland (43)"              = .keep_tips(clade43_tips),
  "Siberia (44)"                = .keep_tips(clade44_tips),
  "Greenland (45)"              = .keep_tips(clade45_tips),
  "Siberia (46)"                = .keep_tips(clade46_tips),
  "Greenland (47)"              = .keep_tips(clade47_tips),
  "Greenland (48)"              = .keep_tips(clade48_tips),
  
  "Siberia (51)"                = .keep_tips(clade51_tips),
  "Siberia (52)"                = .keep_tips(clade52_tips),
  "North American Arctic (53)"  = .keep_tips(clade53_tips),
  "North American Arctic (54)"  = .keep_tips(clade54_tips),
  "Greenland (55)"              = .keep_tips(clade55_tips),
  "North American Arctic (56)"  = .keep_tips(clade56_tips),
  "Siberia (57)"                = .keep_tips(clade57_tips)
)

clades <- c(clades, new_clades)

# Compute MRCA per clade (need ≥ 2 tips)
clade_mrca <- lapply(names(clades), function(nm){
  tips <- clades[[nm]]
  if (length(tips) >= 2) data.frame(clade = nm, node = getMRCA(dogtree, tips))
}) |> dplyr::bind_rows()

if (nrow(clade_mrca) > 0) {
  xy <- pg$data |> dplyr::filter(node %in% clade_mrca$node) |> dplyr::select(node, x, y)
  ann_clade <- dplyr::left_join(clade_mrca, xy, by = "node")
  
  # Region tagging → keys in bioGeo_colors
  ann_clade <- ann_clade |>
    dplyr::mutate(region = dplyr::case_when(
      grepl("Greenland", clade, ignore.case = TRUE)             ~ "Greenland",
      grepl("Central North America", clade, ignore.case = TRUE) ~ "Central_North_America",
      grepl("Central America", clade, ignore.case = TRUE)       ~ "Central_America",
      grepl("(N\\.?\\s*American Arctic|North American Arctic)", clade, ignore.case = TRUE)
      ~ "North_American_Arctic",
      grepl("Siberia", clade, ignore.case = TRUE)               ~ "Eastern_Siberia",   # default for “Siberia (...)”
      clade == "PCD ≥17 ka"                                     ~ "Central_North_America", # FIXED typo
      TRUE                                                      ~ "Eastern_Siberia"
    ))
  
  # Plot MRCA dots, coloured by region from bioGeo_colors
  p <- p +
    geom_point(
      data = ann_clade,
      aes(x = x, y = y, colour = region),
      shape = 16, size = 3.8, inherit.aes = FALSE
    ) +
    scale_colour_manual(
      name   = "Ancestral state",
      values = bioGeo_colors,
      guide  = guide_legend(override.aes = list(shape = 16, size = 4, stroke = 0))
    )
}

# ============================================================
# 5) Pairwise MRCA nodes (clade names or raw tip labels)
#    Coloured by inferred region → bioGeo_colors
# ============================================================
pairs_to_mark <- list(
  c("PCD ≥17 ka", "Bottom clade"),
  c("Inuit ≥3 ka","Arctic 7.5 ka"),
  c("Inuit ≥2 ka","Arctic 7.5 ka"),
  c("Greenland (12)","Inuit ≥2 ka"),
  c("N. American Arctic (13)","PCD ≥17 ka"),
  c("Central America (23)","Central North America (22)")
)

resolve_tips <- function(name) {
  if (name %in% names(clades)) return(clades[[name]])
  if (name %in% dogtree$tip.label) return(name)
  message("Warning: '", name, "' not found as clade or tip.")
  character(0)
}

infer_region_for_pair <- function(names_in_pair) {
  # Priority: Greenland > Central_North_America > Central_America > North_American_Arctic > Eastern_Siberia
  if (any(grepl("Greenland", names_in_pair, ignore.case = TRUE))) return("Greenland")
  if (any(grepl("Central North America", names_in_pair, ignore.case = TRUE))) return("Central_North_America")
  if (any(grepl("Central America", names_in_pair, ignore.case = TRUE))) return("Central_America")
  if (any(grepl("(N\\.?\\s*American Arctic|North American Arctic)", names_in_pair, ignore.case = TRUE)) ||
      any(names_in_pair == "PCD ≥17 ka")) return("North_American_Arctic")
  "Eastern_Siberia"
}

pair_rows <- lapply(pairs_to_mark, function(pr){
  tips_union <- unique(unlist(lapply(pr, resolve_tips)))
  if (length(tips_union) < 2) return(NULL)
  nd <- getMRCA(dogtree, tips_union)
  if (is.na(nd)) return(NULL)
  reg <- infer_region_for_pair(pr)
  data.frame(node = nd, region = reg, stringsAsFactors = FALSE)
})

pair_nodes <- dplyr::bind_rows(Filter(Negate(is.null), pair_rows)) |>
  dplyr::distinct(node, .keep_all = TRUE)

if (nrow(pair_nodes)) {
  pair_xy <- pg$data |>
    dplyr::filter(node %in% pair_nodes$node) |>
    dplyr::select(node, x, y) |>
    dplyr::left_join(pair_nodes, by = "node")
  
  p <- p + geom_point(
    data = pair_xy,
    aes(x = x, y = y, colour = region),
    shape = 16, size = 4.2, inherit.aes = FALSE
  )
}
# --- colours for this block (from your palette) ---
russia_col <- unname(bioGeo_colors["Eastern_Siberia"])        # replaces old "#6BAED3"
usa_col    <- unname(bioGeo_colors["Central_North_America"])  # replaces old "#839E2E"

# --- 1) MRCA for each named clade ---
clades <- list(
  "Inuit ≥3 ka"    = clade1_tips,
  "Arctic 7.5 ka"  = clade2_tips,
  "Inuit ≥2 ka"    = clade3_tips,
  "PCD ≥17 ka"     = pcd_tips,
  "Bottom clade"   = bottom_clade
)

clade_mrca <- lapply(names(clades), function(nm){
  tips <- intersect(clades[[nm]], dogtree$tip.label)
  if (length(tips) >= 2) data.frame(clade = nm, node = getMRCA(dogtree, tips))
}) |> dplyr::bind_rows()

xy <- pg$data |> dplyr::filter(node %in% clade_mrca$node) |> dplyr::select(node, x, y)
ann_clade <- dplyr::left_join(clade_mrca, xy, by = "node")

# split PCD vs others (pattern preserved)
ann_pcd    <- ann_clade |> dplyr::filter(clade == "PCD ≥17 ka")
ann_others <- ann_clade |> dplyr::filter(clade != "PCD ≥17 ka")

# draw clade MRCA dots (constant colours; no scale reset)
p <- p +
  geom_point(data = ann_others, aes(x = x, y = y),
             colour = russia_col, shape = 16, size = 3.8, inherit.aes = FALSE) +
  geom_point(data = ann_pcd, aes(x = x, y = y),
             colour = usa_col,    shape = 16, size = 3.8, inherit.aes = FALSE)

# --- 2) Pairwise MRCA nodes (all “Siberia” colour, as in original) ---
pairs_to_mark <- list(
  c("PCD ≥17 ka",  "Bottom clade"),
  c("Inuit ≥3 ka", "Arctic 7.5 ka"),
  c("Inuit ≥2 ka", "Arctic 7.5 ka")
)

# extra pairs you had
if (!exists("CGG3_tips")) {
  CGG3_tips <- grep("^CGG3", dogtree$tip.label, value = TRUE)
}
pairs_to_mark <- c(pairs_to_mark, list(
  c("CGG3",       "Inuit ≥2 ka"),
  c("Inuit ≥2 ka","PCD ≥17 ka")
))

# helper
get_tips <- function(name){
  if (name == "CGG3") return(intersect(CGG3_tips, dogtree$tip.label))
  intersect(clades[[name]], dogtree$tip.label)
}

nodes_pairs <- integer(0)
for (pr in pairs_to_mark) {
  ta <- get_tips(pr[1]); tb <- get_tips(pr[2])
  if (length(unique(c(ta, tb))) >= 2) {
    nodes_pairs <- c(nodes_pairs, getMRCA(dogtree, unique(c(ta, tb))))
  }
}
nodes_pairs <- unique(nodes_pairs[!is.na(nodes_pairs)])

if (length(nodes_pairs)) {
  pair_xy <- pg$data |> dplyr::filter(node %in% nodes_pairs) |> dplyr::select(x, y)
  p <- p + geom_point(
    data = pair_xy, aes(x = x, y = y),
    colour = russia_col, shape = 16, size = 4.2, inherit.aes = FALSE
  )
}
# ============================================================
# 6) Cultural reference lines (optional) + legend + export
# ============================================================
periods <- data.frame(
  name   = c("Thule (ancestors of Inuit)", "Paleo-Eskimo", "First dogs in Americas"),
  age_ka = c(1, 5, 13),
  stringsAsFactors = FALSE
)

# Convert to your plot’s x-scale (Ma from present; older to the left)
# You already computed x1 <- max(pg$data$x) as "present" in Ma.
periods$xpos <- x1 - periods$age_ka/1000

# Get y-range to place labels just below the axis
yr <- range(pg$data$y, na.rm = TRUE)
ymin <- yr[1]
ymax <- yr[2]
y_off <- 0.03 * (ymax - ymin)  # tweak if you need more/less room

# Ensure we can draw outside the panel for bottom labels
p <- p + coord_cartesian(xlim = c(x0, x1 + pad), expand = FALSE, clip = "off") +
  theme(plot.margin = margin(10, 80, 28, 10))  # a bit more bottom margin

# Add thin red dashed lines
p <- p + 
  geom_vline(data = periods, aes(xintercept = xpos),
             colour = "red3", linetype = "dashed", linewidth = 0.3, alpha = 0.9)

# Add small bottom labels with cultural term + age in ka
periods$label <- paste0(periods$name, " (", periods$age_ka, " ka)")
p <- p +
  geom_text(data = periods,
            aes(x = xpos, y = ymax + y_off, label = label),
            angle = 90, vjust = 0, hjust = 0, size = 2.6,
            colour = "red3", inherit.aes = FALSE, fontface = "bold")   # extra bottom space for labels

print(p)

p <- p +
  scale_colour_manual(
    name   = "Ancestral state",
    values = bioGeo_colors,
    labels = c(
      "Eastern_Siberia"       = "Siberia",
      "Northern_Siberia"      = "Northern Siberia",
      "Alaska"                = "Alaska",
      "North_American_Arctic" = "North American Arctic",
      "Central_North_America" = "Central North America",
      "Greenland"             = "Greenland",
      "East_Asia"             = "East Asia",
      "Central_Asia"          = "Central Asia",
      "Europe"                = "Europe",
      "Middle_East"           = "Middle East",
      "Southwest_Asia"        = "Southwest Asia",
      "Central_America"       = "Central America"
    ),
    guide = guide_legend(
      override.aes = list(shape = 16, size = 4, stroke = 0)
    )
  ) +
  theme(
    legend.position      = c(0.02, 0.02),   # bottom-left corner
    legend.justification = c(0, 0),
    legend.background    = element_rect(fill = "white", colour = NA),
    legend.margin        = margin(4, 6, 4, 6),
    legend.key           = element_rect(fill = NA, colour = NA)
  )


outfile <- "dogs_a4.svg"
ggsave(
  filename = outfile,
  plot     = p,
  width    = 8.27,
  height   = 11.69,
  units    = "in",
  device   = "svg"
)

# Auto-open if possible
sys <- Sys.info()[["sysname"]]
if (sys == "Darwin") {
  system2("open", outfile)
} else if (sys == "Windows") {
  shell.exec(normalizePath(outfile))
} else {
  system2("xdg-open", outfile, wait = FALSE)
}







