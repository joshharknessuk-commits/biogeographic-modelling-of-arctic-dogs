# --- deps ---
library(ape)
library(BioGeoBEARS)

# --- load your run objects (edit paths if needed) ---
load("DECJ_results.Rdata")          # gives resDECj
resDECJ <- resDECj
tree     <- read.tree("FcC_supermatrix.coyote1.rooted.fixed.tre")
tipranges <- getranges_from_LagrangePHYLIP("biogeobears_input_2.geog")

# --- read Biogeotip.csv for ID -> label (and optional Key bolding) ---
tips_df <- read.csv("Biogeotip.csv", stringsAsFactors = FALSE, check.names = FALSE)

# Try to find sensible column names automatically
id_col_candidates  <- c("Taxon","Id","tip","tip_id","tip.label","Tip","Tip_ID","TipID")
lab_col_candidates <- c("Tip_Label","TipLabel","label","Label","new_label","New_Label")
key_col_candidates <- c("Key","IsKeySample","Key_Species","key_species","KeySpecies")

id_col  <- id_col_candidates[id_col_candidates %in% names(tips_df)][1]
lab_col <- lab_col_candidates[lab_col_candidates %in% names(tips_df)][1]
key_col <- key_col_candidates[key_col_candidates %in% names(tips_df)][1]

if (is.na(id_col) || is.na(lab_col)) {
  stop("Biogeotip.csv must contain an ID column (e.g. Taxon/Id/Tip) and a label column (e.g. Tip_Label/Label).")
}

# Build maps
lab_map <- setNames(tips_df[[lab_col]], tips_df[[id_col]])

key_vec <- rep(FALSE, nrow(tips_df))
if (!is.na(key_col)) {
  key_vec <- as.logical(tips_df[[key_col]])
  key_vec[is.na(key_vec)] <- FALSE
}
key_map <- setNames(key_vec, tips_df[[id_col]])

# New labels in the tree’s tip order (+ bolding info)
new_labels <- ifelse(tree$tip.label %in% names(lab_map),
                     lab_map[tree$tip.label],
                     tree$tip.label)

tip_font <- ifelse(tree$tip.label %in% names(key_map) & key_map[tree$tip.label], 2, 1) # 2=bold, 1=plain

# --- “square cladogram” phylo for BioGeoBEARS plotting (keeps numeric lengths) ---
clad2 <- ladderize(tree)

clad2 <- compute.brlen(clad2, 1)  # equal edge lengths -> rectangular look
areas_letters <- getareas_from_tipranges_object(tipranges)         # e.g., c("A","B",...,"L")

# For now, if you already defined area_map elsewhere, keep using it; otherwise show letters:
areas_labels <- areas_letters

# --- make PDF with side legend panel ---
pdf("biogeobears_plot.pdf", width = 15, height = 25)





resDECJ$inputs$plotparams$time_axis_tickmarkpos <- NULL



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



area_map <- c(
  "A"="Alaska","B"="Canada","C"="Central_America","D"="Central_Asia",
  "E"="East_Asia","F"="Europe","G"="Greenland","H"="Middle_East",
  "I"="Russia","J"="Siberia","K"="Southwest_Asia","L"="USA"
)

areas <- colnames(tipranges@df)
k <- resDECJ$inputs$max_range_size
include_null <- isTRUE(resDECJ$inputs$include_null_range)

# make all combos up to k, concatenated like "AB", "ACD", etc.
combos_txt <- unlist(lapply(1:min(k, length(areas)), function(m) {
  apply(combn(areas, m), 2, function(v) paste0(v, collapse = ""))
}), use.names = FALSE)

possible_ranges_list_txt <- c(if (include_null) "_" else character(0), combos_txt)

# --- 1) Map codes to colours via descriptive names ---
area_code_cols <- setNames(bioGeo_colors[area_map], names(area_map))
# e.g., area_code_cols["A"] gives Alaska’s colour

# --- 2) Helpers for combos ---
blend_cols <- function(cols){
  if (length(cols) == 1) return(cols)
  m <- col2rgb(cols)/255; a <- rowMeans(m); rgb(a[1], a[2], a[3])
}
lighten <- function(col, f = 0.20){
  x <- col2rgb(col)/255; x <- pmin(1, x + f*(1-x)); rgb(x[1], x[2], x[3])
}

# --- 3) Assign one colour per state string ---
state_to_color <- function(state_str){
  if (state_str == "_") return("#00000000")    # transparent null
  codes <- strsplit(state_str, "")[[1]]        # split into individual codes
  cols  <- area_code_cols[codes]
  out   <- blend_cols(cols)
  if (length(codes) > 1) out <- lighten(out, 0.20)
  out
}

colors_list_for_states <- vapply(
  possible_ranges_list_txt,
  state_to_color,
  FUN.VALUE = character(1)
)
MLstates <- resDECJ$ML_marginal_prob_each_state_at_branch_top_AT_node
# --- 4) Feed into BioGeoBEARS ---
# Example: MLstates <- resDECJ$ML_marginal_prob_each_state_at_branch_top_AT_node
state_colors <- rangestxt_to_colors(
  possible_ranges_list_txt = possible_ranges_list_txt,
  colors_list_for_states   = colors_list_for_states,
  MLstates                 = MLstates
)






par(mar = c(2,2,2,2))



plot_BioGeoBEARS_results(
  results_object    = resDECJ,
  analysis_titletxt = "BiogeBEARS MODEL: DEC+J (cladogram)",
  addl_params       = list("j"),
  plotwhat          = "pie",
  tr                = clad2,
  statecex          = 0.2,   # small node pies
  splitcex          = 0.2,   # small split pies (or set plotsplits=FALSE to hide)
  tipcex            = 0.0000000000000001,# hide default tip letters
  plotsplits        = TRUE,
  plotlegend        = FALSE,
  titlecex = 0.8,
  colors_list_for_states = colors_list_for_states,
  pie_tip_statecex = 0.00000001,
  xlab = NULL,
  plot_stratum_lines = NULL           # suppresses the timescale grid/axi
)




  # Draw scale on left, root at bottom


mtext(
  side = 3, line = -1.5, cex = 0.5, font = 3,
  text = "Ancestral range estimation for Arctic dog mtDNA phylogeny"
)



# 1) Get the single-area letters actually present in the run
areas_letters <- getareas_from_tipranges_object(tipranges)  # e.g., c("A","B",...,"L")

# 2) Remap to your names
areas_labels  <- unname(area_map[areas_letters])



# Align & draw your CSV-based tip labels (bold where Key==TRUE)
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
xx <- pp$xx; yy <- pp$yy
Ntip <- length(tree$tip.label)

tree_width  <- max(xx) - min(xx)
label_x     <- max(xx) + 0.05 * tree_width   # where to draw the squares (optional)
text_x      <- label_x + 0.01 * tree_width   # where the text goes

# (Optional) dotted connector lines
segments(xx[1:Ntip], yy[1:Ntip], label_x, yy[1:Ntip], lty = "dotted", col = "grey70", lwd = 0.5)

# (Optional) little squares column (comment out if not needed)
rect(label_x, yy[1:Ntip]-0.1, label_x+0.002*tree_width, yy[1:Ntip]+0.1, col = NA, border = NA)

# Tip labels from Biogeotip.csv
text(x = text_x, y = yy[1:Ntip], labels = new_labels, pos = 4, cex = 0.6, font = tip_font)

## Panel 2: external legend (letters or remapped names)
legend("bottomleft",
       inset = c(0, 0.05),   # nudge into bottom margin
       legend = area_map,
       fill   = area_code_cols,
       border = NA, bty = "n",
       cex = 0.9,
       title = "Areas")


dev.off()
if (.Platform$OS.type == "windows") shell.exec("biogeobears_plot.pdf")




