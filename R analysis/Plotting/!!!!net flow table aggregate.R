library(dplyr)
library(tidyr)
library(tibble)
library(gt)

# 1) Read the dispersal matrix
disp_matrix <- read.table(
  "32e2b0c9-68de-4d33-b8be-3ff7d13eed01.txt",
  header = TRUE, row.names = 1, check.names = FALSE
)

# 2) Arctic subset
arctic_regions <- c("Alaska","Canada","Greenland","Russia","Siberia","USA")

# 3) Convert to long format, Arctic–Arctic only
edges_arctic <- disp_matrix[arctic_regions, arctic_regions] %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(-from, names_to = "to", values_to = "weight") %>%
  filter(from != to)

# 4) Keep only weights ≥ 1
edges_filtered <- edges_arctic %>%
  filter(weight >= 1)

# 5) Calculate net flow = exports - imports
net_flow <- edges_filtered %>%
  group_by(from) %>%
  summarise(exports = sum(weight), .groups = "drop") %>%
  full_join(
    edges_filtered %>%
      group_by(to) %>%
      summarise(imports = sum(weight), .groups = "drop"),
    by = c("from" = "to")
  ) %>%
  rename(region = from) %>%
  mutate(
    exports = replace_na(exports, 0),
    imports = replace_na(imports, 0),
    net_flow = round(exports - imports, 2)
  ) %>%
  select(region, net_flow) %>%
  arrange(desc(net_flow))

# 6) GT table
net_flow %>%
  gt() %>%
  tab_header(
    title = md("**Arctic & Sub-Arctic Net Dispersal Flow**"),
    subtitle = "Positive = net exporter, Negative = net importer (≥ 1 events only)"
  ) %>%
  cols_label(
    region   = "Region",
    net_flow = md("**Net Flow**")
  ) %>%
  fmt_number(columns = net_flow, decimals = 2) %>%
  tab_options(
    table.border.top.width = px(2),
    table.border.bottom.width = px(2),
    heading.align = "left"
  )

