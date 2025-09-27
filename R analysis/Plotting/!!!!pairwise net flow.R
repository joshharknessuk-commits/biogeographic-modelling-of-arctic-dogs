library(dplyr)
library(tidyr)
library(tibble)
library(gt)

# 1) Read matrix
disp_matrix <- read.table(
  "32e2b0c9-68de-4d33-b8be-3ff7d13eed01.txt",
  header = TRUE, row.names = 1, check.names = FALSE
)

# 2) Arctic subset
arctic_regions <- c("Alaska","Canada","Greenland","Russia","Siberia","USA")

# 3) Long format, Arctic–Arctic only
edges <- disp_matrix[arctic_regions, arctic_regions] %>%
  as.data.frame() %>%
  rownames_to_column("from") %>%
  pivot_longer(-from, names_to = "to", values_to = "weight") %>%
  filter(from != to)

# 4) Keep events >= 1
edges_f <- edges %>% filter(weight >= 1)

# 5) Force Russia to always be Region A
net_pairs <- edges_f %>%
  mutate(
    RegionA = if_else(from == "Russia" | to == "Russia", "Russia", pmin(from, to)),
    RegionB = if_else(RegionA == "Russia", if_else(from == "Russia", to, from), pmax(from, to))
  ) %>%
  group_by(RegionA, RegionB) %>%
  summarise(
    flow_A_to_B = sum(if_else(from == RegionA & to == RegionB, weight, 0)),
    flow_B_to_A = sum(if_else(from == RegionB & to == RegionA, weight, 0)),
    .groups = "drop"
  ) %>%
  mutate(net_flow = round(flow_A_to_B - flow_B_to_A, 2)) %>%
  arrange(desc(net_flow))

# 6) GT table
net_pairs %>%
  select(`Region A` = RegionA, `Region B` = RegionB, `Net Flow` = net_flow) %>%
  gt() %>%
  tab_header(
    title = md("**Arctic & Sub-Arctic Pairwise Net Dispersal Flow (Russia Fixed as Region A)**"),
    subtitle = "Positive = net flow from Russia to Region B; Negative = opposite direction"
  ) %>%
  fmt_number(columns = `Net Flow`, decimals = 2) %>%
  tab_style(
    style = list(cell_fill(color = "#e6f0ff"), cell_text(weight = "bold")),
    locations = cells_body(
      rows = `Region A` == "Russia",
      columns = everything()
    )
  ) %>%
  tab_options(
    table.border.top.width = px(2),
    table.border.bottom.width = px(2),
    heading.align = "left"
  )




