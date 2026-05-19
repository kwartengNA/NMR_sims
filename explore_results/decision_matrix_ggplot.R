# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(forcats)

# 1. DATA PREPARATION
cols <- c("NMA", "CCIE", "ICIE", "IUMIE", "CUMIE")
rows <- c(
  "Model: NMA*", "Model: CCIE*", "Model: ICIE*", "Model: IUMIE*", 
  "Model: CUMIE*", "Sparse Network†", "Multi-Arm Network‡", "High Heterogeneity*"
)

data_matrix <- matrix(c(
  "✔️", "✔️", "✔️", "✔️", "✔️",       # Row 1
  "❌", "✔️", "✔️", "✔️", "❌",       # Row 2,  
  "❌", "❌", "✔️", "✔️", "❌",       # Row 3
  "❌", "❌", "❌", "✔️", "❌",       # Row 4
  "❌", "❌", "❌", "✔️", "✔️",       # Row 5
  "✔️", "✔️", "❌", "❌", "✔️",       # Row 6
  "✔️", "✔️", "✔️", "❌", "❌",       # Row 7
  "⚠", "⚠", "⚠", "⚠", "⚠"        # Row 8
), nrow = 8, byrow = TRUE)

data_matrix_minimalist <- matrix(c(
  "✓", "✓", "✓", "✓", "✓",       # Row 1
  "✕", "✓", "✓", "✓", "✕",       # Row 2
  "✕", "✕", "✓", "✓", "✕",       # Row 3
  "✕", "✕", "✕", "✓", "✕",       # Row 4
  "✕", "✕", "✕", "✓", "✓",       # Row 5
  "✓", "✓", "✕", "✕", "✓",       # Row 6
  "✓", "✓", "✓", "✕", "✕",       # Row 7
  "⚠", "⚠", "⚠", "⚠", "⚠"        # Row 8
), nrow = 8, byrow = TRUE)


# 2.  COLOR PALETTE ############
# Fill colors (Soft/Pastel traffic light)
fill_colors <- c(
  "Recommended" = "#EBFFEE",     # Soft Sky Blue
  "Not Recommended" = "#FFEBEB",  # Ultra Light Gray
  "Use with caution" = "#FFFBEB"   # Soft Cream/Amber
)
corp_colors <- c(
  "Recommended" = "#2C5282",      # Muted Navy
  "Not Recommended" = "#A0AEC0", # Muted Slate Gray
  "Use with caution" = "#B7791F"  # Muted Amber
)

# Text colors (Deep/Contrasted)
text_colors <- c(
  "Recommended" = "darkgreen",     # Deep Navy
  "Not Recommended" = "darkred", # Slate Gray
  "Use with caution" = "#B7791F"  # Deep Amber
)


# choose which input data/style ########
df <- as.data.frame(data_matrix_minimalist)# or data_matrix_minimalist
colnames(df) <- cols
df$RowLabel <- rows


df_long <- df %>%
  pivot_longer(cols = -RowLabel, names_to = "ColumnLabel", values_to = "Symbol") %>%
  mutate(
    Status = case_when(
      Symbol == "✔️" ~ "Recommended",
      Symbol == "❌" ~ "Not Recommended",
      Symbol == "⚠" ~ "Use with caution",
      Symbol == "✓" ~ "Recommended", 
      Symbol == "✕" ~ "Not Recommended"
    ),
    RowLabel = (factor(RowLabel, levels = rev(rows))),
    ColumnLabel = factor(ColumnLabel, levels = cols)
  )

# 3. GGPLOT SYNTAX traffic light #########
ggplot(df_long, aes(x = ColumnLabel, y = RowLabel)) +
  # Filled Cells with light gray borders
  geom_tile(aes(fill = Status), color = "#E2E8F0", size = 0.7) +
  # Symbols with corresponding deep colors
  geom_text(aes(label = Symbol, color = Status), 
            size = 6, fontface = "bold", family = "sans") +
  # Assign the manual color scales
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = text_colors) +
  # Layout & Labels
  labs(
    x = "Implemented Model", 
    y = "Data Generating Attribute", fill = NULL, color = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    
    axis.title = element_text(
      face = "bold", 
      size = 11, 
      hjust = 0.5, 
      margin = margin(b = 10, t=10),
      color = "black"),
    
    # Axis Styling
    axis.text.x = element_text(
      face = "bold", 
      size = 11, 
      color = "black", 
      margin = margin(t = 10)),
    axis.text.y = element_text(
      face = "bold", size = 11, 
      color = "black", hjust = 1, 
      margin = margin(r = 15)),
    
    # Legend Styling
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_blank(),
    
    # Spacing
    plot.margin = margin(30, 30, 30, 30)
  ) +
  scale_x_discrete(position = "top")

# 3. GGPLOT SYNTAX # minimalist ########################
ggplot(df_long, aes(x = ColumnLabel, y = RowLabel)) +
  # Create the cell background (white tiles with light gray borders)
  geom_tile(fill = "white", color = "#E2E8F0", size = 0.5) +
  # Add the symbols
  geom_text(aes(label = Symbol, color = Status), 
            size = 6, fontface = "bold", family = "sans") +
  # Custom colors
  scale_color_manual(values = corp_colors) +
  # Labels and Titles
  labs(
    title = "Implemented Model Recommendations",
    x = NULL, y = NULL, color = NULL
  ) +
  # Refined Theme
  theme_minimal() +
  theme(
    # Background and Grid
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    
    # Title styling
    plot.title = element_text(
      face = "bold", 
      size = 16, 
      hjust = 0, 
      margin = margin(b = 20),
      color = "#2D3748"
    ),
    
    # Axis styling
    axis.text.x = element_text(
      face = "bold", 
      size = 11, 
      color = "#4A5568", 
      margin = margin(t = 10)
    ),
    axis.text.y = element_text(
      size = 11, 
      color = "#4A5568", 
      hjust = 1, 
      margin = margin(r = 15)
    ),
    
    # Legend styling
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10, color = "#4A5568"),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    
    # Spacing
    plot.margin = margin(30, 30, 30, 30)
  ) +
  # Move X axis to top to behave like a table header
  scale_x_discrete(position = "top")
# Save as Landscape
# ggsave("corporate_matrix_filled.png", width = 11, height = 6, dpi = 300)