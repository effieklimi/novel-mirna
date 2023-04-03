library(tidyverse) 
library(patchwork) 
library(ComplexHeatmap)
library(RColorBrewer)

miRNAnames <- c("hsa-miR-1827", "hsa-miR-323a-3p", "hsa-miR-449b-5p", "hsa-miR-4774-3p", "hsa-miR-491-3p", "hsa-miR-5681b", "hsa-miR-892b")
vsmc50Top2 <-
  readRDS("results/rds/p01/targets50top2-vsmc-p01.rds") %>%
  lapply( "[", , 3)

names(vsmc50Top2) <- miRNAnames

comb_mat <- make_comb_mat(vsmc50Top2)
my_names <- set_name(comb_mat)

my_set_sizes <- set_size(comb_mat) %>% 
  as.data.frame() %>% 
  rename(sizes = ".") %>% 
  mutate(Set = row.names(.)) 

p1 <- my_set_sizes %>% 
    mutate(Set = reorder(Set, sizes)) %>% 
    ggplot(aes(x = Set, y= sizes)) +
    geom_bar(stat = "identity", aes(fill = Set), alpha = 0.8, width = 0.7) +
    geom_text(aes(label = sizes), size = 5, angle = 90, hjust = 0, y = 1) +
    scale_fill_manual(values = brewer.pal(7, "Set2"), limits = my_names) +
    labs(x = NULL, y = "Set size", fill = NULL) +
    theme_classic() +
    theme(legend.position = "right",
        text = element_text(size= 14),
        axis.ticks.y = element_blank(),
        axis.text = element_blank()
        ) 

get_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

p2 <- get_legend(p1)


my_overlap_sizes <- comb_size(comb_mat) %>% 
  as.data.frame() %>% 
  rename(overlap_sizes = ".") %>% 
  mutate(category = row.names(.))

p3 <- my_overlap_sizes %>% 
  mutate(category = reorder(category, -overlap_sizes)) %>% 
  ggplot(aes(x = category, y = overlap_sizes)) +
  geom_bar(stat = "identity", fill = "grey80", color = NA, alpha = 0.8, width = 0.7) +
  geom_text(aes(label = overlap_sizes, y = 0),
            size = 1, hjust = 0, vjust = 0.5) +
  labs(y = "Intersect sizes",
       x = NULL) +
  theme_classic() +
  theme(text = element_text(size= 5, color = "black"),
        axis.text =element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(hjust = 0),
        ) +
  coord_flip()


my_overlap_matrix <- str_split(string = my_overlap_sizes$category, pattern = "", simplify = T) %>% 
as.data.frame() 

colnames(my_overlap_matrix) <- my_names

my_overlap_matrix_tidy <- my_overlap_matrix %>% 
  cbind(category = my_overlap_sizes$category) %>% 
  pivot_longer(cols = !category, names_to = "Set", values_to = "value") %>% 
  full_join(my_overlap_sizes, by = "category") %>% 
  full_join(my_set_sizes, by = "Set")

p4 <- my_overlap_matrix_tidy %>% 
  mutate(category = reorder(category, -overlap_sizes)) %>%  
  mutate(Set = reorder(Set, sizes)) %>%  
  ggplot(aes(x = Set, y = category)) +
  geom_tile(aes(fill = Set, alpha = value), color = "grey30", size = 0.2) +
  scale_fill_manual(values = brewer.pal(7, "Set2"), # feel free to use other colors 
                    limits = my_names) +
  scale_alpha_manual(values = c(0.8, 0),  # color the grid for 1, don't color for 0. 
                     limits = c("1", "0")) +
  labs(x = "Sets",  
       y = "Overlap") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(color = "black", size= 14),
        panel.grid = element_blank(),
        axis.text = element_blank()
        )



wrap_plots(p1, p2, p4, p3, 
          nrow = 2, 
          ncol = 2,
          heights = c(1, 2), # the more rows in the lower part, the longer it should be
          widths = c(1, 0.8),
          guides = "collect") &
  theme(legend.position = "none")


ggsave("results/figures/upsetplot.pdf", height = 9, width = 4, bg = "white") 
