#' ednacompare
#'
#' @docType package
#' @name ednacompare
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import reactable
#' @import htmltools
#' @import htmlwidgets
#' @import glue
#' @import rmarkdown
#' @import utils
NULL

theme_colors <- list("one" = "#B7D1DA", "both" = "#95A3A4", "two" = "#E2E8DD", "alert" = "#ffd480")

#' Read dataset
#' @param dataset_path Path to the dataset directory, should contain occurrence.txt and dnaderiveddata.txt
#' @param occurrence_file Occurrence filename
#' @param dna_file DNADerivedData filename
#' @return Data frame
#' @export
read_dataset <- function(dataset_path, occurrence_file = "occurrence.txt", dna_file = "dnaderiveddata.txt") {
  occurrence <- read.table(file.path(dataset_path, occurrence_file), sep = "\t", header = TRUE, na.strings = "", comment.char = "", quote = "")
  dna <- read.table(file.path(dataset_path, dna_file), sep = "\t", header = TRUE, na.strings = "", comment.char = "", quote = "")
  left_join(occurrence, dna, by = "occurrenceID")
}

#' Create common species plot
#' @param ds_1 First dataset
#' @param ds_2 Second dataset
#' @param labels Dataset labels
#' @return A ggplot object
#' @export
create_common_species_plot <- function(ds_1, ds_2, labels) {
  ds_1_species <- ds_1 %>% filter(taxonRank == "species") %>% group_by(target_gene, scientificName) %>% summarize()
  ds_2_species <- ds_2 %>% filter(taxonRank == "species") %>% group_by(target_gene, scientificName) %>% summarize()
  
  df <- bind_rows(
    mutate(ds_1_species, dataset = labels[1]),
    mutate(ds_2_species, dataset = labels[2])
  ) %>% 
    group_by(target_gene, scientificName) %>% 
    summarize(
      dataset = case_when(
        labels[1] %in% dataset & labels[2] %in% dataset ~ "both",
        labels[1] %in% dataset ~ labels[1],
        labels[2] %in% dataset ~ labels[2]
      )
    ) %>% 
    group_by(target_gene, dataset) %>% 
    summarize(n = n()) %>%
    ungroup() %>% 
    complete(target_gene, dataset, fill = list(n = 0))
  
  df_with_both <- filter(df, dataset == "both") %>%
    select(target_gene, both_n = n) %>%
    right_join(df, by = "target_gene")
  
  df_plot <- mutate(df_with_both,
    xmin = case_when(
      dataset == labels[1] ~ -both_n / 2 - n,
      dataset == "both" ~ -both_n / 2,
      dataset == labels[2] ~ both_n / 2
    ),
    xmax = case_when(
      dataset == labels[1] ~ -both_n / 2,
      dataset == "both" ~ both_n / 2,
      dataset == labels[2] ~ both_n / 2 + n
    ),
  )
  
  ggplot(df_plot) +
    geom_vline(xintercept = 0) +
    geom_rect(aes(
      ymin = as.numeric(factor(target_gene)) - 0.4,
      ymax = as.numeric(factor(target_gene)) + 0.4,
      xmin = xmin, xmax = xmax, fill = dataset
    )) +
    scale_y_continuous(
      breaks = seq_along(unique(df_plot$target_gene)),
      labels = unique(df_plot$target_gene)
    ) +
    coord_cartesian(clip = "off") +
    scale_fill_manual(values = setNames(unlist(theme_colors), c(labels[1], "both", labels[2]))) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank()) +
    ylab("Target") +
    xlab("Species")
}

#' Create ASV table
#' @param ds_1 First dataset
#' @param ds_2 Second dataset
#' @param labels Dataset labels
#' @return A reactable table
#' @export
create_asv_table_reactable <- function(ds_1, ds_2, labels) {
  identifications_1 <- ds_1 %>% 
    group_by(target_gene, DNA_sequence, scientificName, taxonRank) %>% 
    summarize() %>% 
    rename(scientificName_1 = scientificName, taxonRank_1 = taxonRank) %>% 
    ungroup()
  
  identifications_2 <- ds_2 %>% 
    group_by(target_gene, DNA_sequence, scientificName, taxonRank) %>% 
    summarize() %>% 
    rename(scientificName_2 = scientificName, taxonRank_2 = taxonRank) %>% 
    ungroup()
  
  identifications <- full_join(identifications_1, identifications_2, by = c("target_gene", "DNA_sequence")) %>% 
    filter(taxonRank_1 == "species" | taxonRank_2 == "species") %>% 
    select(target_gene, DNA_sequence, scientificName_1, scientificName_2, taxonRank_1, taxonRank_2) %>% 
    arrange(target_gene, scientificName_1)
  
  table_data <- mutate(identifications,
    scientificName_1 = ifelse(is.na(scientificName_1), "", scientificName_1),
    scientificName_2 = ifelse(is.na(scientificName_2), "", scientificName_2),
    color_1 = case_when(
      scientificName_1 != "" & scientificName_1 != scientificName_2 & taxonRank_1 == "species" & taxonRank_2 == "species" ~ theme_colors$alert,
      scientificName_1 != "" & scientificName_1 != scientificName_2 & taxonRank_1 == "species" ~ theme_colors$one,
      TRUE ~ ""
    ),
    color_2 = case_when(
      scientificName_2 != "" & scientificName_1 != scientificName_2 & taxonRank_1 == "species" & taxonRank_2 == "species" ~ theme_colors$alert,
      scientificName_2 != "" & scientificName_1 != scientificName_2 & taxonRank_2 == "species" ~ theme_colors$two,
      TRUE ~ ""
    )
  )

  browsable(tagList(
    tags$style(HTML("
    .rt-pagination {
      width: 100%;
      border-top: none !important;
      margin-bottom: 10px;
    }
  ")),
  reactable(
    select(table_data, target_gene, scientificName_1, scientificName_2, DNA_sequence),
    columns = list(
      target_gene = colDef(
        name = "Target",
        width = 80
      ),
      scientificName_1 = colDef(
        name = glue("Scientific name {labels[1]}"),
        style = function(value, index) {
          color <- table_data$color_1[index]
          if (color != "") {
            list(backgroundColor = color)
          } else {
            list()
          }
        }
      ),
      scientificName_2 = colDef(
        name = glue("Scientific name {labels[2]}"),
        style = function(value, index) {
          color <- table_data$color_2[index]
          if (color != "") {
            list(backgroundColor = color)
          } else {
            list()
          }
        }
      ),
      DNA_sequence = colDef(
        name = "Sequence",
        width = 200,
        cell = function(value) {
          div(
            style = "max-width: 200px; overflow-x: auto; white-space: nowrap; font-family: monospace; font-size: small;",
            value
          )
        }
      )
    ),
    pagination = TRUE,
    showPageSizeOptions = TRUE,
    pageSizeOptions = c(10, 20, 100, 1000),
    defaultPageSize = 10,
    searchable = FALSE,
    sortable = TRUE,
    filterable = TRUE
  ) %>%
    onRender("
    function(el, x) {
      const pagination = el.querySelector('.rt-pagination');
      if (pagination && pagination.parentNode) {
        el.insertBefore(pagination, el.firstChild);
      }
    }
  ")
  ))
}

#' Compare datasets
#' @param ds_1 First dataset
#' @param ds_2 Second dataset
#' @param labels Dataset labels
#' @return A list
#' @export
compare_datasets <- function(ds_1, ds_2, labels = c("1", "2")) {
  common_species_plot <- create_common_species_plot(ds_1, ds_2, labels)
  asv_table <- create_asv_table_reactable(ds_1, ds_2, labels)
  list(
    common_species_plot = common_species_plot,
    asv_table = asv_table
  )
}

#' Generate comparison report
#' @param results Results list
#' @param output_file Output filename
#' @return None
#' @export
generate_comparison_report <- function(results, output_file = "output.html") {
  template_path <- system.file("template.Rmd", package = "ednacompare")
  output_dir <- dirname(output_file)
  output_name <- basename(output_file)
  
  rmarkdown::render(
    input = template_path,
    output_file = output_name,
    output_dir = output_dir,
    params = results,
    envir = new.env(parent = globalenv())
  )
  
  message("Report written to: ", normalizePath(file.path(output_dir, output_name)))
}
