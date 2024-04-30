
#' Volcano
#'
#' @param dataset
#' @param type (p-value adjustment. Either p, xiao, or q)
#' @param threshold (cut-off)
#'
#' @return a volcano plot

volcano <- function (dataset, type, threshold) {

    dataset%>%
        mutate(color = case_when(
            logFC >= 0 & {{type}} <= {{threshold}} ~ "Upregulated",
            logFC <= 0 & {{type}} <= {{threshold}} ~ "Downregulated",
            TRUE ~ "Unchanged")) %>%
        ggplot(aes(x=logFC, y=-log10(P.Value), label=row.names(dataset)))+
        geom_point(aes(color = color, alpha=color), size = 3)+
        theme(panel.background = element_rect(fill = "transparent", colour = NA),
              plot.background = element_rect(fill = "transparent", colour = NA),
              panel.grid.minor=element_blank(),
              panel.grid.major = element_blank(),
              axis.line = element_line(colour = "black"),
              legend.text=element_text(size=10),
              text = element_text(size = 12),
              legend.title = element_blank(),
              legend.key = element_blank(),
              plot.title = element_text(hjust = 0.5, face="bold"),
              axis.text.x = element_text(color="black", size=10),
              axis.text.y = element_text(color="black", size=10),
              legend.position = "none")+
        geom_text_repel(point.size=4, size=3, min.segment.length = 0.1, force=0.3)+
        scale_color_manual(breaks = c("Upregulated", "Downregulated", "Unchanged"),
                           values=c("dodgerblue3", "firebrick3", "gray50"))+
        scale_alpha_manual(breaks = c("Upregulated", "Downregulated", "Unchanged"),
                           values=c(1, 1, 0.1))+
        xlab("Log2fold change (post-pre)") + ylab("-log10(p)")+
        xlim(-2.5, 2.5)+
        ylim(0,9)
}


#' Volcano for interaction plots
#'
#' @param dataset
#' @param type (p-value adjustment. Either p, xiao, or q)
#' @param threshold (cut-off)
#'
#' @return a volcano plot

volcano_interaction <- function (dataset, type, threshold) {

    dataset%>%
        mutate(color = case_when(
            logFC >= 0 & {{type}} <= {{threshold}} ~ "Upregulated",
            logFC <= 0 & {{type}} <= {{threshold}} ~ "Downregulated",
            TRUE ~ "Unchanged")) %>%
        ggplot(aes(x=logFC, y=-log10(P.Value), label=row.names(dataset)))+
        geom_point(aes(color = color, alpha=color), size = 3)+
        theme(panel.background = element_rect(fill = "transparent", colour = NA),
              plot.background = element_rect(fill = "transparent", colour = NA),
              panel.grid.minor=element_blank(),
              panel.grid.major = element_blank(),
              axis.line = element_line(colour = "black"),
              legend.text=element_text(size=10),
              text = element_text(size = 12),
              legend.title = element_blank(),
              legend.key = element_blank(),
              plot.title = element_text(hjust = 0.5, face="bold"),
              axis.text.x = element_text(color="black", size=10),
              axis.text.y = element_text(color="black", size=10),
              legend.position = "none")+
        geom_text_repel(point.size=4, size=3, min.segment.length = 0.1, force=0.3)+
        scale_color_manual(breaks = c("Upregulated", "Downregulated", "Unchanged"),
                           values=c("dodgerblue3", "firebrick3", "gray50"))+
        scale_alpha_manual(breaks = c("Upregulated", "Downregulated", "Unchanged"),
                           values=c(1, 1, 0.1))+
        xlab("Log2 difference (type II - type I)") + ylab("-log10(p)")+
        xlim(-2.5, 2.5)+
        ylim(0,9)
}



#' Volcano with continuous color based on type
#'
#' @param dataset
#' @param type
#' @param threshold
#'
#' @return a volcano plot

volcano_con <- function (dataset, type) {

    xiao_gradient <- colorRampPalette(c((colorRampPalette(c("gray","gray", "gray", "gray","gray", "#F94040"))(50)), rev(colorRampPalette(c("gray","gray","gray","gray","gray", "#5757F9"))(50))))

    dataset%>%
        dplyr::mutate(color = ifelse(logFC>0, {{type}}, {{type}}*-1)) %>%
        ggplot(aes(x=logFC, y=-log10(P.Value), label=row.names(dataset)))+
        geom_point(aes(color = color), size = 3)+
        theme(panel.background = element_rect(fill = "transparent", colour = NA),
              plot.background = element_rect(fill = "transparent", colour = NA),
              panel.grid.minor=element_blank(),
              panel.grid.major = element_blank(),
              axis.line = element_line(colour = "black"),
              legend.text=element_text(size=10),
              text = element_text(size = 12),
              legend.title = element_blank(),
              legend.key = element_blank(),
              plot.title = element_text(hjust = 0.5, face="bold"),
              axis.text.x = element_text(color="black", size=10),
              axis.text.y = element_text(color="black", size=10))+
        geom_text_repel(point.size=4, size=3, min.segment.length = 0.1, force=0.3)+
        xlab("Log2fold change") + ylab("-log10(p)")+
        scale_colour_gradientn(colors=xiao_gradient(100))
}


#' Title
#'
#' @param se
#'
#' @return Returns a heatmap visualizing missing values (in white) and valid values (in black) and column annotations from metadata in the form of intervention, time, id, and fiber type.
#' @export Nothing
#'
#' @examples missing_plot(se_ter_i)
missing_plot <- function(se){

#Create a new data frame containing only proteins with missing values
df_missing <- SummarizedExperiment::assay({{se}}) %>%
    dplyr::filter(!complete.cases(.)) %>% #Removes rows with no missing values
    dplyr::mutate_all(~ ifelse(is.na(.), 0, 1)) #NA's replaced with 0 and everything else with 1

#Heatmap
heatmap_missing <- pheatmap(df_missing,
                            cluster_rows = F,
                            cluster_cols = T,
                            annotation_col = dplyr::select(metadata, c("time", "intervention", "id")),
                            annotation_colors=list(time=c(pre="#ff6361", post="#ffa600"),
                                                   intervention=c(fb="#c06c85", lb="#345c7e")),
                            show_rownames = F,
                            color=colorRampPalette(c("white", "black"))(2),
                            #cellwidth =5,
                            border_color = NA,
                            legend = F,
                            annotation_legend = T
)
}


#' Title
#'
#' @param summarized_df e.g. df_long_l2fc_mean (string)
#' @param intervention e.g. terbutaline or resistance (string)
#' @param fiber_type e.g. I or II (string)
#'
#' @return a bar plot of mitochondrial subunits ordered by median l2fc.
#' @export
#'
#' @examples
mito_subunits_median <- function(summarized_df, intervention, fiber_type) {

df_long_l2fc_mean %>%
    dplyr::filter(intervention == {{intervention}} & fiber_type == {{fiber_type}}) %>%
    dplyr::filter(grepl("subunit", mito, ignore.case = T)) %>%
    dplyr::filter(!is.na(l2fc_median)) %>%
    ggplot2::ggplot(aes(x=reorder(protein, l2fc_median), y=l2fc_median, fill=fiber_type))+
    geom_col()+
    theme(
        axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1, color="black", size=10),
        axis.text.y= element_text(color="black", size=10),
        axis.line = element_line(color="black", size = 0.5),
        axis.ticks = element_blank(),
        panel.background = element_rect(color = "black", fill=NA, size = 0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(size = 12),
        text = element_text(family = "Source Sans Pro", size=11),
        legend.position = "none"
    )+
    scale_fill_manual(values = c("I" = "#cb4c52",
                                 "II" = "#7c98ce"))+
    labs(y = "Median l2fc", x = "")+
    ylim(-.8, 0.4)
}



#' Title
#'
#' @param summarized_df
#'
#' @return A violion plot of median l2fc of mitochondrial subunits stacked.
#' @export
#'
#' @examples
mito_subunits_summed_median <- function (summarized_df) {

    df_long_l2fc_mean %>%
    dplyr::filter(grepl("subunit", mito, ignore.case = T)) %>%
    ggplot2::ggplot(aes(x=intervention, y=l2fc_median, fill=fiber_type))+
    geom_violin(aes(group=interaction(fiber_type, intervention)), position =position_dodge(width=1))+
    geom_point(aes(fill=fiber_type), na.rm=TRUE, position = position_dodge(width = 1), size=5, alpha=0.2)+
    geom_hline(yintercept = 0, linetype = "dashed")+
    theme(
        axis.text.x= element_text(color="black", size = 12),
        axis.text.y= element_text(color="black", size = 12),
        axis.line = element_line(color="black", size = 0.5),
        panel.background = element_rect(color = "black", fill=NA, size = 0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(size = 12),
        plot.margin = unit(c(5, 5, 5, 25), "pt"),
        text = element_text(family = "Source Sans Pro", size=12)
    )+
    scale_fill_manual(values=c("#cb4c52", "#7c98ce"),
                      name = "Fiber Type",
                      labels = c("I",
                                 "II"))+
    ggtitle("Mitochondrial subunits")+
    labs(x = "", y = "Summed median l2fc")+
    scale_x_discrete(labels = c(resistance = "Resistance training", terbutaline = expression("Beta"[2]*"-agonist")))
}

#' Title
#'
#' @param summarized_df e.g. df_long_l2fc_mean (string)
#' @param intervention e.g. terbutaline or resistance (string)
#' @param fiber_type e.g. I or II (string)
#'
#' @return a bar plot of mitochondrial subunits ordered by median l2fc.
#' @export
#'
#' @examples
ribo_subunits_median <- function(summarized_df, intervention, fiber_type) {

    df_long_l2fc_mean %>%
        dplyr::filter(intervention == {{intervention}} & fiber_type == {{fiber_type}}) %>%
        dplyr::filter(grepl('cytoplasmic translation', gobp, ignore.case = T)) %>%
        dplyr::filter(!is.na(l2fc_median)) %>%
        ggplot2::ggplot(aes(x=reorder(protein, l2fc_median), y=l2fc_median, fill=fiber_type))+
        geom_col()+
        theme(
            axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1, color="black", size=10),
            axis.text.y= element_text(color="black", size=10),
            axis.line = element_line(color="black", size = 0.5),
            axis.ticks = element_blank(),
            panel.background = element_rect(color = "black", fill=NA, size = 0.5),
            panel.grid.minor=element_blank(),
            panel.grid.major = element_blank(),
            plot.background = element_blank(),
            plot.title = element_text(size = 12),
            text = element_text(family = "Source Sans Pro", size=11),
            legend.position = "none"
        )+
        scale_fill_manual(values = c("I" = "#cb4c52",
                                     "II" = "#7c98ce"))+
        labs(y = "Median l2fc", x = "")+
        ylim(-.5, 0.65)
}

#' Title
#'
#' @param summarized_df
#'
#' @return A violion plot of median l2fc of mitochondrial subunits stacked.
#' @export
#'
#' @examples
ribo_subunits_summed_median <- function (summarized_df) {

    df_long_l2fc_mean %>%
        dplyr::filter(grepl('cytosolic ribosome', gocc, ignore.case = T)) %>%
        ggplot2::ggplot(aes(x=intervention, y=l2fc_median, fill=fiber_type))+
        geom_violin(aes(group=interaction(fiber_type, intervention)), position =position_dodge(width=1))+
        geom_point(aes(fill=fiber_type), na.rm=TRUE, position = position_dodge(width = 1), size=5, alpha=0.2)+
        geom_hline(yintercept = 0, linetype = "dashed")+
        theme(
            axis.text.x= element_text(color="black", size = 12),
            axis.text.y= element_text(color="black", size = 12),
            axis.line = element_line(color="black", size = 0.5),
            panel.background = element_rect(color = "black", fill=NA, size = 0.5),
            panel.grid.minor=element_blank(),
            panel.grid.major = element_blank(),
            plot.background = element_blank(),
            plot.title = element_text(size = 12),
            plot.margin = unit(c(5, 5, 5, 25), "pt"),
            text = element_text(family = "Source Sans Pro", size=12)
        )+
        scale_fill_manual(values=c("#cb4c52", "#7c98ce"),
                          name = "Fiber Type",
                          labels = c("I",
                                     "II"))+
        ggtitle("Ribosomal proteins")+
        labs(x = "", y = "Summed median l2fc")+
        scale_x_discrete(labels = c(resistance = "Resistance training", terbutaline = expression("Beta"[2]*"-agonist")))
}

#' Title
#'
#' @param summarized_df e.g. df_long_l2fc_mean (string)
#' @param intervention e.g. terbutaline or resistance (string)
#' @param fiber_type e.g. I or II (string)
#' @param term string of search word
#' @param term_type string of ontology. Keywords, gobp, gocc, gomf, or mito.
#'
#' @return a bar plot of any term ordered by median l2fc.
#' @export
#'
#' @examples
term_median <- function(summarized_df, intervention, fiber_type, term, term_type) {

    df_long_l2fc_mean %>%
        dplyr::filter(intervention == {{intervention}} & fiber_type == {{fiber_type}}) %>%
        dplyr::filter(grepl({{term}}, {{term_type}}, ignore.case = T)) %>%
        dplyr::filter(!is.na(l2fc_median)) %>%
        ggplot2::ggplot(aes(x=reorder(protein, l2fc_median), y=l2fc_median, fill=fiber_type))+
        geom_col()+
        theme(
            axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1, color="black", size=10),
            axis.text.y= element_text(color="black", size=10),
            axis.line = element_line(color="black", size = 0.5),
            axis.ticks = element_blank(),
            panel.background = element_rect(color = "black", fill=NA, size = 0.5),
            panel.grid.minor=element_blank(),
            panel.grid.major = element_blank(),
            plot.background = element_blank(),
            plot.title = element_text(size = 12),
            text = element_text(family = "Source Sans Pro", size=11),
            legend.position = "none"
        )+
        scale_fill_manual(values = c("I" = "#cb4c52",
                                     "II" = "#7c98ce"))+
        labs(y = "Median l2fc", x = "")
}

#' Title
#'
#' @param summarized_df
#'
#' @return A violion plot of median l2fc of mitochondrial subunits stacked.
#' @param term string of search word
#' @param term_type string of ontology. Keywords, gobp, gocc, gomf, or mito.
#' @export
#'
#' @examples
term_summed_median <- function (summarized_df, term, term_type) {

    df_long_l2fc_mean %>%
        dplyr::filter(grepl({{term}}, {{term_type}}, ignore.case = T)) %>%
        ggplot2::ggplot(aes(x=intervention, y=l2fc_median, fill=fiber_type))+
        geom_violin(aes(group=interaction(fiber_type, intervention)), position =position_dodge(width=1))+
        geom_point(aes(fill=fiber_type), na.rm=TRUE, position = position_dodge(width = 1), size=5, alpha=0.2)+
        geom_hline(yintercept = 0, linetype = "dashed")+
        theme(
            axis.text.x= element_text(color="black", size = 12),
            axis.text.y= element_text(color="black", size = 12),
            axis.line = element_line(color="black", size = 0.5),
            panel.background = element_rect(color = "black", fill=NA, size = 0.5),
            panel.grid.minor=element_blank(),
            panel.grid.major = element_blank(),
            plot.background = element_blank(),
            plot.title = element_text(size = 12),
            plot.margin = unit(c(5, 5, 5, 25), "pt"),
            text = element_text(family = "Source Sans Pro", size=12)
        )+
        scale_fill_manual(values=c("#cb4c52", "#7c98ce"),
                          name = "Fiber Type",
                          labels = c("I",
                                     "II"))+
        labs(x = "", y = "Summed median l2fc")+
        scale_x_discrete(labels = c(resistance = "Resistance training", terbutaline = expression("Beta"[2]*"-agonist")))
}

#' Title
#'
#' @param term
#' @param term_type
#'
#' @return
#' @export
#'
#' @examples
term_overview <- function(term, term_type) {

    ter_i <- term_median(df_long_l2fc_mean, "terbutaline", "I", {{term}}, {{term_type}}) + ggtitle(expression(paste("Beta"[2], "-agonist, type I")), )+theme(plot.title = element_text(vjust = -1))
    ter_ii <- term_median(df_long_l2fc_mean, "terbutaline", "II", {{term}}, {{term_type}}) + ggtitle(expression(paste("Beta"[2], "-agonist, type IIa")), )+theme(plot.title = element_text(vjust = -1))
    res_i <- term_median(df_long_l2fc_mean, "resistance", "I", {{term}}, {{term_type}})+ ggtitle("Resistance training, type I")
    res_ii <- term_median(df_long_l2fc_mean, "resistance", "II", {{term}}, {{term_type}}) + ggtitle("Resistance training, type IIa")

    combined <- term_summed_median(df_long_l2fc, {{term}}, {{term_type}})

    layout <- grid.arrange(ter_i, ter_ii, res_i, res_ii, combined,
                              widths = c(2, 2, 1),
                              layout_matrix = rbind(c(1, 2, 5),
                                                    c(3, 4, 5)))

   ggplot2::ggsave(here::here(paste0("data/figures/protein_groups/", term, ".svg")), plot=layout, width = 25, height = 5)

}



#' GSEA bubble plot for resistance exercise - downregulated
#'
#' @param gsea_object
#'
#' @return Plot
#' @export
#'
#' @examples
gsea_res_down <- function(gsea_object) {

    #Make GSEA object dataframe and add desired labels
    df_gsea_res <- as.data.frame(gsea_res_both_bp) %>%
        dplyr::rename(fiber_type = Cluster,
                      description = Description) %>%
        dplyr::mutate(label = case_when(
            description == "cellular respiration" ~ description,
            description == "cell migration"~ description,
            description == "gluconeogenesis"~ description,
            description == "muscle organ development"~ description,
            description == "muscle contraction"~ description,
            description == "cell adhesion"~ description,
            description == "muscle structure development"~ description,
            description == "cytoskeleton organization"~ description,
            description == "vesicle−mediated transport"~ description,
            description == "positive regulation of growth"~ description)) %>%
        dplyr::mutate(log10p = -log10(pvalue), #Add -log10p column
                      count = str_count(.$core_enrichment, '/')+1) %>%
        dplyr::filter(NES < 0)

    #Plot
    ggplot2::ggplot(df_gsea_res, aes(x = log10p, y = NES, label = label, size = count))+
        geom_point(aes(color = fiber_type))+
        geom_text_repel(point.size=6, size=4, min.segment.length = 0, max.overlaps = Inf, point.padding = 0.2, box.padding = 0.5)+
        theme(
            axis.text.x= element_text(color="black"),
            axis.text.y= element_text(color="black"),
            axis.title.y = element_text(),
            axis.line = element_line(color="black", linewidth = 0.2),
            panel.background = element_rect(color = "black", fill=NA, size = 0.2),
            panel.grid.minor=element_blank(),
            panel.grid.major = element_blank(),
            plot.background = element_blank(),
            #plot.title = element_text(),
            plot.margin = unit(c(5, 5, 5, 25), "pt"),
            text = element_text(family = "Source Sans Pro", size = 15),
            #legend.position = c(0, 1),
            #legend.justification = c(0, 1), #Place legend in top left corner
            legend.background = element_blank(),
            legend.text = element_text(),
            legend.title = element_text()
        )+
        scale_color_manual(values=c("#c54e56", "#8398cd"),
                           name = "Fiber type",
                           labels = c("I" = "Type I",
                                      "II" = "Type IIa"))+
        labs(y = "Normalized enrichment score", x="-Log10(p-value)")
}

#' GSEA bubble plot for resistance exercise - upregulated
#'
#' @param gsea_object
#'
#' @return Plot
#' @export
#'
#' @examples
gsea_res_up <- function(gsea_object) {

    #Make GSEA object dataframe and add desired labels
    df_gsea_res <- as.data.frame(gsea_res_both_bp) %>%
        dplyr::rename(fiber_type = Cluster,
                      description = Description) %>%
        dplyr::mutate(label = case_when(
            description == "cellular respiration" ~ description,
            description == "cell migration"~ description,
            description == "gluconeogenesis"~ description,
            description == "muscle organ development"~ description,
            description == "muscle contraction"~ description,
            description == "cell adhesion"~ description,
            description == "muscle structure development"~ description,
            description == "cytoskeleton organization"~ description,
            description == "vesicle−mediated transport"~ description,
            description == "positive regulation of growth"~ description)) %>%
        dplyr::mutate(log10p = -log10(pvalue), #Add -log10p column
                      count = str_count(.$core_enrichment, '/')+1) %>%
        dplyr::filter(NES > 0)

    #Plot
    ggplot2::ggplot(df_gsea_res, aes(x = log10p, y = NES, label = label, size = count))+
        geom_point(aes(color = fiber_type))+
        geom_text_repel(point.size=6, size=4, min.segment.length = 0, max.overlaps = Inf, point.padding = 0.2, box.padding = 0.5)+
        theme(
            axis.text.x= element_text(color="black"),
            axis.text.y= element_text(color="black"),
            axis.title.y = element_text(),
            axis.line = element_line(color="black", linewidth = 0.2),
            panel.background = element_rect(color = "black", fill=NA, size = 0.2),
            panel.grid.minor=element_blank(),
            panel.grid.major = element_blank(),
            plot.background = element_blank(),
            #plot.title = element_text(),
            plot.margin = unit(c(5, 5, 5, 25), "pt"),
            text = element_text(family = "Source Sans Pro", size = 15),
            #legend.position = c(0, 1),
            #legend.justification = c(0, 1), #Place legend in top left corner
            legend.background = element_blank(),
            legend.text = element_text(),
            legend.title = element_text()
        )+
        scale_color_manual(values=c("#c54e56", "#8398cd"),
                           name = "Fiber type",
                           labels = c("I" = "Type I",
                                      "II" = "Type IIa"))+
        labs(y = "Normalized enrichment score", x="-Log10(p-value)")
}


#' GSEA bubble plot for resistance exercise - upregulated
#'
#' @param gsea_object
#'
#' @return Plot
#' @export
#'
#' @examples
gsea_ter_up <- function(gsea_object) {

    #Make GSEA object dataframe and add desired labels
    df_gsea_ter <- as.data.frame(gsea_ter_both_bp) %>%
        dplyr::rename(fiber_type = Cluster,
                      description = Description) %>%
        dplyr::mutate(label = case_when(
            description == "cellular respiration" ~description,
            description == "mitochondrial respiratory chain complex assembly" ~description,
            description == "cytoskeleton organization" ~description,
            description == "regulation of cell morphogenesis" ~description,
            description == "mitotic cell cycle" ~description,
            description == "proton transmembrane transport" ~description,
            description == "cytoplasmic translation" ~description,
            description == "vesicle−mediated transport" ~description,
            description == "cellular response to lipid" ~description,
            description == "negative regulation of ion transmembrane transporter activity" ~description)) %>%
        dplyr::mutate(log10p = -log10(pvalue), #Add -log10p column
                      count = str_count(.$core_enrichment, '/')+1) %>%
        dplyr::filter(NES > 0)

    #Plot
    ggplot2::ggplot(df_gsea_ter, aes(x = log10p, y = NES, label = label, size = count))+
        geom_point(aes(color = fiber_type))+
        geom_text_repel(point.size=6, size=4, min.segment.length = 0, max.overlaps = Inf, point.padding = 0.2, box.padding = 0.5)+
        theme(
            axis.text.x= element_text(color="black"),
            axis.text.y= element_text(color="black"),
            axis.title.y = element_text(),
            axis.line = element_line(color="black", linewidth = 0.2),
            panel.background = element_rect(color = "black", fill=NA, size = 0.2),
            panel.grid.minor=element_blank(),
            panel.grid.major = element_blank(),
            plot.background = element_blank(),
            #plot.title = element_text(),
            plot.margin = unit(c(5, 5, 5, 25), "pt"),
            text = element_text(family = "Source Sans Pro", size = 15),
            #legend.position = c(0, 1),
            #legend.justification = c(0, 1), #Place legend in top left corner
            legend.background = element_blank(),
            legend.text = element_text(),
            legend.title = element_text()
        )+
        scale_color_manual(values=c("#c54e56", "#8398cd"),
                           name = "Fiber type",
                           labels = c("I" = "Type I",
                                      "II" = "Type IIa"))+
        labs(y = "Normalized enrichment score", x="-Log10(p-value)")
}

#' GSEA bubble plot for resistance exercise - downregulated
#'
#' @param gsea_object
#'
#' @return Plot
#' @export
#'
#' @examples
gsea_ter_down <- function(gsea_object) {

    #Make GSEA object dataframe and add desired labels
    df_gsea_ter <- as.data.frame(gsea_ter_both_bp) %>%
        dplyr::rename(fiber_type = Cluster,
                      description = Description) %>%
        dplyr::mutate(label = case_when(
            description == "cellular respiration" ~description,
            description == "mitochondrial respiratory chain complex assembly" ~description,
            description == "cytoskeleton organization" ~description,
            description == "regulation of cell morphogenesis" ~description,
            description == "mitotic cell cycle" ~description,
            description == "proton transmembrane transport" ~description,
            description == "cytoplasmic translation" ~description,
            description == "vesicle−mediated transport" ~description,
            description == "cellular response to lipid" ~description,
            description == "negative regulation of ion transmembrane transporter activity" ~description)) %>%
        dplyr::mutate(log10p = -log10(pvalue), #Add -log10p column
                      count = str_count(.$core_enrichment, '/')+1) %>%
        dplyr::filter(NES < 0)

    #Plot
    ggplot2::ggplot(df_gsea_ter, aes(x = log10p, y = NES, label = label, size = count))+
        geom_point(aes(color = fiber_type))+
        geom_text_repel(point.size=6, size=4, min.segment.length = 0, max.overlaps = Inf, point.padding = 0.2, box.padding = 0.5)+
        theme(
            axis.text.x= element_text(color="black"),
            axis.text.y= element_text(color="black"),
            axis.title.y = element_text(),
            axis.line = element_line(color="black", linewidth = 0.2),
            panel.background = element_rect(color = "black", fill=NA, size = 0.2),
            panel.grid.minor=element_blank(),
            panel.grid.major = element_blank(),
            plot.background = element_blank(),
            #plot.title = element_text(),
            plot.margin = unit(c(5, 5, 5, 25), "pt"),
            text = element_text(family = "Source Sans Pro", size = 15),
            #legend.position = c(0, 1),
            #legend.justification = c(0, 1), #Place legend in top left corner
            legend.background = element_blank(),
            legend.text = element_text(),
            legend.title = element_text()
        )+
        scale_color_manual(values=c("#c54e56", "#8398cd"),
                           name = "Fiber type",
                           labels = c("I" = "Type I",
                                      "II" = "Type IIa"))+
        labs(y = "Normalized enrichment score", x="-Log10(p-value)")
}

#' Title
#'
#' @param df
#' @param protein
#'
#' @return
#' @export
#'
#' @examples
protein_plot <- function(df, protein) {

    df %>%
        dplyr::filter(protein == {{protein}}) %>%
        ggplot(aes(x=fiber_type, y=l2fc, fill=fiber_type))+
        geom_violin(trim = TRUE, width=1, color=NA, alpha=0.5)+
        geom_boxplot(width=0.25, color="black", fill="white", alpha=0.5)+
        geom_jitter(size=5, width=0.03)+
        geom_hline(yintercept=0, linetype="dashed")+
        scale_fill_manual(values=c("#0078b0", "#c41b1b"))+
        scale_x_discrete(labels=c("Type I", "Type II"))+
        theme(panel.background = element_rect(fill = "transparent", colour = NA),
              plot.background = element_rect(fill = "transparent", colour = NA),
              axis.line = element_line(colour = "black"),
              legend.position = "none",
              text = element_text(size = 20),
              axis.text.x= element_text(color="black"),
              axis.text.y= element_text(color="black")
        )+
        ylab("Log2fold change")+
        xlab("Fiber type")+
        facet_wrap(~intervention,
                   labeller=as_labeller(c(resistance="Resistance",
                                          terbutaline="Terbutaline"
                   )))


}

#' Title
#'
#' @param df long format l2fc
#' @param protein protein of choice
#'
#' @return
#' @export
#'
#' @examples
protein_stats <- function(df, protein) {

    #Create dataset
    stats <- df %>%
        dplyr::filter(grepl({{protein}}, protein, ignore.case = T))

    #Fit model
    lm_interaction <- lmer(l2fc ~ fiber_type*intervention + (1|id), data = stats, REML=FALSE)

    fixed_interaction <- anova(lm_interaction)

    qqmath(lm_interaction)

    main_intervention <- emmeans(lm_interaction, specs = pairwise ~ intervention, adjust="none", pbkrtest.limit = 5000, lmerTest.limit = 5000) %>%
        summary(infer=TRUE)

    main_fiber_type <- emmeans(lm_interaction, specs = pairwise ~ fiber_type, adjust="none", pbkrtest.limit = 5000, lmerTest.limit = 5000) %>%
        summary(infer=TRUE)

    interaction <- emmeans::emmeans(lm_interaction, pairwise ~ fiber_type|intervention, adjust="none", pbkrtest.limit = 5000, lmerTest.limit = 5000) %>%
        summary(infer=TRUE)

    #Outputs
    stats <- list("type_iii_fixed"=fixed_interaction,
                  "intervention_em"=as.data.frame(main_fiber_type$emmeans),
                  "intervention_contrasts"=as.data.frame(main_fiber_type$contrasts),
                  "intervention_by_time_em"=as.data.frame(interaction$emmeans),
                  "intervention_by_time_contrasts"=as.data.frame(interaction$contrasts)
    )

    print(stats)
}
