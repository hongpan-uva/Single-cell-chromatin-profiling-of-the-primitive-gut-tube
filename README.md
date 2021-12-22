# Single-cell-chromatin-profiling-of-the-primitive-gut-tube
code for the bioinformatics analyses of "Single-cell chromatin profiling of the primitive gut tube reveals regulatory dynamics underlying lineage fate decisions"

<table>
	<tr>
		<th>figure</th>
		<th>script</th>
 	</tr>
 	<tr>
  		<td>fig 1B</td>
   		<td>scATAC/scATAC_ggplot_umap.R</td>
 	</tr>
	<tr>
		<td rowspan="2">fig 1C</td>
		<td>scATAC/scATAC_ggplot_umap.R</td>
 	</tr>
 	<tr>
   		<td>scRNA/scRNA_markerGene_expression_boxplot.R</td>
 	</tr>
	<tr>
		<td>fig 1D</td>
		<td>scATAC/GREAT_GOterm_dotplot.R</td>
 	</tr>
 	<tr>
  		<td>fig 1E</td>
   		<td>scATAC/Sankey_Diagram.R</td>
 	</tr>
	<tr>
		<td>fig 2A</td>
		<td>Bulk_data/E13.5/E13.5_kmanes_call_difPeak.R</td>
 	</tr>
	<tr>
		<td>fig 2B</td>
		<td>scATAC/scATAC_ArchR_4_peakEnrichment.R</td>
 	</tr>
	<tr>
		<td>fig 2C</td>
		<td>Temporal_Analysis/Association_Score.R</td>
 	</tr>
	<tr>
		<td>fig 2E</td>
		<td>Relative_motif_enrichment_score/motifEnrich_scatterplot_E95E135CrossStages.R</td>
 	</tr>
	<tr>
		<td>fig 2F</td>
		<td>Relative_motif_enrichment_score/motifEnrich_scatterplot_E135E16CrossStages.R</td>
 	</tr>
	<tr>
		<td>fig 3A</td>
		<td>Bulk_data/E13.5/E13.5_difPeak_ATACseq_pattern.R</td>
 	</tr>
	<tr>
		<td>fig 3B, C</td>
		<td>Bulk_data/E13.5/E13.5_difPeak_ChIPseq_pattern.R</td>
 	</tr>
	<tr>
		<td>fig 3D</td>
		<td>Bulk_data/E16.5/E16.5_KO_peak_classification_heatmap.R</td>
 	</tr>
	<tr>
		<td>fig 3E</td>
		<td>scATAC/scATAC_ArchR_5_ChromVar.R</td>
 	</tr>
	<tr>
		<td>fig 3F</td>
		<td>Bulk_data/E16.5/E16.5_KO_peak_classification_barplot.R</td>
 	</tr>
	<tr >
		<td>fig 3H</td>
		<td>Bulk_data/E18.5/4_get_DEG_associatedWithPeaks.R</td>
 	</tr>
	<tr>
		<td rowspan="2">fig 3I</td>
		<td>Bulk_data/E18.5/5_David_GO_dotplot_Cdx2KO.R</td>
 	</tr>
	<tr>
   		<td>Bulk_data/E18.5/6_David_GO_dotplot_WT.R</td>
 	</tr>
	<tr>
		<td>fig 5A, B</td>
		<td>TCGA_data_analysis/TCGA_sigGene_GSEA_refined.R</td>
 	</tr>
	<tr>
		<td>fig 5C, D</td>
		<td>TCGA_data_analysis/TCGA_sigPeak_GSEA_refined.R</td>
 	</tr>
	<tr>
		<td>fig 5E, F, G, H</td>
		<td>TCGA_data_analysis/Bart_result_scatterplot.R</td>
 	</tr>
	<tr>
		<td>fig 6A</td>
		<td>TCGA_data_analysis/TCGA_sigGene_GSEA_refined.R</td>
 	</tr>
	<tr>
		<td>fig 6B</td>
		<td>TCGA_data_analysis/RP_cdfplot.R</td>
 	</tr>
	<tr>
		<td>fig S1B</td>
		<td>scATAC/cell_ranger_ATAC.sh</td>
 	</tr>
	<tr>
		<td>fig S1C</td>
		<td>scATAC/scATAC_ggplot_umap.R</td>
 	</tr>
	<tr>
		<td>fig S1D</td>
		<td>scATAC/scATAC_ArchR_2_recluster.R</td>
 	</tr>
	<tr>
		<td>fig S2A</td>
		<td>scATAC/scATAC_ArchR_3_organSpecificGenes&Peaks.R</td>
 	</tr>
	<tr>
		<td rowspan="2">fig S2B</td>
		<td>scRNA/scRNA_seurat_pipeline.R</td>
 	</tr>
	<tr>
		<td>scATAC/scATAC_ggplot_umap.R</td>
 	</tr>
	<tr>
		<td>fig S2C</td>
		<td>scATAC/scATAC_ArchR_1_cluster&integration.R</td>
 	</tr>
	<tr>
		<td>fig S3A</td>
		<td>scATAC/scATAC_ggplot_umap.R</td>
 	</tr>
	<tr>
		<td>fig S3B</td>
		<td>scATAC/scATAC_ArchR_1_cluster&integration.R</td>
 	</tr>
	<tr>
		<td>fig S3C</td>
		<td>scATAC/cell_label_percentage_stack_barplot.R</td>
 	</tr>
	<tr>
		<td>fig S3D</td>
		<td>scATAC/cell_label_replicate_association_sankey_diagram.R</td>
 	</tr>
	<tr>
		<td>fig S4</td>
		<td>Temporal_Analysis/E95_E135_FPKM_barplot.R</td>
 	</tr>
	<tr>
		<td>fig S5</td>
		<td>scATAC/E13.5_kmeans_call_DEG.R</td>
 	</tr>
	<tr>
		<td>fig S6A</td>
		<td>scATAC/commonPeaks_foregut/3_E95foregut_commonPeaks_composition.R</td>
 	</tr>
	<tr>
		<td>fig S6B</td>
		<td>scATAC/commonPeaks_foregut/4_E95foregut_commonGenes_composition_byFC.R</td>
 	</tr>
	<tr>
		<td>fig S7</td>
		<td>Temporal_Analysis/Relative_motif_enrichment_score/with_expression_data/motifEnrich_scatterplot_E95E135CrossStages_byExpressionlogFC.R</td>
 	</tr>
	<tr>
		<td>fig S8</td>
		<td>Temporal_Analysis/Relative_motif_enrichment_score/with_expression_data/motifEnrich_scatterplot_E135E16CrossStages_byExpressionlogFC.R</td>
 	</tr>
	<tr>
		<td>fig S9A</td>
		<td>Bulk_data/E16.5/E16.5_KO_peak_classification_heatmap.R</td>
 	</tr>
	<tr>
		<td>fig S9B</td>
		<td>scATAC/scATAC_ArchR_5_ChromVar.R</td>
 	</tr>
	<tr>
		<td>fig S9C</td>
		<td>Bulk_data/E16.5/E16.5_KO_peak_classification_barplot.R</td>
 	</tr>
	<tr>
		<td>fig 11A, B</td>
		<td>TCGA_data_analysis/TCGA_sigGene_GSEA.R</td>
 	</tr>
	<tr>
		<td>fig 11C, D</td>
		<td>TCGA_data_analysis/TCGA_sigPeak_GSEA.R</td>
 	</tr>
	<tr>
		<td>fig 11E</td>
		<td>TCGA_data_analysis/Bart_result_scatterplot.R</td>
 	</tr>
	<tr>
		<td>fig 11F</td>
		<td>TCGA_data_analysis/Bart_result_scatterplot.R</td>
 	</tr>
	<tr>
		<td>fig 11G</td>
		<td>TCGA_data_analysis/RP_cdfplot.R</td>
 	</tr>
</table>
