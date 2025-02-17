# NMF-Gene-Clustering
Non-Negative Matrix Factorization for Gene Expression Clustering
Berezin Lab, Washington University in St. Louis, 2024
Overview
This MATLAB function applies Non-Negative Matrix Factorization (NMF) to gene expression data to identify clusters or latent factors that reveal how groups of genes are commonly regulated in response to oxaliplatin treatment. The function enables:

  •	Identification of Gene Clusters (Basis Matrix W): Groups of genes that are co-expressed or have similar mean expression levels across different samples.
  
  •	Identification of Mouse Clusters (Loading Matrix H): How each mouse aligns with gene clusters, distinguishing between treatment and control groups.
  
  •	Interpretation of Biological Mechanisms: The identified gene clusters can be analyzed using pathway enrichment or Gene Ontology (GO) analysis to uncover biological functions affected by oxaliplatin.
  
Input Requirements

  •	The input data should be a matrix of samples (rows) by genes (columns).
  
  •	The function requires an Excel file (e.g., drg_oxa_89_genes.xlsx) where: 
  
    o	The first row contains gene names.
    
    o	The first column contains sample identifiers.
    
    o	The rest of the matrix consists of numerical expression values.
    
Installation
1.	Clone this repository: 
2.	git clone https://github.com/MikhailBerezin/NMF-Gene-Clustering 
3.	Open MATLAB and navigate to the project folder.
4.	Ensure you have the required toolboxes (MATLAB Statistics and Machine Learning Toolbox).
Usage
1.	Run the function in MATLAB: 
2.	nnmf_rna
3.	A file selection dialog will appear. Select the Excel file containing gene expression data.
4.	Enter the number of clusters when prompted (default: 2).
5.	The function performs NMF decomposition and generates:
6.	
  o	Heatmaps for Basis Matrix W (Gene Clusters)

  o	Heatmaps for Loading Matrix H (Mouse Clusters)
  
  o	Bar plots showing treatment association with latent factors
  
  o	Top 10 genes per cluster with the highest loading coefficients
Outputs

  •	Basis Matrix (W): Identifies gene clusters with similar expression patterns.
  
  •	Loading Matrix (H): Shows how each sample aligns with gene clusters.
  
  •	Visualization Plots: Heatmaps and bar plots for easy interpretation.
  
  •	Top Genes per Cluster: Displays the top 10 genes in each cluster.
  
Interpretation
  •	Genes in the same cluster may indicate co-regulated pathways.
  
  •	Mice with similar loading patterns might share transcriptional responses.
  
  •	Pathway enrichment analysis can be performed on clusters to uncover biological relevance (e.g., DNA repair, apoptosis, cell cycle regulation).
  
Dependencies
  •	MATLAB (R2021a or newer recommended)
  
  •	MATLAB Statistics and Machine Learning Toolbox
  
Citation
If you use this function in your research, please cite: Berezin Lab, Washington University in St. Louis, 2024
License
This project is licensed under the BSD 2-Clause License.
________________________________________
For any questions or contributions, feel free to create an issue or submit a pull request!

