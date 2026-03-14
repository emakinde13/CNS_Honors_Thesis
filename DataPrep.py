# Function which turns clusters into bed files
def turn_into_bed(cluster_file,name):
    bed_df = cluster_file[['Chromosome','ClusterStart','ClusterEnd','ClusterID','ReadCount','Strand']]
    bed_df.columns = ['Chromosome', 'Start', 'End','Name','Score','Strand'] # chrom chromStart chromEnd name score strand
    bf = pybed.BedFrame.from_frame(meta=[], data=bed_df)
    bf.to_file(f"{name}.bed")
    return bf
  
# Define a function which calculates the fractions of region type the clusters map onto (introns, 3' UTR, 5' UTR, and CDS)
def __calculate_fractions__(dataframe):
    """
    calculates the fractions of the region type the cluster maps onto
    dataframe: type pandas dataframe
    """
    # Using groupby to collapse to a gene level series based on count of region_type
    region_counts = dataframe.groupby(["Ensembl_ID", "Region_Type"]).size()

    # Converting that into a dataframe. Fill NaNs with 0.
    region_counts_df = region_counts.unstack().fillna(0)

    region_cols = ["utr3","utr5","cds","intron"]
    region_counts_df["n_clusters"] = region_counts_df[region_cols].sum(axis=1)
    
    # Calculate the fractions by doing column wide division. ***Almost used a for loop here.
    region_counts_df["frac_utr3"] = region_counts_df["utr3"] / region_counts_df["n_clusters"]
    region_counts_df["frac_utr5"] = region_counts_df["utr5"] / region_counts_df["n_clusters"]
    region_counts_df["frac_cds"]  = region_counts_df["cds"]  / region_counts_df["n_clusters"]
    region_counts_df["frac_intron"] = region_counts_df["intron"] / region_counts_df["n_clusters"]

    # Safeguard against genes without a given region
    for col in ["utr3","utr5","cds","intron"]:
        if col not in region_counts_df:
            region_counts_df[col] = 0
    # --- Failed Code --- Merging with cluster ID bad after gene level aggregation.
    # ensembl_id_df = dataframe[['ClusterID',"Ensembl_ID"]]
    # print(region_counts_df.head())
    # final_region_counts_df = pd.merge(region_counts_df, ensembl_id_df, on="Ensembl_ID")
    return region_counts_df

# Define a function to aggregate total read count, mean mode score, and n clusters by gene.
def __binding_df__(dataframe):
    binding_df = (
    dataframe
    .groupby("Ensembl_ID")
    .agg(
        total_ReadCount=("ReadCount", "sum"),
        mean_ModeScore=("ModeScore", "mean"),
        n_clusters=("ClusterID", "nunique"))
        )
    return binding_df

# Organize the Columns
def __change_columns__(intersection_data):
    intersection_data.columns = ["Chromosome","Start","Stop",
                                 "ClusterID","Score","Strand",
                                 "Chromosome_ref","Start_ref","Stop_ref",
                                 "Ensembl_ID", "Score_ref","Strand_ref",
                                 "HAVANA","gene","blank","ID","Overlap_Count",
                                 "Region_Type"]
# Define a function which returns the final data matrix
def __create_master_mat__(intersection_data, cluster_data, isoform_name):
    __change_columns__(intersection_data)
    data_orig_file = cluster_data[["ClusterID","ReadCount","ModeScore"]]
    init_data_mat = pd.merge(data_orig_file, intersection_data, on = ["ClusterID"])
    cis_data_mat = __calculate_fractions__(init_data_mat)
    binding_mat = __binding_df__(init_data_mat)
    final_df = pd.merge(cis_data_mat, binding_mat, on="Ensembl_ID")
    final_df["Isoform"] = isoform_name
    return final_df[['n_clusters_x','frac_utr3','frac_utr5','frac_cds','frac_intron','total_ReadCount','mean_ModeScore','Isoform']]

# Define a function that recreate the original complete dataset of TRAP-genes now with the Ensembl IDs
def _merge_ensembls(all_genes, ensembl_id): 
    all_genes = all_genes.rename(columns = {all_genes.columns[0]:'GeneSymbol'})
    return pd.merge(all_genes, ensembl_id,on="GeneSymbol")

# Define a filter function to filter the TRAP genes based on pval or padj
def _filter(merge_matrix,filter='padj',threshold=0.05):
  """
  THIS IS FUNCTION IS ONLY TO BE CALLED TO CREATE A MATRIX OF DATA WITH DEGs (as opposed to DEGs & non-DEGs)
  """
    if filter == 'pvalue':
        df = merge_matrix[["Ensembl_ID","log2FoldChange","pvalue"]]
        filtered_df = df[df['pvalue'] < threshold]
    elif filter == 'padj':
        df = merge_matrix[["Ensembl_ID","log2FoldChange","padj"]]
        filtered_df = df[df['padj'] < threshold]
    else:
        raise ValueError(f"{filter} is an unrecognized. Did you mean 'padj' or 'pvalue'?") # just practicing catching errors
    return filtered_df
# Define a function to create a new column with a 1 or 0 on DEA based on log2FC.
def _binarize_genes(df,threshold=0.5):
    df['Expression'] = np.where((np.abs(df['log2FoldChange']) > threshold), 1, 0)
    new_df = df[["Ensembl_ID","Expression"]]
    new_df["Ensembl_ID"] = new_df["Ensembl_ID"].astype(str).str.split(".").str[0]
    return new_df
  
# Define a function to create a new column with a 1 or 0 on DEA based on PVALUE.
def _binarize_genes_(df, threshold=0.05):
    df = df.copy()  # ← prevent mutation of the caller's dataframe
    df['Expression'] = np.where(df['padj'] < threshold, 1, 0) 
    new_df = df[["Ensembl_ID", "Expression"]].copy()  
    new_df["Ensembl_ID"] = new_df["Ensembl_ID"].astype(str).str.split(".").str[0]
    new_df = new_df.drop_duplicates(subset="Ensembl_ID", keep="first")  # ← remove version conflicts

    return new_df

# Define a function to create a data matrix
def _create_data_mat(clip_mat, trap_mat):
    clip = clip_mat.copy()
    """Create a merged pandas dataframe on the Ensembl IDs between the CLIP-data and TRAP-data

    Parameters:
    clip_mat (pd.Dataframe): CLIP-seq data with Ensembl_ID key. The function also checks for this column label and adds if absent.
    trap_mat (pd.Dataframe): TRAP-seq data with Ensembl_ID key

    Returns:
    pd.Dataframe :Returning value
    """
    # If Ensembl_ID is the index, bring it back as a column
    if "Ensembl_ID" not in clip.columns:
        if clip.index.name == "Ensembl_ID":
            clip = clip.reset_index()
        else:
            # reset anyway and inspect what column name we got
            clip = clip.reset_index()

    # After reset_index, sometimes the column is called 'index'
    if "Ensembl_ID" not in clip.columns and "index" in clip.columns:
        clip = clip.rename(columns={"index": "Ensembl_ID"})

    # Normalize Ensembl IDs on the CLIP side
    clip["Ensembl_ID"] = clip["Ensembl_ID"].astype(str).str.strip().str.replace(r"\.\d+$", "", regex=True)

    # trap_mat should already be normalized, but normalize again safely
    trap = trap_mat.copy()
    trap["Ensembl_ID"] = trap["Ensembl_ID"].astype(str).str.strip().str.replace(r"\.\d+$", "", regex=True)
    return pd.merge(clip, trap, on="Ensembl_ID", how="inner")

# Define a function to run the pipeline
def _run_all(all_genes, ensembl_id,clip_mat ,p_threshold=0.05):
    return _create_data_mat(clip_mat, _binarize_genes_(_merge_ensembls(all_genes, ensembl_id), p_threshold))

# Define a function to create a new column with a 1 or 0 on DEA based on PVALUE.
def _binarize_genes_direction_(df, threshold=0.05):
    df = df.copy()  # ← prevent mutation of the caller's dataframe
    df['Expression'] = np.where(df['padj'] < threshold, 1,0) # Upregulated means 1
    df['Direction'] = "NoChange"
    for gene in range(len(df['Expression'])):
        if df['Expression'][gene] == 1:
            if df['log2FoldChange'][gene] > 0:
                df['Direction'][gene] = "UP"
            if df['log2FoldChange'][gene] < 0 :
                df['Direction'][gene] = "DOWN"
    # labels = df[df.Direction == 'UP' and df.Direction == 'DOWN']
    new_df = df[["Ensembl_ID", "Expression", "Direction","log2FoldChange"]].copy()  
    new_df["Ensembl_ID"] = new_df["Ensembl_ID"].astype(str).str.split(".").str[0]
    new_df = new_df.drop_duplicates(subset="Ensembl_ID", keep="first")  # ← remove version conflicts
    new_df = new_df[~(new_df["Direction"]=="NoChange")]

    return new_df

# Define a function to run the pipeline
def _run_all__(all_genes, ensembl_id,clip_mat ,p_threshold=0.05):
    return _create_data_mat(clip_mat, _binarize_genes_direction_(_merge_ensembls(all_genes, ensembl_id), p_threshold))
