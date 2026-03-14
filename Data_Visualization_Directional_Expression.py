# This is data visualization steps for the directional expression. This is applicable to the differential expression. Instead of indexing the "Direction" column which contains the labels for the direction, you index "Expression" column, or the one which contains the expression

# This data contains all the targets of hnRNPD found in the TRAP-seq set for 60min and 15min timepoints for both shock/hc and shock/box sets.
p45_60_box = pd.read_csv("/Users/emmanuelmakinde/Documents/00_Project_hnRNPD/001_Machine_Learning/Isoform_Specific_Data_Matrices/p45/p45_Camk2a_60_min_shock_box_direction.csv")

# One Hot Encoding the Labels
le = LabelEncoder()
p45_60_box['Direction'] = le.fit_transform(p45_60_box['Direction'])

# Define the feature groups
binding = ['n_clusters_x', 'total_ReadCount', 'mean_ModeScore']
regions = ['frac_utr3','frac_intron', 'frac_utr5', 'frac_cds']
all = ['n_clusters_x', 'total_ReadCount', 'mean_ModeScore','frac_utr3','frac_intron']

# PRINCIPAL COMPONENT ANALYSIS
X = p45_60_box[all]
X = StandardScaler().fit_transform(X)

pca = PCA(n_components=3)

principalComponents = pca.fit_transform(X)

principalDf = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2', 'principal component 3'])

finalDf = pd.concat([principalDf, p45_60_box[['Direction']]], axis = 1)

fig = plt.figure(figsize = (4,4))
ax = fig.add_subplot(1,1,1, projection='3d') 
ax.set_xlabel(f'Principal Component 1:  {pca.explained_variance_ratio_[0].round(3)}', fontsize = 7)
ax.set_ylabel(f'Principal Component 2:  {pca.explained_variance_ratio_[1].round(3)}', fontsize = 7)
ax.set_zlabel(f'Principal Component 3:  {pca.explained_variance_ratio_[2].round(3)}', fontsize = 7)
ax.set_title('PCA: Direction (UP vs. Down)', fontsize = 10)

print(pca.explained_variance_ratio_)

targets = [1, 0]
colors = ['r', 'g']
for target, color in zip(targets, colors):
    indicesToKeep = finalDf['Direction'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2'],
               finalDf.loc[indicesToKeep, 'principal component 3']

               , c = color
               , s = 25)
ax.legend(targets)
ax.grid()

# CORRELATION MATRIX
correlation_matrix = p45_60_box_ALL.corr(numeric_only = True)
plt.figure(figsize=(8, 6))
sns.heatmap(correlation_matrix, annot=False, cmap='viridis', linewidths=0)
plt.title("Correlation Heatmap: Binding Confidence")
plt.show()
