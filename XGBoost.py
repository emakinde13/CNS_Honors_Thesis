# Import Relevant Packages
from xgboost import XGBClassifier
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score, confusion_matrix, accuracy_score, confusion_matrix, precision_recall_curve,average_precision_score, roc_curve, roc_auc_score
from sklearn.impute import KNNImputer
from sklearn.pipeline import Pipeline
from sklearn.base import BaseEstimator, TransformerMixin
import seaborn as sns
from skopt import BayesSearchCV

# Import Data
p45_60_box = pd.read_csv("/Users/emmanuelmakinde/Documents/00_Project_hnRNPD/001_Machine_Learning/Isoform_Specific_Data_Matrices/p45/p45_Camk2a_60_min_shock_box.csv")

# Define Feature Groups
binding = ['n_clusters_x', 'total_ReadCount', 'mean_ModeScore']
regions = ['frac_utr3','frac_intron']
all = ['n_clusters_x', 'total_ReadCount', 'mean_ModeScore', 'frac_utr3','frac_intron', 'frac_utr5', 'frac_cds']
(p45_60_box['Expression'] == 0)

# Outlier Remover
class IQROutlierDetection(BaseEstimator, TransformerMixin):
    def __init__(self, threshold=1.5):
        self.threshold = threshold
    
    def fit(self, X, y=None):
        # Calculate IQR for each column
        self.iqr_ = np.percentile(X, 75, axis=0) - np.percentile(X, 25, axis=0)
        self.lower_bound_ = np.percentile(X, 25, axis=0) - self.threshold * self.iqr_
        self.upper_bound_ = np.percentile(X, 75, axis=0) + self.threshold * self.iqr_
        return self
    
    def transform(self, X):
        # Replace outliers with NaN (anything outside the IQR bounds)
        X_transformed = X.copy()
        for i in range(X.shape[1]):
            outlier_mask = (X[:, i] < self.lower_bound_[i]) | (X[:, i] > self.upper_bound_[i])
            X_transformed[outlier_mask, i] = np.nan
        return X_transformed

# Pipeline
pipeline = Pipeline([
    ('scaler', MinMaxScaler()), # TODO: Why MinMaxScaler()?
    ('imputer', KNNImputer()),
    ('outlier', IQROutlierDetection()),
    ('classifier', XGBClassifier(objective='binary:logistic',
                                 eval_metric='auc',
                                 random_state=42,
                                 scale_pos_weight=4
                                 ))
])

xgb = XGBClassifier
# Hyperparameter tuning with BayesSearchCV
from skopt.space import Real, Integer

search_space = {
    "classifier__n_estimators": Integer(200, 1200),
    "classifier__learning_rate": Real(1e-2, 2e-1, prior="log-uniform"),
    "classifier__max_depth": Integer(2, 8),
    "classifier__min_child_weight": Integer(1, 15),
    "classifier__subsample": Real(0.5, 1.0),
    "classifier__colsample_bytree": Real(0.5, 1.0),
    "classifier__gamma": Real(0.0, 5.0),
    "classifier__reg_alpha": Real(1e-8, 10.0, prior="log-uniform"),
    "classifier__reg_lambda": Real(1e-3, 50.0, prior="log-uniform"),
}

bayes = BayesSearchCV(
    estimator=pipeline,
    search_spaces=search_space,
    n_iter=40,              # bump to 80+ if you can afford it
    scoring="roc_auc",
    cv=5,
    n_jobs=-1,
    random_state=42,
    verbose=1,
    refit=True,
)
X_array = p45_60_box[all].to_numpy()

bayes.fit(X_array, y)

print("Best CV ROC AUC:", bayes.best_score_)
print("Best params:", bayes.best_params_)
best_model = bayes.best_estimator_

# K Folds
n_splits = 10
skf = KFold(n_splits=n_splits, random_state = 42, shuffle=True)
X_array = p45_60_box[all].to_numpy()
y = p45_60_box['Expression']
data_splits_objects = skf.split(X_array,y)

# empty lists to capture the pred and exp
predicted_y = []
expected_y = []
fprs = []
tprs = []

precisions = []
recalls = []

all_y_true = []
all_y_prob = []

average_ap = []

conf_all = []

# Fit iteratively
count = 0
for train_index, test_index in skf.split(X_array, y):
    
    count +=1
    X_train = X_array[train_index]
    X_test  = X_array[test_index]
    y_train = y[train_index]
    y_test  = y[test_index]
    X_train = pd.DataFrame(X_train)
    pipeline.fit(X_train, y_train)
    y_pred = pipeline.predict(X_test)

    predicted_y.extend(y_pred)
    expected_y.extend(y_test)
    accuracy = accuracy_score(y_test, y_pred) 
    conf_m = confusion_matrix(y_test, y_pred)
    # conf_all.extend(conf_m)

    print("---------------------------------------------")
    print(f"Summary for validation {count}/{n_splits}:")
    #Display the accuracy
    print(f'Accuracy: {accuracy:.2f}')
    # Get ROC Metric
    y_pred_prob = pipeline.predict_proba(X_test)[:, 1]
    precision, recall, _ = precision_recall_curve(y_test, y_pred_prob)
    ap = average_precision_score(y_test, y_pred_prob)
    average_ap.append(ap)
    all_y_true.extend(y_test)
    all_y_prob.extend(y_pred_prob)
    print(f'PR-AUC (Average Precision): {ap}')

    precisions.append(precision)
    recalls.append(recall)  
    fpr, tpr, _ = roc_curve(y_test,  y_pred_prob)
    auc = roc_auc_score(y_test, y_pred_prob)
    print(f'AUC: {auc}')
    fprs.append(fpr)
    tprs.append(tpr)

# Plot the PR-AUC
sns.set_theme(style="whitegrid")

plt.figure(figsize=(6,5))
sns.lineplot(x=recall, y=precision)

baseline = np.mean(all_y_true)
plt.axhline(baseline, linestyle="--", color="red",
            label=f"Class prevalence = {baseline:.3f}")

plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title(f"Cross-validated PR curve (AP = {ap:.3f}): DEG vs. nonDEG")

plt.legend()
plt.show()

# Plot the Confusion Matrix
plt.figure(figsize=(6, 6))
sns.heatmap(conf_m, annot=True, fmt="d", cmap="Blues", cbar=True, square=True)
plt.xlabel("Predicted")
plt.ylabel("True")
plt.title("Confusion Matrix: DEG vs. nonDEG")
plt.legend()
plt.show()
