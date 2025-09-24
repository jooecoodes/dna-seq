import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier, export_text
from sklearn.metrics import classification_report, accuracy_score

# Load labeled dataset
data = pd.read_csv("training_data/labeled_best_algorithms.csv")

# Features: only DNA-intrinsic properties
X = data[["length", "gc_content", "entropy"]]

# Target: best algorithm
y = data["best_algorithm"]

# Split into training & test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Train lightweight Decision Tree
clf = DecisionTreeClassifier(max_depth=8, random_state=42)
clf.fit(X_train, y_train)

# Evaluate
y_pred = clf.predict(X_test)
print("✅ Accuracy:", accuracy_score(y_test, y_pred))
print("\nClassification Report:\n", classification_report(y_test, y_pred))

# Print tree rules for interpretability
print("\n=== Decision Tree Rules ===")
print(export_text(clf, feature_names=list(X.columns)))

# Save model (optional)
import joblib
joblib.dump(clf, "best_algo_decision_tree.pkl")
print("\n💾 Model saved to models/best_algo_decision_tree.pkl")
