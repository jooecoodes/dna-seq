import joblib
import numpy as np

# Load trained model
clf = joblib.load("best_algo_decision_tree.pkl")

# Example: new pattern features
# (length, gc_content, entropy)
new_pattern = np.array([[64,0.8,1.9]])

# Predict best algorithm
prediction = clf.predict(new_pattern)
print("🚀 Best algorithm:", prediction[0])