import xgboost as xgb
import numpy as np

train_csv = "models/temp/train.csv"
val_csv = "models/temp/val.csv"

read_data = False

if read_data:
    data_path = "models/temp/m6A_train_other_half_hifi.npz"
    train_val_data = np.load(data_path, allow_pickle=True)

    # Get the dictionary from the containing relevant data
    train_val_data = train_val_data["save_data_dict"][()]

    # Load training and validation features and labels
    X_train = train_val_data["X_train"]
    y_train = train_val_data["y_train"]
    X_val = train_val_data["X_val"]
    y_val = train_val_data["y_val"]

    # Reshape features to one dimension
    X_train = X_train.reshape(X_train.shape[0], X_train.shape[1] * X_train.shape[2])
    X_val = X_val.reshape(X_val.shape[0], X_val.shape[1] * X_val.shape[2])

    # Add labels before feature list
    y_train = y_train[:, 0]
    y_val = y_val[:, 0]
    final_train_data = np.concatenate([y_train[:, np.newaxis], X_train], axis=1)
    final_val_data = np.concatenate([y_val[:, np.newaxis], X_val], axis=1)
    print("saving")
    # Save data as CSV
    np.savetxt(train_csv, final_train_data, delimiter=",")
    np.savetxt(val_csv, final_val_data, delimiter=",")
    print("saved")
elif False:
    train_csv = "models/temp/train_large.csv"
    data_path = "models/temp/large.PS00075.npz"
    train_val_data = np.load(data_path, allow_pickle=True)

    # Load training and validation features and labels
    X_train = train_val_data["features"]
    y_train = train_val_data["labels"]

    # Reshape features to one dimension
    X_train = X_train.reshape(X_train.shape[0], X_train.shape[1] * X_train.shape[2])

    # Add labels before feature list
    final_train_data = np.concatenate([y_train[:, np.newaxis], X_train], axis=1)

    data_path = "models/temp/PS00075_2.val.npz"
    train_val_data = np.load(data_path, allow_pickle=True)
    X_val = train_val_data["features"]
    y_val = train_val_data["labels"]
    X_val = X_val.reshape(X_val.shape[0], X_val.shape[1] * X_val.shape[2])
    final_val_data = np.concatenate([y_val[:, np.newaxis], X_val], axis=1)

    # Read in training and validation data
    dtrain = xgb.DMatrix(X_train, label=y_train)
    dval = xgb.DMatrix(X_val, label=y_val)

    # Save data as CSV
    print("saving")
    np.savetxt(train_csv, final_train_data, delimiter=",")
    np.savetxt(val_csv, final_val_data, delimiter=",")
    print("saved")

# Read in training and validation data
dtrain = xgb.DMatrix(f"{train_csv}?format=csv&label_column=0")
dval = xgb.DMatrix(f"{val_csv}?format=csv&label_column=0")

# Specify parameters via map
param = {
    "max_depth": 6,
    "objective": "binary:logistic",
    "eval_metric": ["aucpr", "auc"],
}
num_round = 100
eval_s = [(dtrain, "train"), (dval, "validation")]

# Train model
bst = xgb.train(param, dtrain, num_round, evals=eval_s)

# Save model
bst.save_model("models/xgboost_large_2.0.81.bin")
