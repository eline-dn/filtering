# goal: compute our set of metrics on experimentally tested  de novo binders and non-binders to compare with our selection metrics.

# load dependencies
import pandas as pd
import json
import numpy as np

# step 1 process proteinbase db
#db=pd.read_csv("proteinbase_all_data.csv")
#db.drop(label=["name","author"], axis=1)
# parse experimental data column

## ------------------------------------------------------------protein database parsing (esp evaluation col)------------------------------------
METRICS_SCHEMA = {
    "binding": {
        "valueType": "boolean",
        "dtype": "boolean",
    },
    "binding_strength": {
        "valueType": "label",
        "dtype": "string",
    },
    "expressed": {
        "valueType": "boolean",
        "dtype": "boolean",
    },
    "kd": {
        "valueType": "numeric",
        "dtype": "float",
    } 
    """"seqidentityafdb50": {
        "valueType": "numeric",
        "dtype": "float",
    },"""
}


df = pd.read_csv("proteinbase.csv")

# Parse JSON array column
df["evaluations"] = df["evaluations"].apply(json.loads)


# helper:
def extract_metric_list(evaluations, metric_name, expected_value_type):
    """
    Returns a list of values for a given metric.
    If the metric is absent, returns np.nan (not an empty list).
    """
    values = [
        ev.get("value")
        for ev in evaluations
        if ev.get("metric") == metric_name
        and ev.get("valueType") == expected_value_type
    ]
    return values if values else np.nan


for metric, spec in METRICS_SCHEMA.items():
    df[metric] = df["evaluations"].apply(
        extract_metric_list,
        metric_name=metric,
        expected_value_type=spec["valueType"]
    )

# handle target name
def extract_binding_target_list(evaluations):
    targets = [
        ev.get("target")
        for ev in evaluations
        if ev.get("metric") == "binding"
        and "target" in ev
    ]
    return targets if targets else np.nan

df["binding_target"] = df["evaluations"].apply(extract_binding_target_list)

### --------------------------clean the dataset-----------------
# #Goal:  handle nan and multiple measurments:
# helpeers:
def is_listlike(x):
    return isinstance(x, (list, tuple, np.ndarray))

def ensure_list(x):
    if pd.isna(x):
        return []
    return list(x) if is_listlike(x) else [x]

def all_equal(values):
    return len(set(values)) == 1 if values else True

# remove binders with no binding information (boolean) 
# remove binders that did not express ("expressed" column), or expressed inconsistentely (true & false) => we only keep binders that expressed all the time
# or with unsure/inconsistent binding to the target: (true & false in the same target/binder pair) => only binds all the time to one target, or doesn't bind at all
# for binders with several measurments of binding (check same number of measurments in binding, binding_target, and KD column), keep common value in booleans, binding target and average the kd values

### -------------------------- clean the dataset -----------------

# Ensure list-like consistency
for col in ["binding", "binding_target", "kd", "expressed"]:
    if col in df.columns:
        df[col] = df[col].apply(ensure_list)


# 1. Remove binders with no binding  or expression information at all
df = df[df["binding"].apply(len) > 0]
df = df[df["expressed"].apply(len) > 0]


# 2. Remove binders that did not express
#    or expressed inconsistently (True & False)
def expressed_consistent_and_true(expr_list):
    return (
        len(expr_list) > 0
        and all_equal(expr_list)
        and expr_list[0] is True
    )

df = df[df["expressed"].apply(expressed_consistent_and_true)]


# 3. Remove binders with inconsistent binding to the same target
#    Rule:
#    - either always binds (True only)
#    - or never binds (False only)
def binding_consistent(binding_list):
    return len(binding_list) > 0 and all_equal(binding_list)

df = df[df["binding"].apply(binding_consistent)]


# 4. Handle multiple measurements
#    - binding: keep the common boolean
#    - binding_target: must be identical → keep value
#    - kd: average values (if present)

def reduce_row(row):
    out = row.copy()

    # binding
    out["binding"] = row["binding"][0] if row["binding"] else np.nan # np.nan case should not exist because of previous filtering

    # binding_target
    if row["binding_target"]:
        if not all_equal(row["binding_target"]):
            return None  # inconsistent targets → drop
        out["binding_target"] = row["binding_target"][0]
    else:
        out["binding_target"] = np.nan

    # kd
    out["kd"] = (
        float(np.mean(row["kd"])) if row["kd"] else np.nan
    )

    # expressed (already consistent & true by construction)
    out["expressed"] = True

    return out


df = df.apply(reduce_row, axis=1)
df = df[df.notna().all(axis=1)]  # drop rows rejected above

      
