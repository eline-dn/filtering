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
# #To handle nan and multiple measurments:
# remove binders with no binding information (boolean) 

# remove binders that did not express "expressed" column, or expressed inconsistentely (true & false)

# or with unsure/inconsistent binding to the target: (true & false in the same target/binder pair)

# for binders with several measurments of binding (check same number of measurments in binding, binding_target, and KD column), keep common value in booleans, binding target and average kd


      
