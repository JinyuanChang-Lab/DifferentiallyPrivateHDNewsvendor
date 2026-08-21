# Processed Data

This directory is used to store the analysis data generated from the original
Corporación Favorita files.

Run the following command from the root directory of the replication package:

```bash
Rscript scripts/real_data/prepare_real_data.R
```

The processing script generates `data/processed/processed_data.rds`.

The generated file contains the objects used as inputs to the empirical
analysis, including:

- `RD_X`: the predictor matrix;
- `RD_Y`: the response variable.

The processed file can be deleted and regenerated from the original data at
any time by rerunning the processing script.