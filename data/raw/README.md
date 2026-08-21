# Raw data

The empirical analysis uses data from the Kaggle competition
“Corporación Favorita Grocery Sales Forecasting.”

The original competition files are not redistributed in this repository.

## Download Instructions

1. Download the competition data from the [Corporación Favorita Grocery Sales Forecasting data page](https://www.kaggle.com/competitions/favorita-grocery-sales-forecasting/data). A Kaggle account and acceptance of the competition rules may be required.


2. Extract the following files into `data/raw/`:
   - `holidays_events.csv`
   - `items.csv`
   - `oil.csv`
   - `stores.csv`
   - `train.csv`
   - `transactions.csv`
   
   After extraction, the directory should contain:

```text
data/raw/
├── holidays_events.csv
├── items.csv
├── oil.csv
├── stores.csv
├── train.csv
├── transactions.csv
└── README.md
```

## Data Processing
From the root directory of the replication repository, run:

```bash
Rscript scripts/real_data/prepare_real_data.R
```

The processing script reads the original files directly and generates `data/processed/processed_data.rds`.

No manual modification of the original CSV
files is required.
