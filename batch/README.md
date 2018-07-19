## Batch Data Processing

The script in this folder can be used to achieve batch data processing with Google's Pipelines API. 

Usage:
```
bash submit_batch.sh batch.json batch.tsv
```

Options that apply to all samples in the batch are read from the JSON file. The tsv file (tab-separated) can be used to modify specific fields for each sample. In `submit_batch.sh`, the variable `n_concurrent` controls the number of concurrent jobs while `polling_interval` controls the polling interval. 

