This is the command used on Google BigQuery to generate all the tables used on the analysis globally. 
They were downloaded as compressed CSV files and combined into included Parquet files. 
The downstream files were processed using DuckDB SQL engine API through R on a laptop on disk
to analyze data that cannot be fit in RAM.
