You'll need zlib installed - `sudo apt install zlib1g-dev` on Debian/Ubuntu.

Then just `make`.

Run as `./fastqst_possible_samp r1_file.fastq.gz r2_file.fastq.gz config.csv <total_reads_in_input_fastq>`

Where config.csv is formatted like:
```
an_output_base_path_,1000
another_output_base_path_,50000
...
```
The above config file would give you 1000 samples in an_output_base_path_r1.fastq.gz and an_output_base_path_r2.fastq.gz,
50000 samples in another_output_base_path_r1.fastq.gz and another_output_base_path_r2.fastq.gz, etc.
