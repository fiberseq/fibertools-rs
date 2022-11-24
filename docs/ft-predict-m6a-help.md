```
Predict m6A positions using HiFi kinetics data and encode the results in the MM and ML bam tags

Usage: ft predict-m6a [OPTIONS] [BAM] [OUT]

Arguments:
  [BAM]
          Bam HiFi file with kinetics

          [default: -]

  [OUT]
          Output bam file with m6A calls in new/extended MM and ML bam tags

          [default: -]

Options:
  -k, --keep
          Keep hifi kinetics data

  -c, --cnn
          Use CNN model for prediction instead of XGBoost

  -f, --full-float
          Add a bam tag (mp) with the full floating point predictions of the ML model

  -b, --batch-size <BATCH_SIZE>
          Number of reads to include in batch prediction

          Increasing improved GPU performance at the cost of memory.

          [default: 1]

  -h, --help
          Print help information (use `-h` for a summary)

  -V, --version
          Print version information

Global-Options:
  -t, --threads <THREADS>
          Threads for decompression

          [default: 8]

Debug-Options:
  -v, --verbose...
          Logging level [-v: Info, -vv: Debug, -vvv: Trace]

      --quiet
          Turn of all logging
```
