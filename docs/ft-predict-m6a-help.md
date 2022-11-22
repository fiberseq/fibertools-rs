```
Predict m6A positions using HiFi kinetics data and encode the results in the MM and ML bam tags

Usage: ft predict-m6a [OPTIONS] [BAM] [OUT]

Arguments:
  [BAM]  Bam HiFi file with kinetics [default: -]
  [OUT]  Output bam file with m6A calls in new/extended MM and ML bam tags [default: -]

Options:
  -k, --keep        Keep hifi kinetics data
  -c, --cnn         Use CNN model for prediction instead of XGBoost
  -f, --full-float  Add a bam tag (mp) with the full floating point predictions of the ML model
  -h, --help        Print help information
  -V, --version     Print version information

GLOBAL-OPTIONS:
  -t, --threads <THREADS>  Threads for decompression [default: 8]

DEBUG-OPTIONS:
  -v, --verbose...  Logging level [-v: Info, -vv: Debug, -vvv: Trace]
      --quiet       Turn of all logging
```
