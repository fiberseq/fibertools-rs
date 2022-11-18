```
ft-predict-m6a 0.0.8
Predict m6A positions using HiFi kinetics data and encode the results in the MM and ML bam tags

USAGE:
    ft predict-m6a [OPTIONS] <BAM> <OUT>

ARGS:
    <BAM>    HiFi file with kinetics
    <OUT>    Output bam file with m6A calls in new/extended MM and ML bam tags

OPTIONS:
    -k, --keep          keep hifi kinetics data
    -c, --cnn           use CNN model for prediction instead of XGB
    -f, --full-float    Add a bam tag (mp) with the full floating point predictions of the ML model
    -h, --help          Print help information
    -V, --version       Print version information
```
