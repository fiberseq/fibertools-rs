```
Predict m6A positions using HiFi kinetics data and encode the results in the MM and ML bam tags.
Also adds nucleosome (nl, ns) and MTase sensitive patches (al, as)

Usage: ft predict-m6a [OPTIONS] [BAM] [OUT]

Arguments:
  [BAM]
          Bam HiFi file with kinetics
          
          [default: -]

  [OUT]
          Output bam file with m6A calls in new/extended MM and ML bam tags
          
          [default: -]

Options:
  -n, --nucleosome-length <NUCLEOSOME_LENGTH>
          Minium nucleosome length
          
          [default: 75]

  -c, --combined-nucleosome-length <COMBINED_NUCLEOSOME_LENGTH>
          Minium nucleosome length when combining over a single m6A
          
          [default: 100]

      --min-distance-added <MIN_DISTANCE_ADDED>
          Minium distance needed to add to an already existing nuc by crossing an m6a
          
          [default: 25]

  -d, --distance-from-end <DISTANCE_FROM_END>
          Minimum distance from the end of a fiber to call a nucleosome or MSP
          
          [default: 45]

  -k, --keep
          Keep hifi kinetics data

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version

Developer-Options:
  -m, --min-ml-score <MIN_ML_SCORE>
          Set a minimum ML score to keep on instead of using the model specific minimum ML score

  -a, --all-calls
          Keep all m6A calls regardless of how low the ML value is

      --xgb
          Use the XGBoost model for prediction

      --cnn
          Use the CNN model for prediction

  -s, --semi
          Use the semi-supervised CNN model for prediction [default: true]

  -f, --full-float
          Add a bam tag (mp) with the full floating point predictions of the ML model
          
          For debugging only.

  -b, --batch-size <BATCH_SIZE>
          Number of reads to include in batch prediction
          
          Increasing improves GPU performance at the cost of memory.
          
          [default: 1]

Global-Options:
  -t, --threads <THREADS>
          Threads
          
          [default: 8]

Debug-Options:
  -v, --verbose...
          Logging level [-v: Info, -vv: Debug, -vvv: Trace]

      --quiet
          Turn of all logging
```
