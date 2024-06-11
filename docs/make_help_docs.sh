#!/usr/bin/env bash

echo $LIBTORCH
export LD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH

out="src/help.md"
echo "# Help pages for fibertools subcommands" >$out
echo "<!-- toc -->" >>$out
echo "" >>$out

for subcommand in "" "predict-m6a" "fire" "extract" "center" "add-nucleosomes" "footprint" "clear-kinetics" "strip-basemods" "track-decorators" "pileup"; do
    echo $subcommand
    #out="docs/ft-${subcommand}-help.md"

    echo "## \`ft $subcommand\`" >>$out
    echo '```console' >>$out
    cargo run -- $subcommand -h >>$out
    echo '```' >>$out
    echo "" >>$out
done
