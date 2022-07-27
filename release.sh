#!/bin/bash
set -euxo pipefail
V=$(grep "^version" Cargo.toml | sed 's/.*= //g' | sed 's/"//g' | xargs)
echo $V

cargo clippy
cargo test
cargo run -- --help

exit
# broken for some reason
# target=aarch64-unknown-linux-musl
# target=x86_64-unknown-linux-gnu

mkdir -p dists

for target in x86_64-unknown-linux-musl; do
    echo $target
    cross build --release --target ${target}
    tar -czvf ./dists/fibertools-rs_v${V}-${target}.tar.gz \
        -C ./target/${target}/release/ \
        ft
done

cargo publish

gh release create \
    "v${V}" \
    -t "v${V}" \
    -n "v${V}" \
    ./dists/fibertools-rs_v${V}-*tar.gz
