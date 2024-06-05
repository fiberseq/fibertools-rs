# With libtorch

Get `libtorch` **v2.2.0** from the [PyTorch website](https://pytorch.org/get-started/) and extract the content of the zip file.

- On Linux/Unix system you can download with:
  - `wget  https://download.pytorch.org/libtorch/cu118/libtorch-shared-with-deps-2.2.0%2Bcu118.zip`
- On macOS you can download with:
  - `wget https://download.pytorch.org/libtorch/cpu/libtorch-macos-2.2.0.zip`
- Windows is not supported and will not be.

Then add the following to your `.bashrc` or equivalent, where `/path/to/libtorch` is the path to the directory that was created when unzipping the file:

```bash
export LIBTORCH_CXX11_ABI=0
export LIBTORCH=/path/to/libtorch # e.g. export LIBTORCH=/Users/mrvollger/lib/libtorch
export LD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH
```

Finally install using `cargo`:
```bash
cargo install --all-features fibertools-rs
```

