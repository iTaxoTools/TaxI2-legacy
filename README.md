# TaxI2
Calculates genetic differences between DNA sequences

## Scores for alignment
The file `data/scores.tab` contains the scores used in the sequence alignment.
Each line has the format:
```
score_identifier<Tab>value
```

The possible scores are:
* `gap penalty`: Score to open a gap in the middle of a sequence
* `gap extend penalty`: Score to extend an existing gap in the middle of a sequence
* `end gap penalty`: Score to create a gap at an end of a sequence.
* `end gap extend penalty`: Score to extend a gap at an end of a sequence.
* `match score`: Score for matching nucleotides
* `mismatch score`: Score for non-matching nucleotides

## Compiling the Rust module

TaxI2 included a module for calculating distances between sequences, which is written in Rust. It might need to be compiled before it can be used.

### Rust installation
[Instructions on the Rust site](https://www.rust-lang.org/tools/install).

Detailed instruction below.

#### Linux installation

```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
rustup update
```

#### Windows installation
Download and run [rustup-init.exe](https://static.rust-lang.org/rustup/dist/i686-pc-windows-gnu/rustup-init.exe).

You might need to install [Microsoft Build Tools for C++](https://visualstudio.microsoft.com/de/visual-cpp-build-tools), which should include the c++ compiler, Windows 10 SDK and English language files.

Find out whether your Python interpreter is 32- or 64-bit. For this run:
```
python
>>> import platform
>>> platfom.architecture
(bits, linkage)
```
The `bits` will contain the required information.

If Python is 32-bit, do:
```
rustup default stable-i686-pc-windows-msvc
```

If Python is 64-bit, do:
```
rustup default stable-x86_64-pc-windows-msvc
```

### Compilation

In the directory `library/calculate_distances` run:
```
cargo build --release
```

To compile for another Python interpreter (for example `python3.6`), do instead:
* On Linux:
```
PYO3_PYTHON=python3.6 cargo build --release
```
* On Windows:
```
$env:PYO3_PYTHON=python3.6 
cargo build --release
```

In the directory `library` run:
```
python copy_rust_lib.py
```

