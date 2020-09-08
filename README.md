# Squiggler: a very tiny converter from DNA string to event stream for MinION sequencer

- author: Bansho Masutani<banmasutani@gmail.com>
- language: Rust

## Short summary

This is a tiny program to convert input DNA string into nanopore-like event stream using 6-mer HMM provided by ONT.


## How to use

0. Clone this repository to local envirnment
```
hg clone https://ban-m@bitbucket.org/ban-m/squiggler
```
1. Add the directly you clone this repository to your Cargo.toml 
```
[dependencies]
squiggler = {path = "[the location]"}
```
2. Import this crate in .rs file and use it.
```
extern crate squiggler;
use squiggler::Squiggler;
use std::path::Path;
fn main(){
   let fasta = "> this is header \n AAAATTGCGCGCGAAAGCGTGGCGTTTGTG";
   let model = Squiggler::new(&Path::new("./path/to/model/file.model")).unwrap();
   let result = model.get_signal_from_fasta(fasta);
   println!("{:?}",result);
}
```

## Requirement

- Rust language with version 1.2 or later.
- Model file
    - You can get 6-mer HMM model from ONT github repository.

## Restrictions

Squiggler do almost nothing than just convert each 6-mer into specific real value.
This might be improved by other ML approach including NN.

In principle, squiggler can be used as a converter for R9.5 flowcell. However, currently there seems
to be no model for R9.5 and thus this crate only serves for R9.2/4 flowcell.
