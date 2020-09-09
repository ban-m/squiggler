use rand;
use rand::distributions::IndependentSample;
use rand::distributions::Normal;
use rand::Rng;
use std;
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::vec::Vec;
//const KMER_SIZE : usize = 6;
const LENGTH: u32 = 3;
const START_TIME: u32 = 0;
#[derive(Debug)]
pub struct SignalModel {
    level_mean: f32,
    level_stdv: f32,
    sd_mean: f32,
    sd_stdv: f32,
}

impl SignalModel {
    pub fn new(level_mean: f32, level_stdv: f32, sd_mean: f32, sd_stdv: f32) -> SignalModel {
        SignalModel {
            level_mean,
            level_stdv,
            sd_mean,
            sd_stdv,
        }
    }
    ///output signal
    ///[inverse gaussian distribution][https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution#cite_note-2]
    pub fn generate(&self) -> (u32, u32, f32, f32) {
        let _normal = Normal::new(self.level_mean as f64, self.level_stdv as f64);
        let _v = _normal.ind_sample(&mut rand::thread_rng()) as f32;
        (START_TIME, LENGTH, self.level_mean, self.sd_mean)
    }
    pub fn getmean(&self) -> f32 {
        self.level_mean
    }
}

#[derive(Debug)]
pub struct Squiggler {
    kmersize: usize,
    pub model: HashMap<String, SignalModel>,
}

#[derive(Debug)]
pub enum SquigglerError {
    NotValidModel,
    NotEnoughData,
    IoError,
    ParseError,
}

impl Squiggler {
    pub fn new(path: &std::path::Path) -> Result<Squiggler, SquigglerError> {
        use self::SquigglerError::{IoError, NotEnoughData, NotValidModel, ParseError};

        let mut model = HashMap::with_capacity(4usize.pow(6));
        let kmersize = 6;
        let file = File::open(path).map_err(|_| IoError)?;
        let bufreader = BufReader::new(file);
        for line in bufreader.lines() {
            let line = line.map_err(|_| IoError)?;
            if line.starts_with("kmer") {
                // its a header line
                continue;
            } else {
                let contents: Vec<_> = line.split('\t').collect();
                let kmer = contents[0].to_string();
                let mean: f32 = contents[1].parse().map_err(|_| ParseError)?;
                let stdv: f32 = contents[2].parse().map_err(|_| ParseError)?;
                let nmean: f32 = contents[3].parse().map_err(|_| ParseError)?;
                let nstdv: f32 = contents[4].parse().map_err(|_| ParseError)?;
                let kmermodel = SignalModel::new(mean, stdv, nmean, nstdv);
                if model.insert(kmer, kmermodel).is_some() {
                    return Err(NotValidModel);
                }
            }
        }
        if model.keys().count() < (4usize.pow(6)) {
            println!("not good:{:?}", path);
            Err(NotEnoughData)
        } else {
            Ok(Squiggler { kmersize, model })
        }
    }
    pub fn get_signal_from_path(
        &self,
        fasta_path: &std::path::Path,
    ) -> std::io::Result<Vec<(u32, u32, f32, f32)>> {
        use std::io::prelude::*;
        let mut file = File::open(fasta_path)?;
        let mut fasta = String::new();
        file.read_to_string(&mut fasta)?;
        Ok(self.get_signal_from_fasta(&fasta))
    }
    pub fn get_signal_from_fasta(&self, fasta: &str) -> Vec<(u32, u32, f32, f32)> {
        let mut seq = String::new();
        for line in fasta.to_string().split('\n') {
            if line.starts_with('>') || line.starts_with(';') {
                continue;
            } else {
                seq.push_str(line);
            }
        }
        let range = seq.len() - self.kmersize + 1;
        let mut res = vec![];
        for index in 0..range {
            res.push(
                self.model
                    .get(&seq[index..index + self.kmersize])
                    .unwrap()
                    .generate(),
            );
        }
        res
    }
    pub fn get_all_kmer(&self) -> Vec<String> {
        self.model.keys().map(|key| key.to_string()).collect()
    }

    /// convert the file(fasta) to squiggle and
    /// output to the standard output
    pub fn convert(&self, fasta_path: &str) {
        let fasta_path = std::path::Path::new(fasta_path);
        let signals = match self.get_signal_from_path(fasta_path) {
            Ok(ok) => ok,
            Err(_) => {
                println!("invalid file");
                vec![]
            }
        };
        println!("start,length,mean,stdv");
        let mut start = START_TIME;
        for signal in signals {
            println!("{},{},{},{}", start, signal.1, signal.2, signal.3);
            start += signal.1;
        }
    }
    pub fn get_size(&self) -> usize {
        self.model.len()
    }
}

/// clipping given events length for a given probability
pub fn clipping<T>(events: &mut Vec<T>, prob: f32) {
    let mut rng = rand::thread_rng();
    events.retain(|_| rng.gen_range(0., 1.) < prob);
}

/// remove duplicate in given range
pub fn dedup(events: &[f32], range: f32) -> Vec<f32> {
    events
        .iter()
        .fold((vec![], events[0], 1.), |(mut acc, prev, num), &x| {
            if (prev - x).abs() < range {
                // dedup
                (acc, (prev * num + x) / (num + 1.), num + 1.)
            } else {
                //acc.push(x); wrong implementaion
                acc.push(prev);
                (acc, x, 1.)
            }
        })
        .0
}

/// scrappie column
#[derive(Debug, Clone, Copy)]
pub struct Signal {
    /// base
    pub base: char,
    /// current
    pub current: f32,
    /// sd
    pub sd: f32,
    /// dwell
    pub dwell: f32,
}

impl Signal {
    /// create new Signal
    pub fn new(base: char, current: f32, sd: f32, dwell: f32) -> Signal {
        Signal {
            base,
            current,
            sd,
            dwell,
        }
    }
}

/// convert scrappie prediction into vector
pub fn parse(path: &std::path::Path) -> std::io::Result<Vec<Signal>> {
    let res = BufReader::new(File::open(path)?)
        .lines()
        .filter_map(|e| e.ok())
        .filter(|e| !e.starts_with('#'))
        .filter_map(|line| {
            let contents: Vec<_> = line.split('\t').collect();
            let base = match contents[1].chars().next() {
                Some(c) => c,
                None => return None,
            };
            let current: f32 = match contents[2].parse() {
                Ok(res) => res,
                Err(_) => return None,
            };
            let sd: f32 = match contents[3].parse() {
                Ok(res) => res,
                Err(_) => return None,
            };
            let dwell: f32 = match contents[4].parse() {
                Ok(res) => res,
                Err(_) => return None,
            };
            Some(Signal::new(base, current, sd, dwell))
        })
        .collect();
    Ok(res)
}
