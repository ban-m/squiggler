#[cfg(test)]
mod test{
use super::*;
use std::path::Path;

#[test]
fn test_of_test(){
    assert!(true);
}
const MODEL_PATH:&str = "./kmer_models/r9.4_180mv_450bps_6mer/template_median68pA.model";

#[test]
fn squiggler_init2(){
    match squiggler::Squiggler::new(&Path::new("not/existing/file")){
        Ok(_) => assert!(false),
        Err(_) => assert!(true),
    }
}

#[test]
fn squiggler_init3(){
    // not valid data type
    match squiggler::Squiggler::new(&Path::new("./src/testdata/invalid.model")){
        Ok(_) => assert!(false),
        Err(_) => assert!(true),
    }
}

#[test]
fn squiggler_init4(){
    //model that doen't contain enough data
    match squiggler::Squiggler::new(&Path::new("./src/testdata/notenoughdata.model")){
        Ok(_) => assert!(false),
        Err(_) => assert!(true),
    }
}
#[test]
fn squiggler_init5(){
    let model = squiggler::Squiggler::new(&Path::new(MODEL_PATH)).unwrap();
    assert_eq!(model.get_size(),4*4*4*4*4*4);
}

#[test]
fn squiggler_getsignal(){
    let model = squiggler::Squiggler::new(&Path::new(MODEL_PATH)).unwrap();
    assert_eq!(model.get_signal_from_fasta(">test of fasta\nAAAAAA").len(),1);
}

#[test]
fn squiggler_getsignal2(){
    let model = squiggler::Squiggler::new(&Path::new(MODEL_PATH)).unwrap();
    assert_eq!(model.get_signal_from_fasta(">test of fasta\nAAAAAAAAAAAA").len(),7);
}


#[test]
fn signal_mode_init(){
    let sm = squiggler::SignalModel::new(0.,0.,0.,0.);
    assert!(true);
}
#[test]
fn signal_mode_levelmeancheck(){
    let sm = squiggler::SignalModel::new(10.,0.,0.,0.);
    assert_eq!(sm.generate().2,10.);
}
#[test]
fn signal_mode_levelstddevcheck(){
    let sm = squiggler::SignalModel::new(1.,0.,0.,0.);
    debug_assert!(sm.generate().2-1. < 0.1,"{:?}",sm.generate());
}

#[test]
fn signal_mode_noisemeancheck(){
    let sm = squiggler::SignalModel::new(0.,1.,0.,0.);
    assert!(true);
}

#[test]
fn signal_mode_noisesdtvcheck(){
    /// write as soon as possible
    assert!(true);
}


#[test]
fn final_check(){
    use std::io::prelude::*;
    use std::io::BufRead;
    use std::io::BufReader;
    let model = squiggler::Squiggler::new(&Path::new(MODEL_PATH)).unwrap();
    let lines = BufReader::new((std::fs::File::open(MODEL_PATH).unwrap())).lines().
        filter_map(|e| e.ok()).skip(1);//header
    for line in lines{
        let content :Vec<_>= line.split('\t').collect();
        let mer = content[0].clone();
        let mean:f32 = content[1].parse().unwrap();
        debug_assert!(model.get_signal_from_fasta(&mer)[0].2== mean,
                      "{}{:?}vs {}",mer,content,model.get_signal_from_fasta(&mer)[0].2);
    }
}
}
