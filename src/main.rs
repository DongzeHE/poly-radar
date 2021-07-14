use bio::io::fasta;
use bio::pattern_matching::myers::Myers;
use std::fs::File;
use std::io::Write;

fn main() {
    let fasta_file_path = std::env::args().nth(1).expect("no pattern given");
let out_file_path = std::env::args().nth(2).expect("no path given");
    // read in file
    let fasta_file = std::path::Path::new(
        &fasta_file_path
    );
    let mut out_file = File::create(
        &out_file_path
    )
    .unwrap();
    let reader = fasta::Reader::from_file(fasta_file).unwrap();

    // set up probes
    let mut poly_a = Myers::<u64>::new(vec![65u8; 11]);
    let mut poly_t = Myers::<u64>::new(vec![84u8; 11]);

    // iterate over records
    for result in reader.records() {
        let record = result.expect("Error during fasta record parsing");

        // get matches
        let mut poly_occ: Vec<_> = poly_a
            .find_all(record.seq(), 1)
            .map(|(a, _, _)| a as u64)
            .collect();
        let mut poly_t_occ: Vec<_> = poly_t
            .find_all(record.seq(), 1)
            .map(|(a, _, _)| a as u64)
            .collect();
        poly_occ.append(&mut poly_t_occ);
        poly_occ.sort();
        poly_occ.dedup();
        // if no hit, next record
        if poly_occ.is_empty() {
            writeln!(out_file, "{}\t{}", record.id(), 0).unwrap();
            continue;
        }

        // if there are hits, process polyA then polyT
        if poly_occ.len() > 1 {
            // process the hit results in the format Vec<(start, end, edit distance)>
            // a long polyA will produce many hits, so remove continuous hits
            let mut keep = vec![true; poly_occ.len()];
            for idx in 1..poly_occ.len() {
                if poly_occ[idx] <= (poly_occ[idx - 1] + 11) {
                    keep[idx] = false;
                }
            }
            let mut i = 0;
            poly_occ.retain(|_| (keep[i], i += 1).0);
        }
        // println!("{:?}", type_of(poly_t_occ.first()));
        // writeln!(out_file, "{}\t{}\t{}", record.id(), poly_occ.len(), poly_a_occ.join("\t"))
        writeln!(out_file, "{}\t{}", record.id(), poly_occ.len())
        // writeln!(
        //     &mut out_file,
        //     "{}\t{}\t{}",
        //     record.id(),
        //     poly_occ.len(),
        //     poly_occ
        //         .into_iter()
        //         .map(|a| a.to_string())
        //         .collect::<Vec<String>>()
        //         .join("\t")
        // )
        .unwrap();
    }
}
