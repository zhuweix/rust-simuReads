use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use bio::io::fasta;
use bio::alphabets::dna::revcomp;
use clap::Parser;

#[derive(Parser)]
#[command(name="SimuReads", version, 
        about="Simulate reads from a fasta file", long_about = None,
        author = env!("CARGO_PKG_AUTHORS"))]
struct Cli {
    /// Input reference file in fasta format
    #[arg(short='i', long="input")]
    input: PathBuf,
    /// Output prefix for the generated reads
    #[arg(short='o', long="output")]
    output: String,
    /// Read length
    #[arg(short='l', long="readlength", default_value_t=50)]
    read_length: usize,
    /// Insert size
    #[arg(short='s', long="insertsize", default_value_t=300)]
    insert_size: usize,
}


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    if cli.read_length * 2 > cli.insert_size {
        return Err("Read length should be less than half of insert size".into());
    }

    let reader = fasta::Reader::from_file(cli.input)?;
    let mut writer1 = BufWriter::new(File::create(format!("{}.r1.fasta", cli.output))?);
    let mut writer2 = BufWriter::new(File::create(format!("{}.r2.fasta", cli.output))?);

    for record in reader.records() {
        let record = record?;
        let seq = record.seq();
        let seq_len = seq.len();
        let mut start = 0;
        while start + cli.insert_size <= seq_len {
            let end = start + cli.insert_size;
            let read = &seq[start..end];
            let read_rc = revcomp(read);
            let read1 = &read[..cli.read_length];
            let read2 = &read_rc[..cli.read_length];
            writeln!(writer1, ">{}-{}:W", record.id(), start)?;
            writeln!(writer1, "{}", std::str::from_utf8(read1)?)?;
            writeln!(writer2, ">{}-{}:C", record.id(), end)?;
            writeln!(writer2, "{}", std::str::from_utf8(read2)?)?;
            start += 1;
        }
    }

    Ok(())
}