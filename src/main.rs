mod bam_io;
mod cli;
mod pipeline;
mod reference;
mod stats;
mod tags;

use anyhow::Result;
use clap::Parser;

fn main() -> Result<()> {
    pipeline::run(cli::Args::parse())
}
