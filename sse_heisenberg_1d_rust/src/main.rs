/**************************************************
 *  SSE for 1D S=1/2 Heisenberg model
 *  Author: Yi-Ming Ding
 *  Email: dingyiming@westlake.edu.cn
 *  Updated: Nov 28, 2024
 *************************************************/
mod hm;
mod aux;
use std::{env, time::Instant, io::Write, fs::OpenOptions};
use crate::aux::report_time;

fn main() {
    let start_time = Instant::now();
    let para_seed: u32 = aux::get_time_as_seed();

    // ::::::::::::::::::::::::::::::::::::
    //  Collect params from shell
    // ::::::::::::::::::::::::::::::::::::
    let args: Vec<String> = env::args().collect();
    let para_l: usize = args[1].parse().unwrap();
    let para_beta: f64 = args[2].parse().unwrap();
    let num_thm: usize = args[3].parse().unwrap();
    let num_stat: usize = args[4].parse().unwrap();
    let num_bins: usize = args[5].parse().unwrap();

    // ::::::::::::::::::::::::::::::::::::
    //  Report the environment
    // ::::::::::::::::::::::::::::::::::::
    println!("■ SSE for S=1/2 Heisenberg model, L = {para_l}, beta = {para_beta}");
    println!("■ num_thm = {num_thm}, num_stat = {num_stat}, num_bins = {num_bins}");

    // ::::::::::::::::::::::::::::::::::::
    //  Preparing for writing the results
    // ::::::::::::::::::::::::::::::::::::
    let mut file_energy = OpenOptions::new().write(true).append(true)
        .create(true).open("./data/energy.dat").unwrap();
    let mut file_zz = OpenOptions::new().write(true).append(true)
        .create(true).open("./data/zz_correlation.dat").unwrap();
    let mut file_heat_cap = OpenOptions::new().write(true).append(true)
        .create(true).open("./data/heat_capacity.dat").unwrap();

    // ::::::::::::::::::::::::::::::::::::
    //  Monte Carlo simulations
    // ::::::::::::::::::::::::::::::::::::
    let mut model = hm::HeisenbergModel::new(para_l, para_beta, para_seed);
    model.initialization();

    println!("\t---> Thermalization...");
    for _ in 0..num_thm {
        model.mc_update();
        model.adjust_m();
    }
    println!("\t---> Maximum cut-off = {}", model.m);

    println!("\t---> Measuring...");
    for _ in 0..num_bins {
        model.init_measure();
        for __ in 0..num_stat {
            model.mc_update();
            model.measure();
        }
        model.statisticize(num_stat);

        // ------------------------------------
        //  Saving the data
        // ------------------------------------
        file_energy.write_all(format!("{:<16.10}\n", model.energy).as_bytes()).unwrap();
        file_zz.write_all(format!("{}\n",
                                  model.zz_correlation.iter().map(|x| format!("{:<16.10}", x)).collect::<Vec<_>>().join("\t")
        ).as_bytes()).unwrap();
        file_heat_cap.write_all(format!("{:<16.10}\n", model.heat_capacity).as_bytes()).unwrap();
    }

    // ::::::::::::::::::::::::::::::::::::
    //  Report the runtime
    // ::::::::::::::::::::::::::::::::::::
    for _ in 0..50 { print!("-"); } print!("\n");
    report_time(start_time);
}