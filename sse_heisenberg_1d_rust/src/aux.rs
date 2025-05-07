use std::time::{Instant, SystemTime, UNIX_EPOCH};

pub fn get_time_as_seed() -> u32 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("Time went backwards")
        .as_secs() as u32
}

pub fn report_time(start_time: Instant) {
    let elapsed = start_time.elapsed().as_secs() as u32;
    let hours = elapsed / 3600;
    let minutes = (elapsed % 3600) / 60;
    let seconds = elapsed % 60;
    println!("■ Runtime: {}h-{}min-{}s", hours, minutes, seconds);
}