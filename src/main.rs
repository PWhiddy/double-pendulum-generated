
use indicatif::ProgressBar;
use image::{ImageBuffer, Rgba};
use std::fs::create_dir_all;
use std::path::Path;

struct PendulumState {
    theta1: f64, // Angle of the first pendulum
    theta2: f64, // Angle of the second pendulum
    omega1: f64, // Angular velocity of the first pendulum
    omega2: f64, // Angular velocity of the second pendulum
}

impl PendulumState {
    fn new(theta1: f64, theta2: f64, omega1: f64, omega2: f64) -> Self {
        PendulumState { theta1, theta2, omega1, omega2 }
    }

    fn update(&mut self, dt: f64) {
        let g = 9.81; // acceleration due to gravity, in m/s^2
        let l1 = 1.0; // length of the first pendulum rod, in meters
        let l2 = 1.0; // length of the second pendulum rod, in meters
        let m1 = 1.0; // mass of the first pendulum bob, in kg
        let m2 = 1.0; // mass of the second pendulum bob, in kg

        // Helper function to calculate the derivatives of the state variables
        let derivatives = |theta1: f64, theta2: f64, omega1: f64, omega2: f64| -> (f64, f64, f64, f64) {
            let delta_theta = theta2 - theta1;
            let den1 = (m1 + m2) * l1 - m2 * l1 * delta_theta.cos().powi(2);
            let den2 = (l2 / l1) * den1;

            let dtheta1 = omega1;
            let dtheta2 = omega2;

            let domega1 = (m2 * g * theta2.sin() * delta_theta.cos()
                           - m2 * delta_theta.sin() * (l1 * omega1.powi(2) * delta_theta.cos() + l2 * omega2.powi(2))
                           - (m1 + m2) * g * theta1.sin()) / den1;

            let domega2 = ((m1 + m2) * (l1 * omega1.powi(2) * delta_theta.sin() - g * theta2.sin() + g * theta1.sin() * delta_theta.cos())
                           + m2 * l2 * omega2.powi(2) * delta_theta.sin() * delta_theta.cos()) / den2;

            (dtheta1, dtheta2, domega1, domega2)
        };

        // Runge-Kutta 4th order method to update the state
        let (k1_theta1, k1_theta2, k1_omega1, k1_omega2) = derivatives(self.theta1, self.theta2, self.omega1, self.omega2);
        let (k2_theta1, k2_theta2, k2_omega1, k2_omega2) = derivatives(self.theta1 + dt / 2.0 * k1_theta1, self.theta2 + dt / 2.0 * k1_theta2, self.omega1 + dt / 2.0 * k1_omega1, self.omega2 + dt / 2.0 * k1_omega2);
        let (k3_theta1, k3_theta2, k3_omega1, k3_omega2) = derivatives(self.theta1 + dt / 2.0 * k2_theta1, self.theta2 + dt / 2.0 * k2_theta2, self.omega1 + dt / 2.0 * k2_omega1, self.omega2 + dt / 2.0 * k2_omega2);
        let (k4_theta1, k4_theta2, k4_omega1, k4_omega2) = derivatives(self.theta1 + dt * k3_theta1, self.theta2 + dt * k3_theta2, self.omega1 + dt * k3_omega1, self.omega2 + dt * k3_omega2);

        self.theta1 += dt / 6.0 * (k1_theta1 + 2.0 * k2_theta1 + 2.0 * k3_theta1 + k4_theta1);
        self.theta2 += dt / 6.0 * (k1_theta2 + 2.0 * k2_theta2 + 2.0 * k3_theta2 + k4_theta2);
        self.omega1 += dt / 6.0 * (k1_omega1 + 2.0 * k2_omega1 + 2.0 * k3_omega1 + k4_omega1);
        self.omega2 += dt / 6.0 * (k1_omega2 + 2.0 * k2_omega2 + 2.0 * k3_omega2 + k4_omega2);
    }
}

fn main() {
    let n_pendulums = 10000; // Number of pendulums to simulate
    let dt = 0.01; // Time step for the simulation
    let total_time_steps = 3000; // Total number of simulation steps

    let mut pendulums = Vec::new();
    let offset = 0.000005;
    // Initialize pendulums with slightly differing starting conditions
    for i in 0..n_pendulums {
        pendulums.push(PendulumState::new(1.0 + i as f64 * offset, 1.0 + i as f64 * offset, 0.0, 0.0));
    }

    // Setup directory for frame output
    let frames_dir = Path::new("frames");
    create_dir_all(&frames_dir).expect("Failed to create frames directory");

    // Setup progress bar
    let progress_bar = ProgressBar::new(total_time_steps);

    // Simulation and rendering loop
    for t in 0..total_time_steps {
        for pendulum in pendulums.iter_mut(){
            pendulum.update(dt);
        }
        render_pendulums_to_image(&pendulums, t);
        progress_bar.inc(1);
    }

    println!("Simulation complete. Use ffmpeg to compile frames into a video.");
}

fn render_pendulums_to_image(pendulums: &Vec<PendulumState>, frame_number: u64) {
    let width = 1024u32;
    let height = 1024u32;
    let center_x = width as f64 / 2.0;
    let center_y = height as f64 / 2.0;
    let scale = 200.0; // Scale for drawing
    let mut img = ImageBuffer::new(width, height);

    for pendulum in pendulums {
        // Convert pendulum angles to end points
        let x1 = center_x + pendulum.theta1.sin() * scale;
        let y1 = center_y - pendulum.theta1.cos() * scale;
        let x2 = x1 + pendulum.theta2.sin() * scale;
        let y2 = y1 - pendulum.theta2.cos() * scale;

        // Draw the rods and bobs for each pendulum
        draw_line(&mut img, center_x as i32, center_y as i32, x1 as i32, y1 as i32, Rgba([255, 0, 0, 255]));
        draw_circle(&mut img, x1 as i32, y1 as i32, 10, Rgba([0, 0, 255, 255]));
        draw_line(&mut img, x1 as i32, y1 as i32, x2 as i32, y2 as i32, Rgba([255, 0, 0, 255]));
        draw_circle(&mut img, x2 as i32, y2 as i32, 10, Rgba([0, 255, 0, 255]));
    }

    // Save the layered image for the current frame
    img.save(format!("frames/frame_{:04}.png", frame_number)).unwrap();
}

// Bresenham's line algorithm for drawing lines
fn draw_line(img: &mut ImageBuffer<Rgba<u8>, Vec<u8>>, x0: i32, y0: i32, x1: i32, y1: i32, color: Rgba<u8>) {
    let mut x = x0;
    let mut y = y0;
    let dx = (x1 - x0).abs();
    let dy = -(y1 - y0).abs();
    let sx = if x0 < x1 { 1 } else { -1 };
    let sy = if y0 < y1 { 1 } else { -1 };
    let mut err = dx + dy; // error value e_xy

    loop {
        img.put_pixel(x as u32, y as u32, color);
        if x == x1 && y == y1 { break; }
        let e2 = 2 * err;
        if e2 >= dy { // e_xy+e_x > 0
            err += dy;
            x += sx;
        }
        if e2 <= dx { // e_xy+e_y < 0
            err += dx;
            y += sy;
        }
    }
}

// Midpoint circle algorithm for drawing circles
fn draw_circle(img: &mut ImageBuffer<Rgba<u8>, Vec<u8>>, cx: i32, cy: i32, radius: i32, color: Rgba<u8>) {
    let mut x = -radius;
    let mut y = 0;
    let mut err = 2 - 2 * radius; // second parameter of "decision criterion"
    while x < 0 {
        img.put_pixel((cx - x) as u32, (cy + y) as u32, color); // I. Quadrant
        img.put_pixel((cx - y) as u32, (cy - x) as u32, color); // II. Quadrant
        img.put_pixel((cx + x) as u32, (cy - y) as u32, color); // III. Quadrant
        img.put_pixel((cx + y) as u32, (cy + x) as u32, color); // IV. Quadrant
        let r = err;
        if r <= y {
            err += y * 2 + 1;
            y += 1; // Increment y after using it in the calculation
        }
        if r > x || err > y {
            err += x * 2 + 1;
            x += 1; // Increment x after using it in the calculation
        }
    }
}

