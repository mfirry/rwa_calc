extern crate statrs;

use crate::statrs::distribution::ContinuousCDF;
use statrs::distribution::Normal;
use std::f64::consts::*;

pub fn rwa(pd: f64, lgd: f64, residual_maturity: u8, exposure: f64) -> f64 {
    let n = Normal::new(0.0, 1.0).unwrap();
    let corr_cl_cs1 = 0.12 * (1.0 - E.powf(pd * -50.0));
    let corr_cl_cs2 =
        (1.0 - E.powf(-50.0) + 0.24 * (1.0 - (1.0 - E.powf(-50.0 * pd)))) / (1.0 - E.powf(-50.0));
    let corr_cl_cs = corr_cl_cs1 / corr_cl_cs2;
    let factor_ = f64::powf(1.0 - corr_cl_cs, -0.5);
    let w = factor_ * n.inverse_cdf(pd)
        + f64::powf(corr_cl_cs / (1.0 - corr_cl_cs), 0.5) * n.inverse_cdf(0.999);
    let rw1 = lgd * n.cdf(w) - pd * lgd;
    let b = f64::powf(0.11852 - 0.05475 * f64::log(pd, E), 2.0);
    let rw2 = (1.0 + ((residual_maturity as f64) - 2.5) * b) / (1.0 - 1.5 * b) * 12.5 * 1.06;
    return rw1 * rw2 * exposure;
}
