extern crate statrs;

use crate::statrs::distribution::ContinuousCDF;
use statrs::distribution::Normal;
use std::f64::consts::*;

// #[derive(Serialize)]
// enum Brw {
//     CS,
//     CL,
//     RE,
//     RM,
//     RR,
//     XX,
// }

pub fn rwa(
    pd: f64,
    lgd: f64,
    residual_maturity: f64,
    exposure: f64,
    turnover: f64,
    correlation_multiplier: f64,
    sme_factor: f64,
    key_methodical_approach: String,
    risk_weight_transaction: f64,
    brw: String,
) -> f64 {
    let n = Normal::new(0.0, 1.0).unwrap();

    let corr_cl_cs1 = 0.12 * (1.0 - E.powf(pd * -50.0));
    let corr_cl_cs2 =
        (1.0 - E.powf(-50.0) + 0.24 * (1.0 - (1.0 - E.powf(-50.0 * pd)))) / (1.0 - E.powf(-50.0));
    let corr_cl_cs = corr_cl_cs1 / corr_cl_cs2;

    let corr_re = (0.03 * (1.0 - E.powf(-35.0) * pd)) / (1.0 - E.powf(-35.0))
        + 0.16 * (1.0 - (1.0 - E.powf(-35.0) * pd)) / (1.0 - E.powf(-35.0));
    let corr_rm = 0.15;
    let corr_rr = 0.04;

    let mut k2 = 0.0;

    if key_methodical_approach == "STC" {
        return exposure * risk_weight_transaction;
    } else {
        match brw.as_str() {
            "CL" => {}
            "XX" => {
                k2 = corr_cl_cs;
            }
            "CS" => {
                let turnover_mln = turnover / 1000000.0;
                if turnover_mln >= 5.0 && turnover_mln <= 50.0 {
                    k2 = corr_cl_cs
                        - 0.04
                            * (1.0 - ((f64::min(f64::max(5.0, turnover_mln), 50.0) - 5.0) / 45.0));
                }
            }
            "RE" => {
                k2 = corr_re;
            }
            "RM" => {
                k2 = corr_rm;
            }
            "RR" => {
                k2 = corr_rr;
            }
            _ => {}
        }

        let tmp = f64::sqrt(1.0 - k2) * n.inverse_cdf(pd)
            + f64::sqrt(k2 / (1.0 - k2)) * n.inverse_cdf(0.999);
        let rw1 = lgd * n.inverse_cdf(tmp) - pd * lgd;
        let b = f64::powf(0.11852 - 0.05475 * f64::log(pd, E), 2.0);

        if brw.eq("CS") || brw.eq("CL") || brw.eq("XX") {
            let rw2 = (1.0 + (residual_maturity - 2.5) * b) / (1.0 - 1.5 * b) * 12.5 * 1.06;
            return rw1 * rw2 * exposure * correlation_multiplier * sme_factor;
        } else {
            let rw2 = 1.0 * 12.5 * 1.06;
            return rw1 * rw2 * exposure * correlation_multiplier * sme_factor;
        }
    }
}
