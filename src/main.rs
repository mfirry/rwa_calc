use rwa_calc::rwa;

fn main() {
    let correlation_multiplier = 1.0;
    let sme_factor = 1.0;
    let key_methodical_approach = "IRA";
    let risk_weight_transaction = 0.8632362;

    let x = rwa(
        0.0242119,
        0.3585257,
        1.0,
        140472.9,
        81427000.0,
        correlation_multiplier,
        sme_factor,
        key_methodical_approach.to_string(),
        risk_weight_transaction,
        "CL".to_string(),
    );

    println!("{}", x);
}
