use rand::Rng;

fn main() {
    println!("trying to create file");

    let file = roost::file::create("my_file.root".into());

    let tree = roost::tree::new_tree("T".into(), "T".into());

    let hist = roost::hist::new_h1d("hist".into(), ";arb. units;".into(), 20, -5., 6.);
    let hist_2d = roost::hist::new_h2d(
        "hist_2d".into(),
        ";arb. units;arb. units".into(),
        1_000,
        -5.0,
        5.0,
        1_000,
        -5.0,
        5.0,
    );

    let mut branch_x = tree.make_branch("x".into());
    let mut branch_y = tree.make_branch("y".into());

    for _idx in 0..100_000 {
        hist.fill(rand::thread_rng().gen_range(-0.3..3.0));
        hist_2d.fill(
            rand::thread_rng().sample(rand_distr::StandardNormal),
            rand::thread_rng().sample(rand_distr::StandardNormal),
        );

        branch_x
            .pin_mut()
            .set_value(rand::thread_rng().sample(rand_distr::StandardNormal));
        branch_y
            .pin_mut()
            .set_value(rand::thread_rng().sample(rand_distr::StandardNormal));

        tree.fill();
    }

    hist.write();
    hist_2d.write();
    tree.write();
    file.close();
}
