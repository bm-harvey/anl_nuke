#[cxx::bridge]
pub mod file {
    #[namespace = "file"]
    unsafe extern "C++" {
        include!("roost/cpp/src/my_bridge.hpp");

        type RtFile;
        fn open(name: String) -> SharedPtr<RtFile>;
        fn create(name: String) -> SharedPtr<RtFile>;
        fn write(&self) -> ();
        fn close(&self) -> ();
        fn cd(&self) -> ();
    }
}
#[cxx::bridge]
pub mod hist {
    #[namespace = "h1d"]
    unsafe extern "C++" {
        include!("roost/cpp/src/my_bridge.hpp");

        type RtH1D;
        fn fill(&self, value: f64);
        fn weighted_fill(&self, value: f64, weight: f64);
        fn write(&self) -> ();
        fn new_h1d(name: String, title: String, bins: i64, low: f64, high: f64)
            -> SharedPtr<RtH1D>;
    }

    #[namespace = "h2d"]
    unsafe extern "C++" {
        include!("roost/cpp/src/my_bridge.hpp");

        type RtH2D;
        fn fill(&self, value_x: f64, value_y: f64);
        fn weighted_fill(&self, value_x: f64, value_y: f64, weight: f64);
        fn new_h2d(
            name: String,

            title: String,
            bins_x: i64,
            low_x: f64,
            high_x: f64,
            bins_y: i64,
            low_y: f64,
            high_y: f64,
        ) -> SharedPtr<RtH2D>;
        fn write(&self) -> ();
    }
}

#[cxx::bridge]
pub mod tree {
    #[namespace = "tree"]
    unsafe extern "C++" {
        include!("roost/cpp/src/my_bridge.hpp");

        type RtTree;
        type RtBranch;
        type RtIntBranch;
        fn new_tree(name: String, title: String) -> SharedPtr<RtTree>;
        fn write(self: &RtTree) -> ();
        fn fill(self: &RtTree) -> ();
        fn make_branch(self: &RtTree, name: String) -> UniquePtr<RtBranch>;
        fn make_int_branch(self: &RtTree, name: String) -> UniquePtr<RtIntBranch>;

        fn set_value(self: Pin<&mut RtBranch>, value: f64);

        fn set_value(self: Pin<&mut RtIntBranch>, value: i64);
    }
}
