fn main() {
    cxx_build::bridge("src/lib.rs")
        .file("cpp/src/my_bridge.cpp")
        .include("D:\\root6\\root_install\\include")
        .include("cpp\\src\\")
        .include("roost\\cpp")
        .flag("-std:c++17")
        //.flag("-nologo")
        .flag("-Zc:__cplusplus")
        //.flag("-MD")
        //.flag("-GR")
        ////.flag("-D_WIN32")
        .flag("-O2")
        .compile("roost-cpp");

    println!("cargo:rustc-link-search=D:\\root6\\root_install\\lib");
    println!("cargo:rerun-if-changed=src/lib.rs");
    println!("cargo:rerun-if-changed=cpp/src/*");

    println!("cargo:rustc-link-lib=libCore");
    println!("cargo:rustc-link-lib=libImt");
    println!("cargo:rustc-link-lib=libRIO");
    println!("cargo:rustc-link-lib=libNet");
    println!("cargo:rustc-link-lib=libHist");
    println!("cargo:rustc-link-lib=libGraf");
    println!("cargo:rustc-link-lib=libGraf3d");
    println!("cargo:rustc-link-lib=libROOTVecOps");
    println!("cargo:rustc-link-lib=libTree");
    println!("cargo:rustc-link-lib=libTreePlayer");
    println!("cargo:rustc-link-lib=libRint");
    println!("cargo:rustc-link-lib=libPostscript");
    println!("cargo:rustc-link-lib=libMatrix");
    println!("cargo:rustc-link-lib=libPhysics");
    println!("cargo:rustc-link-lib=libMathCore");
    println!("cargo:rustc-link-lib=libThread");
}
