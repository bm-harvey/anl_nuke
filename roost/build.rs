fn main() {
    cxx_build::bridge("src/lib.rs")
        .file("cpp/src/my_bridge.cpp")
        .include("/opt/new/include/root")
        //.include("C:\\Users\\bmhar\\root_6\\root_install\\include")
        //.include("D:\\root6\\root_install\\include")
        .include("roost\\cpp")
        .include("cpp\\src\\")
        //.flag("-std:c++17")
        .flag("-std=c++14")
        //.flag("-nologo")
        //.flag("-Zc:__cplusplus")
        //.flag("-MD")
        //.flag("-GR")
        ////.flag("-D_WIN32")
        .flag("-O2")
        .compile("roost-cpp");

    println!("cargo:rustc-link-search=/opt/new/lib");

    println!("cargo:rerun-if-changed=src/lib.rs");
    println!("cargo:rerun-if-changed=cpp/src/*");

    println!("cargo:rustc-link-lib=Core");
    println!("cargo:rustc-link-lib=Imt");
    println!("cargo:rustc-link-lib=RIO");
    println!("cargo:rustc-link-lib=Net");
    println!("cargo:rustc-link-lib=Hist");
    println!("cargo:rustc-link-lib=Graf");
    println!("cargo:rustc-link-lib=Graf3d");
    println!("cargo:rustc-link-lib=ROOTVecOps");
    println!("cargo:rustc-link-lib=Tree");
    println!("cargo:rustc-link-lib=TreePlayer");
    println!("cargo:rustc-link-lib=Rint");
    println!("cargo:rustc-link-lib=Postscript");
    println!("cargo:rustc-link-lib=Matrix");
    println!("cargo:rustc-link-lib=Physics");
    println!("cargo:rustc-link-lib=MathCore");
    println!("cargo:rustc-link-lib=Thread");
}
