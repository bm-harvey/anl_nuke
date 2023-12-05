use std::{fs::File, io::BufRead, io::BufReader, path::Path};

use faust::phys_vec::PhysVec;

#[derive(Default)]
struct DetectorGeometry {
    normal_vector: PhysVec,
    corners: Vec<PhysVec>,
}

pub struct FaustFilter {
    geometries: Vec<DetectorGeometry>,
}

impl FaustFilter {
    pub const NUM_DETECTORS: usize = 68;
    pub const DETECTOR_AREA_CM2: f64 = 4.;
    pub const DETECTOR_AREA_TOLERANCE: f64 = 0.005;

    pub fn new(path_to_input: &Path) -> Self {
        let input = File::open(&path_to_input).unwrap();
        let mut buff = BufReader::new(input);

        let num_corners = 4;
        let num_detectors = FaustFilter::NUM_DETECTORS;

        let mut line = String::new();
        buff.read_line(&mut line).unwrap();

        let mut geometries = Vec::new();

        for _det_idx in 0..num_detectors {
            let mut geom = DetectorGeometry::default();
            for _corner_idx in 0..num_corners {
                line.clear();

                buff.read_line(&mut line).unwrap();
                let words = line.split_whitespace().collect::<Vec<_>>();

                let x = words[3].parse::<f64>().unwrap();
                let y = words[4].parse::<f64>().unwrap();
                let z = words[2].parse::<f64>().unwrap();

                let vec_to_corner = PhysVec::new([x, y, z]);

                geom.corners.push(vec_to_corner);
            }

            // the following takes advantage of the fact that the faust detectors are squares
            // with normal vector pointing at the target (origin)
            let mut norm = PhysVec::default();
            for idx in 0..num_corners {
                norm += &geom.corners[idx];
            }

            norm.scale(0.25);

            geom.normal_vector = norm;

            geometries.push(geom);
        }

        FaustFilter { geometries }
    }

    pub fn detector_hit(&self, trajectory: &PhysVec) -> Option<usize> {
        for detector_idx in 0..FaustFilter::NUM_DETECTORS {
            let geom = &self.geometries[detector_idx];

            let mut normalized_trajectory = trajectory.clone();

            let scaling_factor =
                geom.normal_vector.mag() / trajectory.inner_angle_rad(&geom.normal_vector).cos();


            normalized_trajectory.scale(scaling_factor / trajectory.mag());

            let mut area = 0.;

            for corner_combo in 0..4 {
                area += FaustFilter::triangle_area(
                    &normalized_trajectory,
                    &geom.corners[(corner_combo + 1) % 4],
                    &geom.corners[corner_combo],
                );
            }
            //dbg!(&area);
            if area <= FaustFilter::DETECTOR_AREA_CM2 + FaustFilter::DETECTOR_AREA_TOLERANCE {
                return Some(detector_idx);
            }
        }

        None
    }

    /// calculate the area of a triangle described by 3 points with Heron's formula
    pub fn triangle_area(v1: &PhysVec, v2: &PhysVec, v3: &PhysVec) -> f64 {
        let s1 = (v2 - v1).mag();
        let s2 = (v3 - v2).mag();
        let s3 = (v1 - v3).mag();

        let p = 0.5 * (s1 + s2 + s3);

        (p * (p - s1) * (p - s2) * (p - s3)).sqrt()
    }
}
