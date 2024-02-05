use super::bio_io::PbChem;
use super::predict_m6a::PredictOptions;
use burn::tensor::{Shape, Tensor};
// burn backend
/*

use burn::backend::candle::CandleDevice;
use burn::backend::Candle;
pub type BurnBackend = Candle;
pub type BurnDevice = CandleDevice;

use burn::backend::ndarray::NdArrayDevice;
use burn::backend::NdArray;
pub type BurnBackend = NdArray;
pub type BurnDevice = NdArrayDevice;
*/

use burn::backend::wgpu::WgpuDevice;
use burn::backend::Wgpu;
pub type BurnBackend = Wgpu;
pub type BurnDevice = WgpuDevice;

use super::predict_m6a::{LAYERS, WINDOW};

pub mod two_zero {
    include!(concat!(env!("OUT_DIR"), "/m6a_burn/two_zero.rs"));
}
pub mod two_two {
    include!(concat!(env!("OUT_DIR"), "/m6a_burn/two_two.rs"));
}
pub mod three_two {
    include!(concat!(env!("OUT_DIR"), "/m6a_burn/three_two.rs"));
}
pub mod revio {
    include!(concat!(env!("OUT_DIR"), "/m6a_burn/revio.rs"));
}

#[derive(Debug)]
pub struct BurnModels {
    pub two_zero: two_zero::Model<BurnBackend>,
    pub two_two: two_two::Model<BurnBackend>,
    pub three_two: three_two::Model<BurnBackend>,
    pub revio: revio::Model<BurnBackend>,
    pub device: BurnDevice,
}

impl BurnModels {
    pub fn new() -> Self {
        let two_zero = two_zero::Model::default();
        let two_two = two_two::Model::default();
        let three_two = three_two::Model::default();
        let revio = revio::Model::default();
        let device = BurnDevice::default();

        // log info about the device used
        log::info!("Using {:?} for Burn device.", device);
        Self {
            two_zero,
            two_two,
            three_two,
            revio,
            device,
        }
    }

    pub fn forward(&self, opts: &PredictOptions, windows: &[f32], count: usize) -> Vec<f32> {
        // predict
        if count == 0 {
            return vec![];
        }
        let shape = Shape::new([count, LAYERS, WINDOW]);
        let flat = Tensor::<BurnBackend, 1>::from_floats(windows, &self.device);
        let input = flat.reshape(shape);

        let forward = match opts.polymerase {
            PbChem::Two => self.two_zero.forward(input),
            PbChem::TwoPointTwo => self.two_two.forward(input),
            PbChem::ThreePointTwo => self.three_two.forward(input),
            PbChem::Revio => self.revio.forward(input),
        };
        let z: Vec<f32> = forward.to_data().value.chunks(2).map(|c| c[0]).collect();
        assert_eq!(z.len(), count);
        z
    }
}

impl Default for BurnModels {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_revio_onnx() {
        let device = BurnDevice::default();
        let model: revio::Model<BurnBackend> = revio::Model::default();
        let input = Tensor::<BurnBackend, 3>::zeros([2, LAYERS, WINDOW], &device);
        let output = model.forward(input);
        let z: Vec<f32> = output.to_data().value.chunks(2).map(|c| c[0]).collect();
        println!("{:?}", z);
    }
}
