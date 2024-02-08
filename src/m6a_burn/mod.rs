use super::bio_io::PbChem;
use super::predict_m6a::PredictOptions;
use burn::tensor::backend::Backend;
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

use burn::backend::wgpu::WgpuDevice;
use burn::backend::Wgpu;
type BurnBackend = Wgpu;
type BurnDevice = WgpuDevice;
*/

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

/// B is for the burn backend and D is for the device
#[derive(Debug)]
pub struct BurnModels<B>
where
    B: Backend,
{
    pub two_zero: two_zero::Model<B>,
    pub two_two: two_two::Model<B>,
    pub three_two: three_two::Model<B>,
    pub revio: revio::Model<B>,
    pub device: B::Device,
}

impl<B> BurnModels<B>
where
    B: Backend,
{
    pub fn new() -> Self {
        let two_zero = two_zero::Model::default();
        let two_two = two_two::Model::default();
        let three_two = three_two::Model::default();
        let revio = revio::Model::default();
        let device = B::Device::default();

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

    pub fn forward(&self, opts: &PredictOptions<B>, windows: &[f32], count: usize) -> Vec<f32> {
        // predict
        if count == 0 {
            return vec![];
        }
        let shape = Shape::new([count, LAYERS, WINDOW]);
        let flat = Tensor::<B, 1, burn::tensor::Float>::from_floats(windows, &self.device);
        let input = flat.reshape(shape);

        let forward: Tensor<B, 2, burn::tensor::Float> = match opts.polymerase {
            PbChem::Two => self.two_zero.forward(input),
            PbChem::TwoPointTwo => self.two_two.forward(input),
            PbChem::ThreePointTwo => self.three_two.forward(input),
            PbChem::Revio => self.revio.forward(input),
        };

        let z: Vec<f32> = forward
            .into_data()
            .convert()
            .value
            .chunks(2)
            .map(|c| c[0])
            .collect();
        z
    }
}

impl<B> Default for BurnModels<B>
where
    B: Backend,
{
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use burn::backend::ndarray::NdArrayDevice;
    use burn::backend::NdArray;
    pub type BurnBackend = NdArray;
    pub type BurnDevice = NdArrayDevice;

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
