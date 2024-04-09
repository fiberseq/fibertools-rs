use super::bio_io::PbChem;
use super::predict_m6a::PredictOptions;
use super::predict_m6a::{LAYERS, WINDOW};
use burn::module::Module;
use burn::tensor::backend::Backend;
use burn::tensor::{Shape, Tensor};

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

#[cfg(feature = "tch")]
pub type BurnDevice = burn::backend::libtorch::LibTorchDevice;
#[cfg(not(feature = "tch"))]
pub type BurnDevice = burn::backend::candle::CandleDevice;

/// B is for the burn backend and D is for the device
#[derive(Debug)]
pub struct BurnModels<B>
where
    B: Backend<Device = BurnDevice>,
{
    pub two_zero: two_zero::Model<B>,
    pub two_two: two_two::Model<B>,
    pub three_two: three_two::Model<B>,
    pub revio: revio::Model<B>,
    pub device: BurnDevice,
}

impl<B> BurnModels<B>
where
    B: Backend<Device = BurnDevice>,
{
    pub fn new() -> Self {
        #[cfg(not(feature = "tch"))]
        let device = B::Device::default();
        #[cfg(feature = "tch")]
        let device = Self::get_libtorch_device();

        let two_zero = two_zero::Model::default().to_device(&device);
        let two_two = two_two::Model::default().to_device(&device);
        let three_two = three_two::Model::default().to_device(&device);
        let revio = revio::Model::default().to_device(&device);

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

    #[cfg(feature = "tch")]
    fn get_libtorch_device() -> BurnDevice {
        use burn::backend::libtorch::LibTorchDevice;
        let device = if tch::utils::has_cuda() {
            LibTorchDevice::Cuda(0)
        } else if tch::utils::has_mps() {
            LibTorchDevice::Mps
        } else {
            LibTorchDevice::Cpu
        };
        log::info!("Number of threads for Torch: {}", tch::get_num_threads());
        device
    }

    pub fn forward(&self, opts: &PredictOptions<B>, windows: &[f32], count: usize) -> Vec<f32> {
        // predict
        if count == 0 {
            return vec![];
        }

        let shape = Shape::new([count, LAYERS, WINDOW]);
        let input =
            Tensor::<B, 1, burn::tensor::Float>::from_floats(windows, &self.device).reshape(shape);

        let forward: Tensor<B, 2, burn::tensor::Float> = match opts.polymerase {
            PbChem::Two => self.two_zero.forward(input),
            PbChem::TwoPointTwo => self.two_two.forward(input),
            PbChem::ThreePointTwo => self.three_two.forward(input),
            PbChem::Revio => self.revio.forward(input),
        };

        forward
            .into_data()
            .convert()
            .value
            .chunks(2)
            .map(|c| c[0])
            .collect()
    }
}

impl<B> Default for BurnModels<B>
where
    B: Backend<Device = BurnDevice>,
{
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use burn::backend::candle::CandleDevice;
    use burn::backend::Candle;
    pub type BurnBackend = Candle;
    pub type BurnDevice = CandleDevice;

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
