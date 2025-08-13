use super::subcommands::predict_m6a::PredictOptions;
use super::subcommands::predict_m6a::{LAYERS, WINDOW};
use crate::utils::bio_io::PbChem;
use burn::module::Module;
use burn::tensor::backend::Backend;
use burn::tensor::{Shape, Tensor};
use std::sync::Once;

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

static LOG_DEVICE_ONCE: Once = Once::new();

/// B is for the burn backend and D is for the device
#[derive(Debug, Clone)]
pub struct BurnModels<B>
where
    B: Backend<Device = BurnDevice>,
{
    pub two_zero: Option<two_zero::Model<B>>,
    pub two_two: Option<two_two::Model<B>>,
    pub three_two: Option<three_two::Model<B>>,
    pub revio: Option<revio::Model<B>>,
    pub device: BurnDevice,
}

impl<B> BurnModels<B>
where
    B: Backend<Device = BurnDevice>,
{
    pub fn new(polymerase: &PbChem) -> Self {
        #[cfg(not(feature = "tch"))]
        let device = B::Device::default();
        #[cfg(feature = "tch")]
        let device = Self::get_libtorch_device();

        // log info about the device used
        LOG_DEVICE_ONCE.call_once(|| {
            log::info!("Using {device:?} for Burn device.");
        });

        match polymerase {
            PbChem::Two => {
                let two_zero = two_zero::Model::default().to_device(&device);
                Self {
                    two_zero: Some(two_zero),
                    two_two: None,
                    three_two: None,
                    revio: None,
                    device,
                }
            }
            PbChem::TwoPointTwo => {
                let two_two = two_two::Model::default().to_device(&device);
                Self {
                    two_zero: None,
                    two_two: Some(two_two),
                    three_two: None,
                    revio: None,
                    device,
                }
            }
            PbChem::ThreePointTwo => {
                let three_two = three_two::Model::default().to_device(&device);
                Self {
                    two_zero: None,
                    two_two: None,
                    three_two: Some(three_two),
                    revio: None,
                    device,
                }
            }
            PbChem::Revio => {
                let revio = revio::Model::default().to_device(&device);
                Self {
                    two_zero: None,
                    two_two: None,
                    three_two: None,
                    revio: Some(revio),
                    device,
                }
            }
        }
    }

    #[cfg(feature = "tch")]
    fn get_libtorch_device() -> BurnDevice {
        use burn::backend::libtorch::LibTorchDevice;
        let device = if tch::utils::has_cuda() {
            LibTorchDevice::Cuda(0)
        }
        //else if tch::utils::has_mps() {
        //LibTorchDevice::Mps
        //}
        else {
            LibTorchDevice::Cpu
        };
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

        // allow fake predictions for testing speed of other parts of the code
        // I have moved this chunk around and it is indeed just the model.forward
        // step that is slow
        if opts.fake {
            return vec![0.0; count];
        }
        let forward: Tensor<B, 2, burn::tensor::Float> = match opts.polymerase {
            PbChem::Two => self.two_zero.as_ref().unwrap().forward(input),
            PbChem::TwoPointTwo => self.two_two.as_ref().unwrap().forward(input),
            PbChem::ThreePointTwo => self.three_two.as_ref().unwrap().forward(input),
            PbChem::Revio => self.revio.as_ref().unwrap().forward(input),
        };
        forward
            .into_data()
            .to_vec::<f32>()
            .unwrap()
            .chunks(2)
            .map(|c| c[0])
            .collect()
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
        let z: Vec<f32> = output
            .to_data()
            .to_vec::<f32>()
            .unwrap()
            .chunks(2)
            .map(|c| c[0])
            .collect();
        println!("{z:?}");
    }
}
