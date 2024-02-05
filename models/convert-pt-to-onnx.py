#!/usr/bin/env python3
import torch
import onnx
import glob

input_example = (1, 6, 15)
device = "cpu"

model_files = glob.glob("models/*semi*.pt")
for model_file in model_files:
    if "revio" in model_file.lower():
        continue
    print(f"Converting {model_file} to ONNX")
    model = torch.jit.load(model_file)
    print(model)

    example = torch.rand(input_example).float().to(device)
    input_names = ["actual_input_1"]
    output_names = ["output1"]

    torch.onnx.export(
        model,
        example,
        f"{model_file}.onnx",
        verbose=True,
        input_names=input_names,
        output_names=output_names,
    )
