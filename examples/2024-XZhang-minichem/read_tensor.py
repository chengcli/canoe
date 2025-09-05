import torch

# load the first file to infer shapes
module = torch.jit.load(pt_files[0][1])
data = {name: param for name, param in module.named_parameters()}

print(data.keys())
print(data["hydro_u/0"].shape)
print(data["scalar_s/1"].shape)
