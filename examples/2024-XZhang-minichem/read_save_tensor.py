import torch


def save_tensors(tensor_map: dict[str, torch.Tensor], filename: str):
    class TensorModule(torch.nn.Module):
        def __init__(self, tensors):
            super().__init__()
            for name, tensor in tensors.items():
                self.register_buffer(name, tensor)

    module = TensorModule(tensor_map)
    scripted = torch.jit.script(module)  # Needed for LibTorch compatibility
    scripted.save(filename)


if __name__ == "__main__":
    tensors = {
        "foo": torch.randn(3, 4),
        "bar": torch.randn(5, 6),
    }
    save_tensors(tensors, "foo_bar.pt")

    # load the first file to infer shapes
    module = torch.jit.load(pt_files[0][1])
    data = {name: param for name, param in module.named_parameters()}

    print(data.keys())
    print(data["hydro_u/0"].shape)
    print(data["scalar_s/1"].shape)
