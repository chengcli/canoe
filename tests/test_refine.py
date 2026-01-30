import torch
import torch.nn.functional as F

def refine_last_two(tensor, method):
    # Save original shape
    orig_shape = list(tensor.shape)
    n, m = orig_shape[-2], orig_shape[-1]

    # Collapse all leading dims into a "batch" dim for F.interpolate
    # Shape becomes (-1, 1, n, m)
    flattened = tensor.reshape(-1, 1, n, m)

    # Upsample
    refined = F.interpolate(flattened, scale_factor=2, mode=method)

    # Reshape back to (..., 2n, 2m)
    new_shape = orig_shape[:-2] + [2*n, 2*m]
    return refined.reshape(new_shape)

def coarsen_last_two(tensor):
    # Save original shape
    orig_shape = list(tensor.shape)
    n, m = orig_shape[-2], orig_shape[-1]

    # Collapse all leading dims into a "batch" dim for F.interpolate
    # Shape becomes (-1, 1, n, m)
    flattened = tensor.reshape(-1, 1, n, m)

    # Downsample
    coarsened = F.interpolate(flattened, scale_factor=0.5, mode='area')

    # Reshape back to (..., n/2, m/2)
    new_shape = orig_shape[:-2] + [n//2, m//2]
    return coarsened.reshape(new_shape)

if __name__ == "__main__":
    # Test the function
    x = torch.randn(3, 3)  # Example tensor with shape (3, 4, 5, 6)
    for i in range(3):
        for j in range(3):
            x[i, j] = i + j

    y1 = refine_last_two(x, "bilinear")
    x1 = coarsen_last_two(y1)
    dy = refine_last_two(x - x1, "area")
    y = y1 + dy
    z = coarsen_last_two(y)

    print("Original x :", x)
    print("Refined y :", y)
    print("Coarsened z :", z)
