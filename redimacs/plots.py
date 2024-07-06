import matplotlib.pyplot as plt
from astropy.visualization import simple_norm


def show_image(array, output_path):
    plt.figure()
    norm = simple_norm(array, stretch='asinh', percent=99.9, asinh_a=0.05)
    plt.imshow(array, origin='lower', norm=norm)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, format='jpg', dpi=400)
    plt.close()
