import os
import numpy as np
import matplotlib.pyplot as plt
import imageio




with imageio.get_writer('mygif.gif', mode='I') as writer:
    n = 20
    Names = ['2.png'] * 2*n + ['3.png'] * n + ['4.png'] * n + ['5.png'] * n + \
            ['6.png'] * n + ['7.png'] * n + ['8.png'] * 3*n
    for filename in Names:
        image = imageio.imread(filename)
        writer.append_data(image)
