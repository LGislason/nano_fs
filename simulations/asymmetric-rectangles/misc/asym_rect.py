import meep as mp
import matplotlib.pyplot as plt

resolution = 40
cell = mp.Vector3(1, 1, 0)

geometry = [
    mp.Block(
        size=mp.Vector3(0.3, 0.15, mp.inf),
        center=mp.Vector3(0.2, -0.1, 0),
        material=mp.Medium(epsilon=4)
    )
]

sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    resolution=resolution
)

sim.init_sim()

eps = sim.get_array(
    center=mp.Vector3(),
    size=cell,
    component=mp.Dielectric
)

plt.imshow(eps.T, origin="lower")
plt.colorbar()
plt.title("Asymmetric Rectangle ε")
plt.show()

