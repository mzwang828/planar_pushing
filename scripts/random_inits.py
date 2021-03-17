import numpy as np
import yaml
from pathlib import Path

n = 10

# random perturbed initials for [x,y,theta,phi] in range +_[0.03, 0.03, 0.4, 0.5]
random_inits = np.random.rand(n,4)
random_inits[:, 0] -= 0.50
random_inits[:, 0] *= 0.06
random_inits[:, 1] -= 0.50
random_inits[:, 1] *= 0.06
random_inits[:, 2] -= 0.50
random_inits[:, 2] *= 0.8
random_inits[:, 3] -= 0.50

random_inits = np.round(random_inits, 6)

data = {"init": random_inits.tolist()}
fname = Path(__file__).parent.parent.absolute() / 'Config' / 'initials.yaml'
with open(fname, "w") as f:
    yaml.dump(data, f, default_flow_style = None)