import numpy as np


for name, value in {
    "object": object,
    "bool": bool,
    "int": int,
    "float": float,
}.items():
    if name not in np.__dict__:
        setattr(np, name, value)
