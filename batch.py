import os

model_dirs = [d for d in os.scandir() if os.path.isdir(d.name)]
input_files = [f.path for d in model_dirs for f in os.scandir(d)]

