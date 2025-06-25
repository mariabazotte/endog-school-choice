import os

general_folder = "../results/"
folder_type   = ["deterministic/","stochastic/"]
folder_transf = ["nonstrategic/","compact/","enumeration/","simple/"]
folder_solver = ["modlconstraints/","cutoffscore/","asgheuristic/","lagrangian/","localsearch","integerlshaped","simulannealing"]

# List of folder names to create
folders = [general_folder + ty + tr + sol for ty in folder_type for tr in folder_transf for sol in folder_solver]

# Loop through the folder names and create each one
for folder in folders:
    try:
        os.makedirs(folder, exist_ok=True)
        print(f"Folder '{folder}' created successfully.")
    except FileExistsError:
        print(f"Folder '{folder}' already exists.")
