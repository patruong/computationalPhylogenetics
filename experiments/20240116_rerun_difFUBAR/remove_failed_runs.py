import pandas as pd
import os 
import shutil
os.chdir("/home/patrick/git/computationalPhylogenetics/experiments/20240116_rerun_difFUBAR/output")

failed_folders = []
success_folders = []
for folder in os.listdir():
    res_file = os.path.join(folder, "patrick_max_child", "_posteriors.csv")
    if os.path.exists(res_file):
        print(res_file)
        success_folders.append(folder)
    else:
        failed_folders.append(folder)
    
def delete_folder(folder_to_delete):
    try:
        shutil.rmtree(folder_to_delete)
        print(f'The folder {folder_to_delete} has been deleted successfully.')
    except FileNotFoundError:
        print(f'The folder {folder_to_delete} does not exist.')
    except Exception as e:
        print(f'An error occurred: {e}')

for folder in failed_folders:
    delete_folder(folder)

print(len(failed_folders)) 
print(len(success_folders)) #2415
