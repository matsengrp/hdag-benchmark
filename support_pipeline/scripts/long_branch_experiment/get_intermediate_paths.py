import os
import glob

clade = "AY.34.2"
base_dir = "/n/fs/ragr-research/users/wh8114/prev/hdag-benchmark/data"
out_path = f"{base_dir}/sim_models/AY.34.2/intermediate_file_paths.txt"


def main():
    latest_intermediate_paths = []
    for trial in [1, 2, 3, 4, 5]:
        for bm in [1, 16, 256]:
            intermediate_dir = f"{base_dir}/sim_models/{clade}/{trial}/branch_multiplier_{bm}/results/historydag"

            # All files and directories ending with .txt and that don't begin with a dot:
            file_paths = glob.glob(f"{intermediate_dir}/final_opt_dag_trimmed.pb.intermediate.*")

            if len(file_paths) == 0:
                print(f"Skipping\ttrial: {trial}\tbm: {bm}...")
                continue

            # print(file_paths)
            latest_intermediate = ""
            latest_iter = 0
            for path in file_paths:
                last_idx = path.rfind('.')  # find last instance of `.`
                if last_idx == -1 or path.endswith("_dir"):
                    continue
                else:
                    iteration = int(path[last_idx+1:])
                    latest_iter = max(iteration, latest_iter)
                    if latest_iter == iteration:
                        latest_intermediate = path

            print(f"\ttrial: {trial}\tbm: {bm}\tlatest_iter: {latest_iter}")
            latest_intermediate_paths.append(f"{latest_intermediate}\n")
        
    with open(out_path, "w") as fp:
        fp.writelines(latest_intermediate_paths)



if __name__ == '__main__':
    main()