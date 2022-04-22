import os
import pandas as pd
import numpy as np
import shutil
import argparse
from sklearn.model_selection import train_test_split

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Split train_test_valid')
    parser.add_argument('--data', required=True)
    parser.add_argument('--data_sep', default='\t')
    parser.add_argument('--valid_size', type=int, default=0.2)
    parser.add_argument('--test_size', type=int, default=0.2)
    parser.add_argument('--save_dir', default='./')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

data = pd.read_csv(args.data, sep=args.data_sep)
file_name = args.data.split('/')[-1].split('.')[0]

valid_cal_size = args.valid_size / (1 - args.test_size)

protrain, test = train_test_split(data, test_size=args.test_size, random_state=args.seed, shuffle=True)
train, val = train_test_split(protrain, test_size=valid_cal_size, random_state=args.seed, shuffle=True)

train.to_csv(f"{args.save_dir}{file_name}_train.txt", sep='\t', index=False, header=False)
val.to_csv(f"{args.save_dir}{file_name}_val.txt", sep='\t', index=False, header=False)
test.to_csv(f"{args.save_dir}{file_name}_test.txt", sep='\t', index=False, header=False)

