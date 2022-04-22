import pandas as pd
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Divide Test Set')
    parser.add_argument('--data', required=True)
    parser.add_argument('--data_sep', default='\t')
    parser.add_argument('--header', default=None)
    parser.add_argument('--save_dir', default='./')
    args = parser.parse_args()

data = pd.read_csv(args.data, sep=args.data_sep, header=args.header)
file_name = args.data.split('/')[-1].split('.')[0]

data[data.columns[0]].to_csv(f'{args.save_dir}{file_name}_X.txt', sep='\t', index=False, header=False)
data[data.columns[1]].to_csv(f"{args.save_dir}{file_name}_y.txt", sep='\t', index=False, header=False)

