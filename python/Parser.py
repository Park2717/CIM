import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parser practice')
    parser.add_argument('--train', required=True)
    parser.add_argument('--validation', required=True)
    parser.add_argument('--test', default='test.txt')
    parser.add_argument('--save_dir', required=True)
    parser.add_argument('--load_model', default=None)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

