file = 'C:/Users/PARK/Downloads/HPK1_pose.sd'

with open(file, 'r') as f:
    lines = f.readlines()
    name_count = {}
    new_lines = ''
    name_idx = 0

    for idx, line in enumerate(lines):

        if idx == name_idx:
            name = line.rstrip()

            if name in name_count.keys():
                name_count[name] += 1
            else:
                name_count[name] = 1

            name = f'{name}_{name_count[name]}'
            new_lines += line
            continue

        if line.startswith('> <Name>'):
            print(lines[idx + 1].rstrip())
            lines[idx + 1] = f'{name}\n'

        new_lines += line

        if line == '$$$$\n':
            name_idx = idx+1
            with open(f'C:/Users/PARK/Downloads/HPK1_ligands/{name}.sdf', 'w') as f:
                f.write(new_lines)
            new_lines = ''

