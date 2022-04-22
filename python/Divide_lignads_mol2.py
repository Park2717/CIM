file = 'C:/Users/PARK/Downloads/HPK1_pose_mol2.mol2'

with open(file, 'r') as f:
    lines = f.readlines()
    name_count = {}
    new_lines = ''
    block = []

    for idx, line in enumerate(lines):
        if line.startswith('@<TRIPOS>MOLECULE'):
            name = lines[idx + 1].rstrip()

            if name == '6NFY_A':
                continue

            block += [idx]

            if name in name_count.keys():
                name_count[name] += 1
            else:
                name_count[name] = 1

            name = f'{name}_{name_count[name]}'
            lines[idx + 1] = f'{name}\n'

    for idx, mol in enumerate(block):
        if idx == len(block) - 1:
            new_lines = lines[mol:]
        else:
            new_lines = lines[mol : block[idx + 1]]
        with open(f'C:/Users/PARK/Downloads/HPK1_ligands_mol2/{new_lines[1].rstrip()}.mol2', 'w') as f:
            for ln in new_lines:
                f.write(ln)
