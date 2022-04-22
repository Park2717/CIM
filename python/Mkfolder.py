import os
import glob
import shutil

'''
Files deleted with 'rmtree' in shutil module cannot be recovered.
So.. Please be careful when using this code!!
Check & Recheck your new folder path.
'''

# search directory
path = glob.glob('./*')
print(path)

# mk new folder
def newfolder(directory):
    try:
        if not os.path.exists(directory):
            print(f'Does not exist output folder in {directory}')
            check = input('Do you want create new folder? (y / n) :')
            if check == 'y' or check == 'Y':
                os.makedirs(directory)
                print(f'Create New Folder {directory}.')
            else:
                print('Code Stop.')
        else:
            print(f'Already exist output folder in {directory}')
            rm_check = input('Do you want recreate new folder? < delete file > (y / n) :')
            if rm_check == 'y' or rm_check == 'Y':
                try:
                    shutil.rmtree(directory)
                    os.mkdir(directory)
                    print(f'Recreate New Folder {directory}.')
                except Exception as e:
                    print(e)
            else:
                print('Code Stop.')
    except OSError:
        print('Error')

# Set Directory & Run
new_folder = 'C:/Users/PARK/Desktop/test2/new'
newfolder(new_folder)


