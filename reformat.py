import Tkinter, tkFileDialog
from os import listdir, rename
from os.path import join

root = Tkinter.Tk()
root.withdraw()
top_dir = tkFileDialog.askdirectory(parent=root, initialdir=dir, 
                                    title='Pick a directory...')
if top_dir is '':
    print 'No directory provided, quitting...'
    sys.exit(0)

files = listdir(top_dir)

for f in files:
    if 'Background' in f:
        fn = f.replace('Background', 'Reference')
        rename(join(top_dir, f), join(top_dir, fn))
    else:
        pass

print 'Done!'