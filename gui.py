# Standard packages
import Tkinter, tkFileDialog

# Third party packages
import matplotlib.pyplot as plt

def pickdir():
    """GUI directory picker.
    
    Uses the Tk module to provide the user with a simple way of selecting
    a directory for processing. Quits if no directory is chosen. Takes no
    inputs and returns a unicode string with the selected directory.
    """
    root = Tkinter.Tk()
    root.withdraw()
    top_dir = tkFileDialog.askdirectory(parent=root, initialdir="..", 
                                        title='Pick a directory')
    if top_dir is '':
        print 'No directory provided, quitting...'
        sys.exit(0)
    return top_dir
    
def contour(data):
    plt.contourf(data)
    show()
    
def slice(**sets):
    plt.hold(True)
    for k in sets.keys():
        plt.plot(sets[k])
    plt.legend(sets.keys())
    plt.show()
    plt.hold(False)