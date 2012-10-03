# Standard packages
import sys
import Tkinter, tkFileDialog

# Third party packages
import matplotlib.pyplot as plt
import numpy as np

def pickdir(text, dir='..'):
    """GUI directory picker.
    
    Uses the Tk module to provide the user with a simple way of selecting
    a directory for processing. Quits if no directory is chosen. Takes no
    inputs and returns a unicode string with the selected directory.
    """
    root = Tkinter.Tk()
    root.withdraw()
    top_dir = tkFileDialog.askdirectory(parent=root, initialdir=dir, 
                                        title=text)
    if top_dir is '':
        print 'No directory provided, quitting...'
        sys.exit(0)
    return top_dir
    
def contour(data):
    levels = np.linspace(np.min(data), np.max(data), 100)
    plt.contourf(data, levels=levels)
    plt.colorbar()
    plt.show()
    
def slice(**sets):
    """ Plots a variable number of dataset slices.
    
    Takes an unrolled dictionary of datasets (keywords are the legend
    titles) and plots them all on a single axis. Useful for comparing
    different slices or processing methods.
    """
    plt.hold(True)
    for k in sets.keys():
        plt.plot(sets[k])
    plt.legend(sets.keys(), loc=4)
    plt.show()
    plt.hold(False)