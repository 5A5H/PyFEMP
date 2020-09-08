# Setup script for PyFEMP
import os, sys, shutil

if __name__ == "__main__":
    print('Execute setup.py')
    
    major, minor, micro, releaselevel, serial = sys.version_info
    
    # find installation directory
    installdir = os.path.join(sys.prefix, "lib", "python"+str(major)+"."+str(minor), "site-packages")
    if os.path.exists(installdir): 
        installdir = os.path.join(installdir, "PyFEMP")
    else:
        installdir = os.path.join(sys.path[2], "PyFEMP")

    # check if module is here
    module = os.path.join(os.path.dirname(__file__), "PyFEMP_src")
    if not os.path.exists(module): 
        raise NameError('\033[91m'+"Error! Redownload and try again!")


    # check the environment
    existingPyFEMP = ""
    try:
        import numpy
    except ImportError:
        raise NameError('\033[91m'+"Error! Install the python package numpy first.")
    try:
        import matplotlib
    except ImportError:
        raise NameError('\033[91m'+"Error! Install the python package matplotlib first.")
    try:
        import PyFEMP
        existingPyFEMP = os.path.dirname(PyFEMP.__file__)
    except:
        pass

    print('Current path     : ',os.getcwd())
    print('Python executable: ',sys.executable)
    print('Python path      : ',sys.prefix)
    print('Python version   : ',str(major)+"."+str(minor)+"."+str(micro))
    print('Install dir      : ','\033[92m'+installdir+'\033[0m')

    # check what to do
    if (len(sys.argv) <2):
        
        # if present delete PyFEMP
        if (existingPyFEMP!=""):
            print('... uninstalling PyFEMP')
            shutil.rmtree(existingPyFEMP)
        
        # install = copy the module
        print('... installing PyFEMP')
        shutil.copytree(module, installdir)
        try:
            import PyFEMP
        except ImportError:
            raise NameError('\033[91m'+"Error! Installation fialed! Redownload and try again.")
        print('\033[92m'+"Installation successful!!"+'\033[0m')


    elif (len(sys.argv) ==2):
        # check the argument
        if (sys.argv[1]=='install'):
        
            # if present delete PyFEMP
            if (existingPyFEMP!=""):
                print('... uninstalling PyFEMP')
                shutil.rmtree(existingPyFEMP)

            print('... installing PyFEMP')
            # install = copy the module
            shutil.copytree(module, installdir)
            try:
                import PyFEMP
            except ImportError:
                raise NameError('\033[91m'+"Error! Installation fialed! Redownload and try again.")
            print('\033[92m'+"Installation successful!!"+'\033[0m')
        
        elif (sys.argv[1]=='uninstall'):

            # if present delete PyFEMP
            if (existingPyFEMP!=""):
                print('... uninstalling PyFEMP')
                shutil.rmtree(existingPyFEMP)
            else:
                raise NameError('\033[91m'+"Error! PyFEMP could not be found and thus not uninstalled!")

        else:
            raise NameError('\033[91m'+"Error! This script supports arguments [install, uninstall]!")
    else:
        raise NameError('\033[91m'+"Error! This script needs one argument only [install, uninstall]!")
    