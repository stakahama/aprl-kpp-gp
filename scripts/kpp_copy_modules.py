
import os
import shutil

modulepath = os.path.join(os.path.dirname(__file__),'modules')

for mod in os.listdir(modulepath):
    shutil.copy2(os.path.join(modulepath,mod),'.')
    print mod+' copied'
