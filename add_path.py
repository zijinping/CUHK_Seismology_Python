import os
from distutils.sysconfig import get_python_lib

#Set path
print(">>> Add path to python library ...")

pwd = os.getcwd()
lib_path = get_python_lib()
path_file = os.path.join(lib_path,'relocation.pth')
with open(path_file,'w') as f:
    f.write(pwd)
f.close()

print("Done!")
