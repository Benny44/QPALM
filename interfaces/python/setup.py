import os
import platform

if __name__ == '__main__':
    print("Running setup")
    python_interface_dir = os.getcwd()
    #os.chdir("../")
    #root_dir = os.getcwd()
    #print(root_dir)
    #qpalm_static_lib_dir = root_dir + "/build/lib/static/"
    if not os.path.isdir('build'):
        os.mkdir('build')
        os.mkdir('build/lib')
        os.mkdir('build/include')

    cwd = os.getcwd()
    
    os.chdir("../")
    if not os.path.isdir('build'):
        os.mkdir('build')
        os.mkdir('build/lib')
        os.mkdir('build/include')

    os.chdir(cwd + '/build')

    os.system('cmake .. -DCMAKE_BUILD_TYPE=release -DCOVERAGE=OFF')
    os.system('make')

