import os
import platform
import subprocess

def build_library():
    cpp_file = "dna_matchers.cpp"
    
    if platform.system() == "Windows":
        lib_name = "dna_matchers.dll"
        cmd = ["g++", "-shared", "-o", lib_name, cpp_file]
    else:
        lib_name = "dna_matchers.so"
        cmd = ["g++", "-shared", "-fPIC", "-o", lib_name, cpp_file]
    
    try:
        subprocess.check_call(cmd)
        print(f"Successfully compiled {lib_name}")
    except subprocess.CalledProcessError as e:
        print(f"Compilation failed: {e}")
    except FileNotFoundError:
        print("g++ compiler not found. Please install MinGW on Windows or build-essential on Linux.")

if __name__ == "__main__":
    build_library()