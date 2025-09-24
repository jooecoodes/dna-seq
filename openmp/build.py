import os
import subprocess

# List of C++ source files to compile
cpp_files = [
    "BM.cpp",
    "BP.cpp",
    "KMP.cpp",
]

# Compiler and flags
compiler = "g++"
flags = ["-fopenmp", "-O3"]

# Compile each file
for cpp in cpp_files:
    if not os.path.exists(cpp):
        print(f"⚠️ Skipping {cpp} (file not found)")
        continue

    exe = os.path.splitext(cpp)[0]  # "BM.cpp" -> "BM"
    cmd = [compiler] + flags + [cpp, "-o", exe]

    print(f"🚀 Compiling {cpp} -> {exe} ...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode == 0:
        print(f"✅ {exe} built successfully!\n")
    else:
        print(f"❌ Failed to compile {cpp}:\n{result.stderr}\n")
