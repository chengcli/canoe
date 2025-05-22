import os

def replace_in_file(file_path):
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()

    new_content = content.replace('configure.hpp', 'configure.h')

    if new_content != content:
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(new_content)
        print(f"Updated: {file_path}")

def replace_in_project(root_dir="."):
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.endswith(".cpp") or filename.endswith(".hpp"):
                full_path = os.path.join(dirpath, filename)
                replace_in_file(full_path)

if __name__ == "__main__":
    replace_in_project()
