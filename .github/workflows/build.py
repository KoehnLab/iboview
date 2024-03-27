#!/usr/bin/env python3

import os
import argparse
import re
import subprocess

def execute_bash(script: str):
    with open("script.sh", "w") as file:
        file.write("set -e\nset -x\n")
        file.write(script)

    subprocess.check_call(["bash", "./script.sh"])

    os.remove("script.sh")

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--os")
    parser.add_argument("--only-install-dependencies", action="store_true", default=False)
    parser.add_argument("--only-compile", action="store_true", default=False)

    args = parser.parse_args()

    # We use docker image names as OS spec -> strip it down to the pure OS name
    os_name = args.os
    if "/" in os_name:
        os_name = os_name[ : os_name.index("/") ]
    if ":" in os_name:
        os_name = os_name[ : os_name.index(":")]
    os_name = os_name.strip()
    assert len(os_name) > 0

    script_dir = os.path.dirname(os.path.realpath(__file__))
    repo_root = os.path.join(script_dir, "..", "..")
    readme_path = os.path.join(repo_root, "README.md")

    assert os.path.isfile(readme_path)

    contents = open(readme_path, "r").read()

    if not args.only_compile:
        match = re.search(r"#\s*install(ing)?\b\s*dependenc(y|ies)", contents, re.IGNORECASE)
        assert match is not None

        dependency_instruction_start = match.start()
        os_idx = contents.lower().index(os_name, dependency_instruction_start)

        install_command_start = contents.find("```", os_idx)
        install_command_start = contents.find("\n", install_command_start)
        install_command_end = contents.find("```", install_command_start)

        install_command = contents[install_command_start : install_command_end].strip()

        # Prevent any interactive yes/no prompts
        install_command = re.sub(r"\binstall\b", "install -y", install_command)
        install_command = re.sub(r"\bremove\b", "remove -y", install_command)

        execute_bash(install_command)

    if not args.only_install_dependencies:
        match = re.search(r"#\s*compil(e|ing)", contents, re.IGNORECASE)
        assert match is not None

        compile_section_start = match.start()

        compile_command_start = contents.find("```", compile_section_start)
        compile_command_start = contents.find("\n", compile_command_start)
        compile_command_end = contents.find("```", compile_command_start)

        compile_commands = contents[compile_command_start : compile_command_end].strip()

        execute_bash(compile_commands)




if __name__ == "__main__":
    main()
