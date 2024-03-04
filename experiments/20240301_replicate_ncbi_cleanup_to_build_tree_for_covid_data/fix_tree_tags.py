import argparse
import os

def fix_special_characters(input_file, output_file):
    with open(input_file, "r") as file:
        file_contents = file.read()

    #updated_contents = file_contents.translate(str.maketrans("ÖÄ", "{}"))
    updated_contents = file_contents.replace("-_-", "{").replace("._.", "}")

    with open(output_file, "w") as output_file:
        output_file.write(updated_contents)

def main():
    parser = argparse.ArgumentParser(description='Fix special characters in a file.')
    parser.add_argument('input_file', help='Path to the input file')
    parser.add_argument('output_file', help='Path to the output file')

    args = parser.parse_args()

    fix_special_characters(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
