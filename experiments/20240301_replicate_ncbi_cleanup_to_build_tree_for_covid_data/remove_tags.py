import re
import argparse
import os

def remove_tags(input_file, output_file):
    with open(input_file, "r") as file:
        file_contents = file.read()

    #updated_contents = file_contents.translate(str.maketrans("ÖÄ", "{}"))
    pattern = r'\{.*?\}'

    cleaned_text = re.sub(pattern, '', file_contents)

    with open(output_file, "w") as output_file:
        output_file.write(cleaned_text)

def main():
    parser = argparse.ArgumentParser(description='Remove tags from file.')
    parser.add_argument('input_file', help='Path to the input file')
    parser.add_argument('output_file', help='Path to the output file')

    args = parser.parse_args()

    remove_tags(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
