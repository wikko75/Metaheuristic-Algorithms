import os
import re

def extract_numeric_data(content):
    numeric_data = re.findall(r'\b\d+\b', content)
    return int(numeric_data[0])

def safe_data(file_path, data):
    with open(file_path, 'a') as file:
        for number in data:
            file.write(str(number) + '\n')


def process_file(input_file_path):

    with open(input_file_path, 'r') as input_file:
        lines = input_file.readlines()

    for line in lines:
        if "Best GenTSP" in line:
            GENOPC = extract_numeric_data(line)
    
    return (GENOPC)


def main():
    script_directory = os.path.dirname(os.path.abspath(__file__))
    input_folder = os.path.join(script_directory, "..", "res/resultsTPC")
    output_folder = os.path.join(script_directory, "..", "res/organizedTPC")


    print(input_folder)
    print(output_folder)

    # if output folder does not exist, create it
    os.makedirs(output_folder, exist_ok=True)

    input_files = [file for file in os.listdir(input_folder) if os.path.isfile(os.path.join(input_folder, file))]


    GenOPC_DATA = []


    for input_file_name in input_files:
        input_file_path = os.path.join(input_folder, input_file_name)
        print("Reading from: " + input_file_path)

        GENOPC = process_file(input_file_path)

        GenOPC_DATA.append(GENOPC)

    GenOPC_DATA.sort()

    output_file = "GEN_TPC_sorted.txt"
    output_file_path = os.path.join(output_folder, output_file)
    safe_data(output_file_path, GenOPC_DATA)

if __name__ == "__main__":
    main()