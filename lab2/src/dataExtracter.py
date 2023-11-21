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
        if "MST" in line:
            MST = extract_numeric_data(line)
        elif "Best LS" in line:
            best_LS = extract_numeric_data(line)
        elif "Average LS" in line:
            avg_LS = extract_numeric_data(line)
        elif "inverts" in line:
            avg_inverts = extract_numeric_data(line)
    
    return (MST, best_LS, avg_LS, avg_inverts)


def main():
    script_directory = os.path.dirname(os.path.abspath(__file__))
    input_folder = os.path.join(script_directory, "..", "res\experiment2Data")
    output_folder = os.path.join(script_directory, "..", "res\experiment2Data\organizedData")


    print(input_folder)
    print(output_folder)

    # if output folder does not exist, create it
    os.makedirs(output_folder, exist_ok=True)

    input_files = [file for file in os.listdir(input_folder) if os.path.isfile(os.path.join(input_folder, file))]


    MST_data_list = []
    best_LS_list = []
    avg_LS_list = []
    avg_inverts_list = []

    for input_file_name in input_files:
        input_file_path = os.path.join(input_folder, input_file_name)
        print("Reading from: " + input_file_path)

        MST, best_LS, avg_LS, avg_inverts = process_file(input_file_path)

        MST_data_list.append(MST)
        best_LS_list.append(best_LS)
        avg_LS_list.append(avg_LS)
        avg_inverts_list.append(avg_inverts)

    MST_data_list.sort()
    best_LS_list.sort()
    avg_LS_list.sort()
    avg_inverts_list.sort()


    # output_file = "MST_sorted.txt"
    # output_file_path = os.path.join(output_folder, output_file)
    # safe_data(output_file_path, MST_data_list)

    output_file = "Best_LS_sorted.txt"
    output_file_path = os.path.join(output_folder, output_file)
    safe_data(output_file_path, best_LS_list)

    output_file = "Avg_LS_sorted.txt"
    output_file_path = os.path.join(output_folder, output_file)
    safe_data(output_file_path, avg_LS_list)

    output_file = "Avg_inverts_sorted.txt"
    output_file_path = os.path.join(output_folder, output_file)
    safe_data(output_file_path, avg_inverts_list)

if __name__ == "__main__":
    main()
