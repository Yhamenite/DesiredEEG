import os
import sys
import csv
import pyfiglet
import pandas as pd
from time import strftime, localtime


## Global variables and dictionaries
NUMBER_OF_SESSIONS = 16
NUMBER_OF_TRAILS_PER_SESSION = 5
NUMBER_OF_RECORDS_PER_TRIAL = 2
NUMBER_OF_MARKERS_PER_TRIAL = 4
NUMBER_OF_MARKERS_PER_SESSION = NUMBER_OF_MARKERS_PER_TRIAL * NUMBER_OF_TRAILS_PER_SESSION
OUTPUT_DIRECTORY_NAME = "\\Output"

subject = {
    "mj": 0,
    "yamen": 1,
    "nuw":2
}

age = {
    "mj":24,
    "yamen":24,
    "nuw":23
}

gender = {
    "mj":0,
    "yamen":0,
    "nuw":0
}

word = {
    "yes":0,
    "no":1,
    "hellfire":2,
    "helfire":2,
    "paradise":3
}

data_class = {
    "spoken":0,
    "imagined":1
}


def main():
   
   # The program starts executing from here
   ascii_art()
   path = get_path()
   print(to_color("[*]","blue"), "Console Log:")
   print(to_color("[+]","green"), "Path is " + path)
   files = list_files(path)
   marker_files, data_files = seperate_marker_files(files)
   processing(marker_files, data_files,path)
   

def to_color(string, color):
    color_code = {  'blue': '\033[34m',
                    'yellow': '\033[33m',
                    'green': '\033[32m',
                    'red': '\033[31m'
                    }
    return color_code[color] + str(string) + '\033[0m'


def ascii_art():
    art = pyfiglet.figlet_format("DesiredEEG")
    print(art)
    print("EEG data preprocessing tool made by:", to_color("GitHub.con/yhamenite",'blue'))
    print()
    

def get_path():
    # Get the path of where the files are located, even though it's best to put the script within the files that will be processed
    path = input(f"Please provide the path of the folder where 'intervalMarker' and main data files reside, {to_color('and please note','red')} that the output folder will be made in this provided path: ")
    print()
    return path


def list_files(path):
    # This function check to see if there are any files within the previously specified directory, if no files are detected it will exit
    # If there are files detected, all files' names will be saved in a list and then returned
    files = os.listdir(path)
    if (len(files) == 0):
        sys.exit("There are no files in the provided directory")
    print(to_color("[+]",'green'), "Files found!")
    return files


def seperate_marker_files(files):
    # This function splits marker files (files ending in "intervalMarkers.csv") and actual data files (ending in ".bp.csv")
    data_files = []
    marker_files = []
    for file in files:
        if file.lower().endswith('intervalmarker.csv'):
            marker_files.append(file)
        elif not file.lower().endswith('.json'):
            data_files.append(file)
    print(to_color("[+]",'green'), "Files are separated successfully")
    return marker_files, data_files


def processing(marker_files,data_files,path):
    # Set the pandas library to expand the nubmers instead of showing them like "1.725716e+09", this will also add extra precision and it will show
    # up to six digits after the dot. This will be needed for later processing
    pd.set_option("display.float_format", "{:.6f}".format)

    # Create output directory (if not created)
    output_path = create_output_directory(path)

    for file in marker_files:
        # Get the file_id which is composed as the following: SUBJECT-WORD-SESSION, and decompose the file_id to three different vairables, and these are
        # subject (has the subject name), word (the imagined and spoken word), and session (session ID)
        file_id = file.split("_")[1]
        subject_as_is, word_as_is, session_as_is = file_id.split("-")
        subject_lowercase = subject_as_is.lower()
        word_lowercase = word_as_is.lower()
        session_lowercase = session_as_is.lower()
        print(to_color("[+]","green"), "Currently working on: " + file_id)
        
        # Open a marker file in a dataframe, and then save the whole "timestamps" column in a variable
        df_marker = pd.read_csv(path+"\\"+file)
        timestamps = df_marker["timestamp"]
        print(to_color("[+]","green"), "Timestamps successfuly extracted from: " + file)

        # Find the actual data file pair of the currently markers file being processed
        data_file_pair = find_pair(file_id,data_files)

        # After identifying the pair of the current interval marker file, we now load the data file to start extracting the EEG data
        # Note how header = 1 is used to skip the informational header that is in each
        # The "k" is used for labeling the trial number within the session.
        # Getting iloc[0,0] of the dataframe is for output time labeling
        df_data = pd.read_csv(path+"\\"+data_file_pair, header=1)
        record_time_epoch_timestamp = df_data.iloc[0,0]
        print(to_color("[+]","green"), "Epoch time extracted: " + str(record_time_epoch_timestamp))
        k = 0

        # Convert epoch timestamp to GMT
        gmt_time = convert_epoch_to_gmt(record_time_epoch_timestamp)
        print(to_color("[+]","green"), "Epoch time converted to GMT time: " + gmt_time)

        # Now for the main algorithm that will separate recordings and save them to the output directory. This is will explained in the report, and it's also
        # relatively easy to track. Basically all the algorithm work relies on utilizing the indices (i,j,k) and the +2 increment each loop.
        for i in range(0, NUMBER_OF_MARKERS_PER_SESSION - 1, 2):
            j = i + 1
            start_value = timestamps.iloc[i]
            finish_value = timestamps.iloc[j]
            # TEST PURPOSES: print(OUTPUT_DIRECTORY_NAME)
            # TEST PURPOSES: print(start_value, finish_value)
            extracted_data = df_data[(df_data["Timestamp"] >= start_value) & (df_data["Timestamp"] <= finish_value)]


            # Now the data is extracted, it's time to throw the output file and name the session accordingly, we start by checking to see whether "i" is
            # multiple of 4 or not, if yes, then the extracted data is for a spoken signal. Otherwise it's an imagined one.
            if (is_multiple(i,4) == True):
                k = (k+1) % 6
                extracted_data.to_csv(
                    f"{output_path}\\S{subject[subject_lowercase]}_A{age[subject_lowercase]}_W{word[word_lowercase]}_C{data_class['spoken']}_D{gmt_time}_SID{session_lowercase}_T{k}.csv"
                    , index = False)
                print(to_color("[+]","green"), "File saved: " + f"{output_path}\\S{subject[subject_lowercase]}_A{age[subject_lowercase]}_W{word[word_lowercase]}_C{data_class['spoken']}_D{gmt_time}_SID{session_lowercase}_T{k}.csv")
            else:
                extracted_data.to_csv(
                    f"{output_path}\\S{subject[subject_lowercase]}_A{age[subject_lowercase]}_W{word[word_lowercase]}_C{data_class['imagined']}_D{gmt_time}_SID{session_lowercase}_T{k}.csv"
                    , index = False
                )
                print(to_color("[+]","green"), "File saved: " + f"{output_path}\\S{subject[subject_lowercase]}_A{age[subject_lowercase]}_W{word[word_lowercase]}_C{data_class['imagined']}_D{gmt_time}_SID{session_lowercase}_T{k}.csv")
        print(to_color("[+]","green"), "Finished working on: " + file_id)
        print()

            

        
def create_output_directory(path):
    output_directory = path+OUTPUT_DIRECTORY_NAME
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        print(to_color("[+]","green"), "Output files created successfuly at: " + output_directory)
        print()
    else:
        print(to_color("[-]",'yellow'),"Output directory already exists at " + output_directory)
        print()
    return output_directory


def convert_epoch_to_gmt(recod_time_epoch_timestamp):
    gmt_time = strftime("%d-%m-%Y-%H-%M", localtime(recod_time_epoch_timestamp))
    return gmt_time



def find_pair(file_id,data_files):
    # Gets the data file of a interval markers file
    for data_file in data_files:
        if file_id in data_file:
            print(to_color("[+]","green"), "Found pair file!: " + data_file)
            return data_file


def is_multiple(x,y):
    return x % y == 0

if __name__ == "__main__":
    main()