import csv
from pathlib import Path
import pandas as pd
import numpy
import statistics
import os

#################### README ##########################
# In order for this program to run, the following items must be in the same directory as this file
# 1 - a directory named "screening_colon" that contains colon patient data in .csv files
# 2 - a directory named "screening_ileum" that contains colon patient data in .csv files
# 3 - a .csv file named "seavue_anonymized_subject_endo_metadata" that contains label matching with patient IDs
# please replace the following value with the path to your directory that contains these items in the same format
path = "C://UROP"
# The final file will be named "final_clean_data.csv" in a directory named "CombinedData"


################### DO NOT EDIT #######################
your_path = Path(path)
# make directories
folder1 = "CleanData-ileum"
if not os.path.exists(folder1):
    # If it doesn't exist, create the folder
    os.makedirs(folder1)
folder2 = "CleanData-colon"
if not os.path.exists(folder2):
    # If it doesn't exist, create the folder
    os.makedirs(folder2)
folder3 = "CombinedData"
if not os.path.exists(folder3):
    # If it doesn't exist, create the folder
    os.makedirs(folder3)

# mode: i -> ileum, c -> colon
def clean(mode, directory_path, output_directory):


    # Keep track of skipped files
    skipped_files = []

    # Iterate through every file in the folder
    for file_path in directory_path.iterdir():
        # Check if it's a file (not a subdirectory)
        if file_path.is_file() and file_path.suffix == '.csv':

            print("Testing: ", file_path)
            fileGood = True

            if mode == 'i':
                # List of column indices you want to read from the CSV file (0-based index)
                column_indices_to_read = [1, 3, 4, 18, 19, 20, 21, 23, 24]
            elif mode == 'c':
                # List of column indices you want to read from the CSV file (0-based index)
                column_indices_to_read = [1, 3, 4, 19, 20, 21, 22, 24, 25]


            # Read specific columns from the CSV file into a DataFrame
            try:
                df = pd.read_csv(file_path, usecols=column_indices_to_read)
            except ValueError as e:
                print(f"Error: {e}")
                print("Skipping: ", file_path)
                skipped_files.append(file_path)
                continue

            # Create a Frame class
            class Frame:
                def __init__(self, row):
                    self.quality = row[0]
                    self.severity = row[1]
                    self.coord = row[2]
                    self.blood = row[3]
                    self.ulcerMild = row[4]
                    self.ulcerModSev = row[5]
                    self.erythma = row[6]
                    self.normal = row[7]
                    self.poorQuality = row[8]

            # Create a list to store frames
            frames = []

            # Iterate over each row in the DataFrame
            for index, row in df.iterrows():
                rowc = [row.iloc[0], row.iloc[1], row.iloc[2], row.iloc[3], row.iloc[4], row.iloc[5], row.iloc[6], 
                        row.iloc[7], row.iloc[8]]
                frame = Frame(rowc)
                frames.append(frame)

            ################# clean low-quality frames ###########################

            # Sort based on quality (low to high)
            sorted_frames = sorted(frames, key=lambda x: x.quality)

            # Fill a list with the total pixel numbers in each relevant frame (no concern)
            totals = []
            for frame in sorted_frames:
                if frame.blood == 0 and frame.ulcerMild == 0 and frame.ulcerModSev == 0 and frame.erythma == 0:
                    total = frame.normal + frame.poorQuality
                    totals.append(total)

            # account for saturated files
            if(totals == []) :
                for frame in sorted_frames:
                    if (frame.erythma == 0):
                        total = frame.normal + frame.poorQuality + frame.blood + frame.ulcerMild + frame.ulcerModSev
                        totals.append(total)
                    elif frame.ulcerMild == 0 and frame.ulcerModSev == 0 and frame.blood == 0:
                        total = frame.normal + frame.poorQuality + frame.erythma
                        totals.append(total)
            
            # Calculate the average number of total pixels
            average = sum(totals) / len(totals)

            # Find average quality
            q = df.iloc[:,0].mean()

            # Determine cutoff - needs tuning at end, for now we assume that the average video has a mean quality of .5, bell curve
            # quality = % to cut
            # .1 = .5
            # .9 = .1
            cut = (-0.5 * q) + 0.55

            # Calculate the number of items to remove 
            num_items_to_remove = int(len(sorted_frames) * cut)

            # Remove the items using slicing
            good_quality_frames = sorted_frames[num_items_to_remove:]

            ################### limit to one frame per coordinate ##################

            #################### method 1 - take medians ###############

            # Dictionary to store coordinates as keys and corresponding info as values
            dict = {}

            # Iterate through the list and fill dictionary
            for frame in good_quality_frames:
                if frame.coord not in dict:
                    dict[frame.coord] = {'quality': [frame.quality], 'severity': [frame.severity], 
                                        'blood': [frame.blood], 'ulcerMild': [frame.ulcerMild], 'ulcerModSev': [frame.ulcerModSev], 
                                        'erythma': [frame.erythma], 'normal': [frame.normal], 'poorQuality': [frame.poorQuality]}
                else:
                    dict[frame.coord]['quality'].append(frame.quality)
                    dict[frame.coord]['severity'].append(frame.severity)
                    dict[frame.coord]['blood'].append(frame.blood)
                    dict[frame.coord]['ulcerMild'].append(frame.ulcerMild)
                    dict[frame.coord]['ulcerModSev'].append(frame.ulcerModSev)
                    dict[frame.coord]['erythma'].append(frame.erythma)
                    dict[frame.coord]['normal'].append(frame.normal)
                    dict[frame.coord]['poorQuality'].append(frame.poorQuality)

            # Create new objects with unique coords and their corresponding median values
            unique_frames = []
            for coord, info in dict.items():
                median_quality = statistics.median(info['quality'])
                median_severity = statistics.median(info['severity'])
                median_blood = statistics.median(info['blood'])
                median_ulcerMild = statistics.median(info['ulcerMild'])
                median_ulcerModSev = statistics.median(info['ulcerModSev'])
                median_erythma = statistics.median(info['erythma'])
                median_normal = statistics.median(info['normal'])
                median_poorQuality = statistics.median(info['poorQuality'])
                data = [median_quality, median_severity, coord, median_blood, median_ulcerMild, median_ulcerModSev, 
                        median_erythma, median_normal, median_poorQuality]
                unique_frame = Frame(data)
                unique_frames.append(unique_frame)

            ################### method 2 - take highest quality frame for each coord
            '''

            # Dictionary to store coordinates as keys and corresponding info as values
            dict = {}

            # Iterate through the list and fill dictionary
            for frame in good_quality_frames:
                if frame.coord not in dict:
                    dict[frame.coord] = {'quality': frame.quality, 'severity': frame.severity, 
                                        'blood': frame.blood, 'ulcerMild': frame.ulcerMild, 'ulcerModSev': frame.ulcerModSev, 
                                        'erythma': frame.erythma, 'normal': frame.normal, 'poorQuality': frame.poorQuality}
                else:
                    if frame.quality > dict[frame.coord]['quality'] :
                        dict[frame.coord]['quality'] = frame.quality
                        dict[frame.coord]['severity'] = frame.severity
                        dict[frame.coord]['blood'] = frame.blood
                        dict[frame.coord]['ulcerMild'] = frame.ulcerMild
                        dict[frame.coord]['ulcerModSev'] = frame.ulcerModSev
                        dict[frame.coord]['erythma'] = frame.erythma
                        dict[frame.coord]['normal'] = frame.normal
                        dict[frame.coord]['poorQuality'] = frame.poorQuality

            # Create new objects with unique coords and their corresponding median values
            unique_frames = []
            for coord, info in dict.items():
                data = [info['quality'], info['severity'], coord, info['blood'], info['ulcerMild'], info['ulcerModSev'], 
                        info['erythma'], info['normal'], info['poorQuality']]
                unique_frame = Frame(data)
                unique_frames.append(unique_frame)

            '''

            ################## normalize pixel number ############################
            ############ everything will now switch in units   pixel # -> % of pixels

            # Convert units pixel# -> % pixels
            for frame in unique_frames:
                frame.blood = frame.blood / average
                frame.ulcerMild = frame.ulcerMild / average
                frame.ulcerModSev = frame.ulcerModSev / average
                frame.erythma = frame.erythma / average
                frame.normal = frame.normal / average
                frame.poorQuality = frame.poorQuality / average

            ################### divide into 5 or 200 sections ############################

            if mode == 'i':
                num_sec = 5
            elif mode == 'c':
                num_sec = 200

            # Sort by coord (low to high)
            sorted_byCoord_frames = sorted(unique_frames, key=lambda x: x.coord)

            #  Calculate interval length
            start_coord = sorted_byCoord_frames[0].coord
            last_coord = sorted_byCoord_frames[-1].coord
            range1 = last_coord - start_coord
            interval = range1 / num_sec

            # Separate into intervals
            condensed_frames = []
            for i in range(num_sec):

                # Calculate coord threshold
                threshold = ((i+1) * interval) + start_coord

                # Make dictionary to calculate medians
                temp_dict = {
                    'quality': [],
                    'severity': [],
                    'blood': [],
                    'ulcerMild': [],
                    'ulcerModSev': [],
                    'erythma': [],
                    'normal': [],
                    'poorQuality': []
                }

                # Fill dictionary using frames in this section
                for frame in sorted_byCoord_frames:
                    if frame.coord < threshold:
                        temp_dict['quality'].append(frame.quality)
                        temp_dict['severity'].append(frame.severity)
                        temp_dict['blood'].append(frame.blood)
                        temp_dict['ulcerMild'].append(frame.ulcerMild)
                        temp_dict['ulcerModSev'].append(frame.ulcerModSev)
                        temp_dict['erythma'].append(frame.erythma)
                        temp_dict['normal'].append(frame.normal)
                        temp_dict['poorQuality'].append(frame.poorQuality)
                        sorted_byCoord_frames.pop(0)
                    else:
                        break

                # Calculate medians for each feature
                try:
                    median_quality = statistics.median(temp_dict['quality'])
                    median_severity = statistics.median(temp_dict['severity'])
                    median_blood = statistics.median(temp_dict['blood'])
                    median_ulcerMild = statistics.median(temp_dict['ulcerMild'])
                    median_ulcerModSev = statistics.median(temp_dict['ulcerModSev'])
                    median_erythma = statistics.median(temp_dict['erythma'])
                    median_normal = statistics.median(temp_dict['normal'])
                    median_poorQuality = statistics.median(temp_dict['poorQuality'])
                except statistics.StatisticsError as e:
                    # Account for the sections with no frames
                    print(f"Error: {e}")
                    print("Skipping: ", file_path)
                    skipped_files.append(file_path)
                    fileGood = False
                    break

                # Make new frame for the whole section and add to list
                data = [median_quality, median_severity, (i+1), median_blood, median_ulcerMild, median_ulcerModSev, 
                        median_erythma, median_normal, median_poorQuality]
                condensed_frame = Frame(data)
                condensed_frames.append(condensed_frame)

            # check to skip file
            if ( not fileGood):
                continue

            #################### condense again ##################################

            # create dictionary for sections
            temp_dict = {
                    'quality': [],
                    'severity': [],
                    'blood': [],
                    'ulcerMild': [],
                    'ulcerModSev': [],
                    'erythma': [],
                    'normal': [],
                    'poorQuality': []
                }
            

            twice_condensed_frames = []
            count = 0
            for frame in condensed_frames:
                count += 1

                # fill dictionary
                temp_dict['quality'].append(frame.quality)
                temp_dict['severity'].append(frame.severity)
                temp_dict['blood'].append(frame.blood)
                temp_dict['ulcerMild'].append(frame.ulcerMild)
                temp_dict['ulcerModSev'].append(frame.ulcerModSev)
                temp_dict['erythma'].append(frame.erythma)
                temp_dict['normal'].append(frame.normal)
                temp_dict['poorQuality'].append(frame.poorQuality)
            
                # condense
                if mode == 'i':
                    cnum = 1
                elif mode == 'c':
                    cnum = 10
                if (count % cnum == 0):
                    # Calculate medians for each feature
                    try:
                        median_quality = statistics.median(temp_dict['quality'])
                        median_severity = statistics.median(temp_dict['severity'])
                        median_blood = statistics.median(temp_dict['blood'])
                        median_ulcerMild = statistics.median(temp_dict['ulcerMild'])
                        median_ulcerModSev = statistics.median(temp_dict['ulcerModSev'])
                        median_erythma = statistics.median(temp_dict['erythma'])
                        median_normal = statistics.median(temp_dict['normal'])
                        median_poorQuality = statistics.median(temp_dict['poorQuality'])
                    except (numpy.core._exceptions._UFuncNoLoopError, TypeError) as e:
                        # Account for the sections with no frames
                        print(f"Error: {e}")
                        print("Skipping: ", file_path)
                        skipped_files.append(file_path)
                        fileGood = False
                        break
                        
                    # Make new frame for the whole section and add to list
                    data = [median_quality, median_severity, 0, median_blood, median_ulcerMild, median_ulcerModSev, 
                        median_erythma, median_normal, median_poorQuality]
                    twice_condensed_frame = Frame(data)
                    twice_condensed_frames.append(twice_condensed_frame)

                    # Reset
                    temp_dict['quality'] = []
                    temp_dict['severity'] = []
                    temp_dict['blood'] = []
                    temp_dict['ulcerMild'] = []
                    temp_dict['ulcerModSev'] = []
                    temp_dict['erythma'] = []
                    temp_dict['normal'] = []
                    temp_dict['poorQuality'] = []

            # check to skip file
            if ( not fileGood):
                continue

            #################### calculate stats and fill final file ###############
            
            # Specify the CSV file path
            csv_file_path = output_directory / (file_path.stem + '_clean.csv')

            # Write the list of objects to the CSV file
            with open(csv_file_path, mode='w', newline='') as file:

                writer = csv.writer(file)

                # Write the header row
                writer.writerow(["patient_id", "section", "quality1", "severity1", "blood1", "ulcer_mild1", "ulcer_moderate_severe1", 
                                 "erythma1", "normal1", "poor_quality1", "quality2", "severity2", "blood2", "ulcer_mild2", "ulcer_moderate_severe2", 
                                 "erythma2", "normal2", "poor_quality2", "quality3", "severity3", "blood3", "ulcer_mild3", "ulcer_moderate_severe3", 
                                 "erythma3", "normal3", "poor_qualit3y", "quality4", "severity4", "blood4", "ulcer_mild4", "ulcer_moderate_severe4", 
                                 "erythma4", "normal4", "poor_quality4", "quality5", "severity5", "blood5", "ulcer_mild5", "ulcer_moderate_severe5", 
                                 "erythma5", "normal5", "poor_quality5"])

                # Write each frame as a row in the CSV file
                count = 0
                patient_id = os.path.basename(file_path)
                patient_id = str(patient_id)[8:14]
                print(patient_id)
                if mode == 'i':
                    data_to_write = [patient_id, "ileum"]
                    for frame in twice_condensed_frames:
                        data_to_write.append(frame.quality)
                        data_to_write.append(frame.severity)
                        data_to_write.append(frame.blood)
                        data_to_write.append(frame.ulcerMild)
                        data_to_write.append(frame.ulcerModSev)
                        data_to_write.append(frame.erythma)
                        data_to_write.append(frame.normal)
                        data_to_write.append(frame.poorQuality)
                    writer.writerow(data_to_write)
                elif mode == 'c':
                    count = 0
                    data_to_write = [patient_id, "left"]
                    for frame in twice_condensed_frames:
                        count += 1
                        data_to_write.append(frame.quality)
                        data_to_write.append(frame.severity)
                        data_to_write.append(frame.blood)
                        data_to_write.append(frame.ulcerMild)
                        data_to_write.append(frame.ulcerModSev)
                        data_to_write.append(frame.erythma)
                        data_to_write.append(frame.normal)
                        data_to_write.append(frame.poorQuality)

                        if (count % 5 == 0):
                            writer.writerow(data_to_write)
                            if (count == 5):
                                data_to_write = [patient_id, "top"]
                            elif (count == 10):
                                data_to_write = [patient_id, "right"]
                            else:
                                data_to_write = [patient_id, "rectum"]


################# make clean data files #############################

# Specify the directory path you want to iterate through - ileum
directory_path = your_path / 'screening_ileum'

# Specify the output directory - ileum
output_directory = your_path / "CleanData-ileum"

# Clean ileum data
clean('i', directory_path, output_directory)

# Specify the directory path you want to iterate through - colon
directory_path = your_path / "screening_colon"

# Specify the output directory - colon
output_directory = your_path / "CleanData-colon"

# Clean colon data
clean('c', directory_path, output_directory)

###################### condense to one file ##############################

# Update output directory and the name of output file
output_directory = your_path / "CombinedData"
final_file_path = output_directory / "CombinedData.csv"

with open(final_file_path, mode='w', newline='') as file:
    writer = csv.writer(file)

    # Write the header row
    writer.writerow(["patient_id", "section", "quality1", "severity1", "blood1", "ulcer_mild1", "ulcer_moderate_severe1", 
                    "erythma1", "normal1", "poor_quality1", "quality2", "severity2", "blood2", "ulcer_mild2", "ulcer_moderate_severe2", 
                    "erythma2", "normal2", "poor_quality2", "quality3", "severity3", "blood3", "ulcer_mild3", "ulcer_moderate_severe3", 
                    "erythma3", "normal3", "poor_qualit3y", "quality4", "severity4", "blood4", "ulcer_mild4", "ulcer_moderate_severe4", 
                    "erythma4", "normal4", "poor_quality4", "quality5", "severity5", "blood5", "ulcer_mild5", "ulcer_moderate_severe5", 
                    "erythma5", "normal5", "poor_quality5"])

    # First directory
    directory_path = your_path / "CleanData-colon"

    # Iterate through every file in the colon folder
    for file_path in directory_path.iterdir():
        # Check if it's a file (not a subdirectory)
        if file_path.is_file() and file_path.suffix == '.csv':

            print("combining: ", file_path)

            # Open the CSV file and read each row of data
            with open(file_path, newline='') as csvfile:
                csv_reader = csv.reader(csvfile)
                count = 0
                for row in csv_reader:
                    # Skip header
                    if (count != 0):
                        writer.writerow(row)
                    count += 1

    # Second directory
    directory_path = your_path / "CleanData-ileum"

    # Iterate through every file in the ileum folder
    for file_path in directory_path.iterdir():
        # Check if it's a file (not a subdirectory)
        if file_path.is_file() and file_path.suffix == '.csv':

            print("combining: ", file_path)

            # Open the CSV file and read each row of data
            with open(file_path, newline='') as csvfile:
                csv_reader = csv.reader(csvfile)
                count = 0
                for row in csv_reader:
                    # Skip header
                    if (count != 0):
                        writer.writerow(row)
                    count += 1



#################### make ulcer and erythma size features #####################

# read combined data into data frame
file_path = your_path / "CombinedData//CombinedData.csv"
df = pd.read_csv(file_path)

# Add ulcer and erythmas for each section
df['ulcer_size1'] = df.apply(lambda row: row['ulcer_mild1'] + row['ulcer_moderate_severe1'], axis=1)
df['ulcer_size2'] = df.apply(lambda row: row['ulcer_mild2'] + row['ulcer_moderate_severe2'], axis=1)
df['ulcer_size3'] = df.apply(lambda row: row['ulcer_mild3'] + row['ulcer_moderate_severe3'], axis=1)
df['ulcer_size4'] = df.apply(lambda row: row['ulcer_mild4'] + row['ulcer_moderate_severe4'], axis=1)
df['ulcer_size5'] = df.apply(lambda row: row['ulcer_mild5'] + row['ulcer_moderate_severe5'], axis=1)
df['erythma_size1'] = df['erythma1']
df['erythma_size2'] = df['erythma2']
df['erythma_size3'] = df['erythma3']
df['erythma_size4'] = df['erythma4']
df['erythma_size5'] = df['erythma5']

# calculate percentiles
ulcer_values = df[['ulcer_size1', 'ulcer_size2', 'ulcer_size3', 'ulcer_size4', 'ulcer_size5']].values.flatten()
ulcer_p75 = numpy.percentile(ulcer_values, 75)
ulcer_p50 = numpy.percentile(ulcer_values, 50)
ulcer_p25 = numpy.percentile(ulcer_values, 25)
erythma_values = df[['erythma_size1', 'erythma_size2', 'erythma_size3', 'erythma_size4', 'erythma_size5']].values.flatten()
erythma_p75 = numpy.percentile(erythma_values, 75)
erythma_p50 = numpy.percentile(erythma_values, 50)
erythma_p25 = numpy.percentile(erythma_values, 25)

# function that gives a score based on percentile location
def categorize_value(value, p25, p50, p75):
    if value <= p25:
        return 0
    elif value < p50:
        return 1
    elif value < p75:
        return 2
    else:
        return 3
    
# replace values
df['ulcer_size1'] = df['ulcer_size1'].apply(lambda x: categorize_value(x, ulcer_p25, ulcer_p50, ulcer_p75))
df['ulcer_size2'] = df['ulcer_size2'].apply(lambda x: categorize_value(x, ulcer_p25, ulcer_p50, ulcer_p75))
df['ulcer_size3'] = df['ulcer_size3'].apply(lambda x: categorize_value(x, ulcer_p25, ulcer_p50, ulcer_p75))
df['ulcer_size4'] = df['ulcer_size4'].apply(lambda x: categorize_value(x, ulcer_p25, ulcer_p50, ulcer_p75))
df['ulcer_size5'] = df['ulcer_size5'].apply(lambda x: categorize_value(x, ulcer_p25, ulcer_p50, ulcer_p75))
df['erythma_size1'] = df['erythma_size1'].apply(lambda x: categorize_value(x, erythma_p25, erythma_p50, erythma_p75))
df['erythma_size2'] = df['erythma_size2'].apply(lambda x: categorize_value(x, erythma_p25, erythma_p50, erythma_p75))
df['erythma_size3'] = df['erythma_size3'].apply(lambda x: categorize_value(x, erythma_p25, erythma_p50, erythma_p75))
df['erythma_size4'] = df['erythma_size4'].apply(lambda x: categorize_value(x, erythma_p25, erythma_p50, erythma_p75))
df['erythma_size5'] = df['erythma_size5'].apply(lambda x: categorize_value(x, erythma_p25, erythma_p50, erythma_p75))

##################### make total calculations for each of the sections 

# make new columns for each section
df['total_ulcer_average'] = df.apply(lambda row: (row['ulcer_mild1'] + row['ulcer_moderate_severe1'] + row['ulcer_mild2'] + row['ulcer_moderate_severe2'] + row['ulcer_mild3'] + row['ulcer_moderate_severe3'] + row['ulcer_mild4'] + row['ulcer_moderate_severe4'] + row['ulcer_mild5'] + row['ulcer_moderate_severe5']) / 5, axis=1)
df['total_erythma_average'] = df.apply(lambda row: (row['erythma1'] + row['erythma2'] + row['erythma3'] + row['erythma4'] + row['erythma5']) / 5, axis=1)

# calculate percentiles
total_ulcer_values = df['total_ulcer_average'].values.flatten()
tulcer_p75 = numpy.percentile(total_ulcer_values, 75)
tulcer_p50 = numpy.percentile(total_ulcer_values, 50)
tulcer_p25 = numpy.percentile(total_ulcer_values, 25)
total_erythma_values = df['total_erythma_average'].values.flatten()
terythma_p75 = numpy.percentile(total_erythma_values, 75)
terythma_p50 = numpy.percentile(total_erythma_values, 50)
terythma_p25 = numpy.percentile(total_erythma_values, 25)

# replace
df['total_ulcer_average'] = df['total_ulcer_average'].apply(lambda x: categorize_value(x, tulcer_p25, tulcer_p50, tulcer_p75))
df['total_erythma_average'] = df['total_erythma_average'].apply(lambda x: categorize_value(x, terythma_p25, terythma_p50, terythma_p75))

#################### match labels ##########################

# read in labels to new data frame 
selected_columns = ['SubjectID', 'VISIT', 'SESCDI', 'SESCDRC', 'SESCDTC', 'SESCDLC', 'SESCDRE']
labels_file_path = your_path / "seavue_anonymized_subject_endo_metadata.csv"
df_labels = pd.read_csv(labels_file_path, usecols=selected_columns)
df_labels.drop(df_labels[df_labels['VISIT'] != "Screening"].index, inplace=True)
df_labels.rename(columns={'SESCDI': 'ileum', 'SESCDRC': 'right', 'SESCDTC': 'top', 'SESCDLC': 'left', 'SESCDRE': 'rectum'}, inplace=True)
print(df_labels)

# match labels to patient
df_labels.set_index('SubjectID', inplace=True)
# Convert 'SubjectID' column in df to int
df_labels.index = df_labels.index.astype(int)
df['patient_id'] = df['patient_id'].astype(int)
def lookup_label(row):
    return df_labels.loc[row['patient_id'], row['section']]
df['label'] = df.apply(lookup_label, axis=1)

# Remove rows with blank cells in the 'label' column
rows_with_blank_label = df['label'].isna() 
df = df[~rows_with_blank_label]

# switch column order
new_column_order = ["patient_id", "section", "label",  "total_ulcer_average", "total_erythma_average", "quality1", "severity1", "blood1", "ulcer_mild1", "ulcer_moderate_severe1", "ulcer_size1",
                    "erythma1", "erythma_size1", "normal1", "poor_quality1", "quality2", "severity2", "blood2", "ulcer_mild2", "ulcer_moderate_severe2", "ulcer_size2",
                    "erythma2", "erythma_size2", "normal2", "poor_quality2", "quality3", "severity3", "blood3", "ulcer_mild3", "ulcer_moderate_severe3", "ulcer_size3",
                    "erythma3", "erythma_size3", "normal3", "poor_qualit3y", "quality4", "severity4", "blood4", "ulcer_mild4", "ulcer_moderate_severe4", "ulcer_size4",
                    "erythma4", "erythma_size4", "normal4", "poor_quality4", "quality5", "severity5", "blood5", "ulcer_mild5", "ulcer_moderate_severe5", "ulcer_size5",
                    "erythma5", "erythma_size5", "normal5", "poor_quality5"]

# write to final file
df = df[new_column_order]
df.to_csv(your_path / "CombinedData//final_clean_data.csv", index=False)



