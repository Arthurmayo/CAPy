import os

def convert_awd_to_csv(awd_file_path, csv_file_path):
    """
    Convert a .awd file to a cleaned .csv file for use with CAPy.
    
    This converter does the following:
      - Preserves the first 7 lines of meta data but ensures each of these lines has two fields.
          * For the first 4 lines, if only one field is present, a blank second field is added.
          * For lines 5–7, which are to be cleared, a blank row (with two empty fields) is written.
      - For every subsequent line (starting at line 8):
          * Splits the line by whitespace.
          * Parses the first token as the event count.
          * If the event count is negative, the entire line is skipped.
          * Otherwise, if a second token exists, uses it as the light value; if missing, uses 0.
          * Writes the row as two comma-separated values: event_count,light_value.
    
    Parameters:
        awd_file_path (str): Path to the input .awd file.
        csv_file_path (str): Path to the output .csv file.
    
    Returns:
        str: The output CSV file path if conversion is successful.
    
    Raises:
        IOError: If reading or writing the file fails.
    """
    try:
        with open(awd_file_path, 'r') as infile:
            lines = infile.readlines()
    except Exception as e:
        raise IOError(f"Error reading AWD file: {e}")
    
    new_lines = []
    total_lines = len(lines)
    for i, line in enumerate(lines):
        
        if i < 4:
            if ',' in line:
                new_lines.append(line)
            else:
                new_lines.append(line.rstrip('\n') + ",\n")
        elif 4 <= i < 7:
            new_lines.append(",\n")
        else:

            parts = line.strip().split()
            if not parts:
                new_lines.append(",\n")
                continue
            try:
                event_val = float(parts[0])
            except ValueError:
                continue

            if event_val < 0:
                continue

            if len(parts) >= 2:
                try:
                    light_val = float(parts[1])
                except ValueError:
                    light_val = 0
            else:
                light_val = 0
            new_lines.append(f"{event_val},{light_val}\n")
    
    try:
        with open(csv_file_path, 'w') as outfile:
            outfile.writelines(new_lines)
    except Exception as e:
        raise IOError(f"Error writing CSV file: {e}")
    
    return csv_file_path

if __name__ == "__main__":

    import sys
    if len(sys.argv) != 3:
        print("Usage: python awd_converter.py input.awd output.csv")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        try:
            convert_awd_to_csv(input_file, output_file)
            print(f"Conversion successful: {output_file}")
        except Exception as e:
            print(f"Conversion failed: {e}")
